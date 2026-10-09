"use strict";

const vscode = require("vscode");
const crypto = require("node:crypto");
const { buildGroups, wrapIndex } = require("./model");

let currentPanel;

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (char) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  })[char]);
}

async function isDirectory(uri) {
  try {
    return Boolean((await vscode.workspace.fs.stat(uri)).type & vscode.FileType.Directory);
  } catch {
    return false;
  }
}

async function findDataRoot() {
  for (const folder of vscode.workspace.workspaceFolders || []) {
    for (const parts of [
      ["RRS_XCO2", "bottom_layer_XCO2_retrievals"],
      ["bottom_layer_XCO2_retrievals"],
    ]) {
      const root = vscode.Uri.joinPath(folder.uri, ...parts);
      if (await isDirectory(root)) return root;
    }
    if (folder.uri.path.endsWith("/bottom_layer_XCO2_retrievals")) return folder.uri;
  }
  const selection = await vscode.window.showOpenDialog({
    canSelectFiles: false, canSelectFolders: true, canSelectMany: false,
    openLabel: "Use retrieval data",
    title: "Select the bottom_layer_XCO2_retrievals folder",
  });
  return selection && selection[0];
}

function html(webview, extensionUri) {
  const nonce = crypto.randomBytes(24).toString("base64");
  const uri = (...parts) => escapeHtml(webview.asWebviewUri(
    vscode.Uri.joinPath(extensionUri, ...parts)).toString());
  const csp = escapeHtml(webview.cspSource);
  return [
    '<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; img-src ' +
      csp + '; style-src ' + csp + '; script-src \'nonce-' + nonce + '\';">',
    '<link rel="stylesheet" href="' + uri("media", "viewer.css") + '">',
    '<title>RRS Ensemble Comparison</title></head><body tabindex="-1">',
    '<header><div class="toolbar">',
    '<button id="previous" title="Previous group (Left arrow)" aria-label="Previous group">←</button>',
    '<select id="group" aria-label="Select four-state group"></select>',
    '<button id="next" title="Next group (Right arrow)" aria-label="Next group">→</button>',
    '<button id="refresh" title="Recheck plot files">Refresh</button>',
    '<button id="restore" hidden>Show all eight</button><span id="position"></span></div>',
    '<div class="details"><span id="states">Loading plot inventory…</span>',
    '<span>← / →: four-state groups · click image: enlarge · Esc: restore grid</span></div>',
    '<div id="status" role="status" aria-live="polite"></div></header>',
    '<main id="grid" aria-label="Top: aerosol and surface, rounds 3–6. Bottom: CO2 and SIF, rounds 3–6."></main>',
    '<script nonce="' + nonce + '" src="' + uri("model.js") + '"></script>',
    '<script nonce="' + nonce + '" src="' + uri("media", "viewer.js") + '"></script>',
    '</body></html>',
  ].join("\n");
}

function activate(context) {
  const open = async () => {
    if (currentPanel) {
      currentPanel.reveal();
      return;
    }
    const root = await findDataRoot();
    if (!root) return;
    const panel = vscode.window.createWebviewPanel(
      "rrsEnsembleViewer", "RRS Ensembles · Rounds 3–6",
      vscode.ViewColumn.Active,
      {
        enableScripts: true,
        localResourceRoots: [root, context.extensionUri],
        enableFindWidget: false,
      });
    currentPanel = panel;
    let disposed = false;
    let inventoryVersion = 0;
    const stateKey = "rrsEnsemble.group:" + root.toString();

    const refresh = async () => {
      const version = ++inventoryVersion;
      try {
        const groups = buildGroups();
        await Promise.all(groups.flatMap((group) => group.panes.map(async (pane) => {
          const file = vscode.Uri.joinPath(root, ...pane.path.split("/"));
          try {
            const info = await vscode.workspace.fs.stat(file);
            pane.exists = Boolean(info.type & vscode.FileType.File);
            pane.url = panel.webview.asWebviewUri(file).with({
              query: "v=" + info.mtime + "-" + info.size,
            }).toString();
          } catch {
            pane.exists = false;
            pane.url = "";
          }
        })));
        if (disposed || version !== inventoryVersion) return;
        await panel.webview.postMessage({
          type: "inventory", groups,
          index: wrapIndex(context.workspaceState.get(stateKey, 0), groups.length),
        });
      } catch (error) {
        if (!disposed) {
          await panel.webview.postMessage({ type: "error", message: String(error.message) });
        }
      }
    };
    context.subscriptions.push(panel.onDidDispose(() => {
      disposed = true;
      if (currentPanel === panel) currentPanel = undefined;
    }));
    context.subscriptions.push(panel.webview.onDidReceiveMessage(async (message) => {
      if (!message || typeof message !== "object") return;
      if (message.type === "ready" || message.type === "refresh") await refresh();
      if (message.type === "selection" && Number.isInteger(message.index) &&
          message.index >= 0 && message.index < 20) {
        await context.workspaceState.update(stateKey, message.index);
      }
    }));
    panel.webview.html = html(panel.webview, context.extensionUri);
  };

  context.subscriptions.push(
    vscode.commands.registerCommand("rrsEnsemble.open", open),
    vscode.commands.registerCommand("rrsEnsemble.refresh", async () => {
      if (!currentPanel) return open();
      currentPanel.reveal();
      await currentPanel.webview.postMessage({ type: "refreshRequested" });
    }),
  );
}

function deactivate() {
  if (currentPanel) currentPanel.dispose();
  currentPanel = undefined;
}

module.exports = { activate, deactivate };
