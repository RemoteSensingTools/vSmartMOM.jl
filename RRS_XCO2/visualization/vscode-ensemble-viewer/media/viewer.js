"use strict";

(() => {
  const vscode = acquireVsCodeApi();
  const { wrapIndex } = window.rrsModel;
  const grid = document.getElementById("grid");
  const picker = document.getElementById("group");
  const status = document.getElementById("status");
  const restore = document.getElementById("restore");
  let groups = [];
  let index = 0;
  let generation = 0;
  let enlarged = -1;
  const initial = vscode.getState();

  function applyEnlargement() {
    grid.classList.toggle("enlarged", enlarged >= 0);
    restore.hidden = enlarged < 0;
    Array.from(grid.children).forEach((pane, n) => {
      pane.classList.toggle("chosen", n === enlarged);
    });
  }

  function togglePane(n) {
    enlarged = enlarged === n ? -1 : n;
    applyEnlargement();
  }

  async function loadImage(pane) {
    if (!pane.exists) return null;
    const image = new Image();
    image.alt = pane.label;
    return new Promise((resolve) => {
      let settled = false;
      const finish = (value) => {
        if (settled) return;
        settled = true;
        clearTimeout(timer);
        image.onload = null;
        image.onerror = null;
        resolve(value);
      };
      const timer = setTimeout(() => finish(null), 20000);
      image.onload = () => finish(image);
      image.onerror = () => finish(null);
      image.src = pane.url;
    });
  }

  async function showGroup(next) {
    if (!groups.length) return;
    index = wrapIndex(next, groups.length);
    const requested = index;
    const ticket = ++generation;
    const group = groups[index];
    status.textContent = "Loading group " + (index + 1) + "/" + groups.length + "…";
    picker.value = String(index);
    // Commit all eight panes together; stale loads cannot replace a newer group.
    const images = await Promise.all(group.panes.map(loadImage));
    if (ticket !== generation) return;
    const fragment = document.createDocumentFragment();
    group.panes.forEach((pane, n) => {
      const section = document.createElement("section");
      section.className = "pane";
      const heading = document.createElement("h2");
      heading.textContent = pane.label;
      section.appendChild(heading);
      if (images[n]) {
        const button = document.createElement("button");
        button.className = "image-button";
        button.title = "Enlarge / restore " + pane.label;
        button.setAttribute("aria-label", button.title);
        button.appendChild(images[n]);
        button.addEventListener("click", () => togglePane(n));
        section.appendChild(button);
      } else {
        const missing = document.createElement("div");
        missing.className = "missing";
        const message = document.createElement("p");
        message.textContent = pane.exists ? "Image could not be loaded. Try Refresh." : "Plot not available";
        const path = document.createElement("code");
        path.textContent = pane.path;
        missing.append(message, path);
        section.appendChild(missing);
      }
      fragment.appendChild(section);
    });
    grid.replaceChildren(fragment);
    applyEnlargement();
    document.getElementById("position").textContent = "Group " + (index + 1) + " / " + groups.length;
    document.getElementById("states").textContent =
      group.label + " · states " + group.states.join(", ") + " (urban, rural, desert, forest)";
    const failures = images.filter((image) => !image).length;
    status.textContent = failures
      ? failures + " of 8 plots unavailable; group alignment preserved."
      : "All 8 plots loaded. At the last group, → wraps to the first.";
    vscode.setState({ index: requested });
    vscode.postMessage({ type: "selection", index: requested });
  }

  document.getElementById("previous").addEventListener("click", () => showGroup(index - 1));
  document.getElementById("next").addEventListener("click", () => showGroup(index + 1));
  restore.addEventListener("click", () => { enlarged = -1; applyEnlargement(); });
  document.getElementById("refresh").addEventListener("click", () => {
    status.textContent = "Refreshing plot inventory…";
    vscode.postMessage({ type: "refresh" });
  });
  picker.addEventListener("change", () => {
    showGroup(Number(picker.value));
    document.body.focus();
  });
  window.addEventListener("keydown", (event) => {
    if (event.altKey || event.ctrlKey || event.metaKey) return;
    if (event.target && event.target.closest("select, input, textarea, [contenteditable=true]")) return;
    if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
      event.preventDefault();
      showGroup(index + (event.key === "ArrowLeft" ? -1 : 1));
    } else if (event.key === "Escape") {
      enlarged = -1;
      applyEnlargement();
    }
  });
  let firstInventory = true;
  window.addEventListener("message", (event) => {
    const message = event.data;
    if (!message || typeof message !== "object") return;
    if (message.type === "error") {
      status.textContent = message.message;
    } else if (message.type === "refreshRequested") {
      vscode.postMessage({ type: "refresh" });
    } else if (message.type === "inventory") {
      groups = message.groups;
      picker.replaceChildren(...groups.map((group, n) => {
        const option = document.createElement("option");
        option.value = String(n);
        option.textContent = (n + 1) + ". " + group.label;
        return option;
      }));
      const next = firstInventory ? (initial ? initial.index : message.index) : index;
      firstInventory = false;
      showGroup(wrapIndex(next, groups.length));
      document.body.focus();
    }
  });
  vscode.postMessage({ type: "ready" });
})();
