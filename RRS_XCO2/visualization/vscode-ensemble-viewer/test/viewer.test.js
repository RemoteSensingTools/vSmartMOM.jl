"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const { buildGroups, wrapIndex } = require("../model");

class Element {
  constructor(tag) {
    this.tag = tag;
    this.children = [];
    this.listeners = {};
    this.textContent = "";
    this.classes = new Set();
    this.classList = {
      toggle: (name, enabled) => enabled ? this.classes.add(name) : this.classes.delete(name),
    };
  }
  appendChild(child) { this.children.push(child); return child; }
  append(...children) { this.children.push(...children); }
  replaceChildren(...children) {
    this.children = children.flatMap((child) => child.tag === "fragment" ? child.children : [child]);
  }
  setAttribute(key, value) { this[key] = value; }
  addEventListener(key, fn) { this.listeners[key] = fn; }
  closest() { return ["select", "input", "textarea"].includes(this.tag) ? this : null; }
  focus() { this.focused = true; }
}

function harness() {
  const elements = Object.fromEntries(
    ["grid", "group", "status", "restore", "previous", "next", "refresh", "states", "position"]
      .map((id) => [id, new Element(id === "group" ? "select" : "div")]));
  const body = new Element("body");
  const listeners = {};
  const pending = [];
  const messages = [];
  let saved;
  class FakeImage extends Element {
    constructor() { super("img"); }
    set src(value) { this.url = value; pending.push(this); }
  }
  const context = vm.createContext({
    window: { rrsModel: { wrapIndex }, addEventListener: (key, fn) => { listeners[key] = fn; } },
    document: {
      body,
      getElementById: (id) => elements[id],
      createElement: (tag) => new Element(tag),
      createDocumentFragment: () => new Element("fragment"),
    },
    Image: FakeImage,
    setTimeout: () => 1, clearTimeout: () => {},
    acquireVsCodeApi: () => ({
      getState: () => undefined,
      setState: (value) => { saved = value; },
      postMessage: (value) => { messages.push(value); },
    }),
  });
  vm.runInContext(fs.readFileSync(path.join(__dirname, "../media/viewer.js"), "utf8"), context);
  const groups = buildGroups();
  for (const group of groups) for (const pane of group.panes) {
    pane.exists = true;
    pane.url = "local:" + pane.path;
  }
  const key = (value, target = body) => {
    let prevented = false;
    listeners.keydown({ key: value, target, preventDefault: () => { prevented = true; } });
    return prevented;
  };
  const complete = async (images) => {
    for (const image of images) if (image.onload) image.onload();
    await new Promise((resolve) => setImmediate(resolve));
  };
  return {
    elements, groups, pending, messages, key, complete,
    selection: () => saved,
    inventory: () => listeners.message({ data: { type: "inventory", groups, index: 0 } }),
  };
}

test("renders eight panes in the required row and round order", async () => {
  const h = harness();
  assert.equal(h.messages[0].type, "ready");
  h.inventory();
  assert.equal(h.pending.length, 8);
  await h.complete(h.pending);
  assert.equal(h.elements.grid.children.length, 8);
  assert.deepEqual(h.elements.grid.children.map((p) => p.children[0].textContent),
    h.groups[0].panes.map((p) => p.label));
  assert.match(h.elements.states.textContent, /001, 021, 041, 061/);
  assert.match(h.elements.status.textContent, /All 8 plots loaded/);
});

test("arrow keys synchronize all panes and wrap; dropdown retains its own arrows", async () => {
  const h = harness();
  h.inventory();
  await h.complete(h.pending);
  assert.equal(h.key("ArrowLeft"), true);
  await h.complete(h.pending.slice(-8));
  assert.equal(h.selection().index, 19);
  assert.match(h.elements.states.textContent, /020, 040, 060, 080/);
  h.key("ArrowRight");
  await h.complete(h.pending.slice(-8));
  assert.equal(h.selection().index, 0);
  assert.equal(h.key("ArrowRight", h.elements.group), false);
  assert.equal(h.selection().index, 0);
});

test("rapid navigation cannot let a stale load replace the requested group", async () => {
  const h = harness();
  h.inventory();
  await h.complete(h.pending);
  h.key("ArrowRight");
  const older = h.pending.slice(-8);
  h.key("ArrowRight");
  const newer = h.pending.slice(-8);
  await h.complete(newer);
  await h.complete(older);
  assert.equal(h.selection().index, 2);
  assert.match(h.elements.states.textContent, /400 ppm/);
  for (const pane of h.elements.grid.children) {
    assert.match(pane.children[1].children[0].url, /bottom_co2_400_/);
  }
});

test("click enlarges a pane; Escape restores the complete grid", async () => {
  const h = harness();
  h.inventory();
  await h.complete(h.pending);
  h.elements.grid.children[5].children[1].listeners.click();
  assert.ok(h.elements.grid.classes.has("enlarged"));
  assert.ok(h.elements.grid.children[5].classes.has("chosen"));
  assert.equal(h.elements.restore.hidden, false);
  h.key("Escape");
  assert.ok(!h.elements.grid.classes.has("enlarged"));
  assert.equal(h.elements.restore.hidden, true);
});

test("missing plots preserve eight positions and identify the missing round", async () => {
  const h = harness();
  h.groups[0].panes[2].exists = false;
  h.inventory();
  assert.equal(h.pending.length, 7);
  await h.complete(h.pending);
  assert.equal(h.elements.grid.children.length, 8);
  const missing = h.elements.grid.children[2];
  assert.match(missing.children[0].textContent, /Round 5/);
  assert.equal(missing.children[1].children[0].textContent, "Plot not available");
  assert.match(h.elements.status.textContent, /1 of 8 plots unavailable/);
});
