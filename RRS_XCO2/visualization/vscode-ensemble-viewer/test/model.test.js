"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { buildGroups, plotRoot, wrapIndex } = require("../model");

test("20 groups cover all 80 truth states exactly once", () => {
  const groups = buildGroups();
  assert.equal(groups.length, 20);
  assert.deepEqual(groups[0].states, ["001", "021", "041", "061"]);
  assert.deepEqual(groups[19].states, ["020", "040", "060", "080"]);
  assert.deepEqual(groups.flatMap((g) => g.states).sort(),
    Array.from({ length: 80 }, (_, n) => String(n + 1).padStart(3, "0")));
  assert.deepEqual(groups.slice(0, 5).map((g) => g.ppm), [360, 380, 400, 420, 440]);
  assert.equal(groups[5].sif, true);
  assert.equal(groups[10].aerosol, "with_aerosol");
});

test("eight panes always have product rows and round columns", () => {
  const groups = buildGroups();
  for (const g of groups) {
    assert.equal(g.panes.length, 8);
    assert.deepEqual(g.panes.map((p) => p.round), [3, 4, 5, 6, 3, 4, 5, 6]);
    assert.deepEqual(g.panes.map((p) => p.product),
      Array(4).fill("profiles_surface").concat(Array(4).fill("co2_profiles")));
    for (const p of g.panes) {
      assert.ok(p.path.includes("_" + g.ppm + "_"));
      assert.ok(p.path.includes("_" + g.aerosol + "_"));
      assert.ok(!p.path.includes(".."));
      assert.ok(!path.isAbsolute(p.path));
    }
  }
  assert.equal(new Set(groups.flatMap((g) => g.panes.map((p) => p.path))).size, 160);
});

test("round and SIF routing stay isolated", () => {
  assert.match(plotRoot(3, true), /correlation_sif\/plots/);
  assert.match(plotRoot(3, false), /correlation_nosif\/plots/);
  for (const round of [4, 5, 6]) {
    assert.ok(plotRoot(round, true).startsWith("round" + round + "_"));
    assert.match(plotRoot(round, true), /plots_sif\//);
    assert.match(plotRoot(round, false), /plots_nosif\//);
  }
  assert.throws(() => plotRoot(7, true), RangeError);
});

test("navigation wraps without skipping groups", () => {
  assert.equal(wrapIndex(-1, 20), 19);
  assert.equal(wrapIndex(20, 20), 0);
  assert.equal(wrapIndex(43, 20), 3);
  assert.equal(wrapIndex(NaN, 20), 0);
});

test("all 160 current campaign PNGs exist and have PNG signatures", (t) => {
  const root = path.resolve(__dirname, "../../../bottom_layer_XCO2_retrievals");
  const generatedRoot = path.join(root, plotRoot(3, false));
  if (!fs.existsSync(generatedRoot)) {
    return t.skip("Generated campaign PNG workspace is not available");
  }
  for (const group of buildGroups()) {
    for (const pane of group.panes) {
      const file = path.join(root, pane.path);
      assert.ok(fs.statSync(file).isFile(), file);
      const fd = fs.openSync(file, "r");
      const header = Buffer.alloc(8);
      try { fs.readSync(fd, header, 0, 8, 0); } finally { fs.closeSync(fd); }
      assert.deepEqual(header, Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]), file);
    }
  }
});
