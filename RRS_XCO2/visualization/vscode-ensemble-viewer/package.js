"use strict";

// Build a dependency-free VSIX using the standard VSIX container format.
// Nothing is published, and the generated archive contains code only, no plots.
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const { execFileSync } = require("node:child_process");
const pkg = require("./package.json");
const stage = fs.mkdtempSync(path.join(os.tmpdir(), "rrs-viewer-vsix-"));
const outputDir = path.join(__dirname, "dist");
const output = path.join(outputDir, pkg.name + "-" + pkg.version + ".vsix");
const xml = (text) => String(text).replace(/[&<>"']/g, (c) => ({
  "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&apos;",
})[c]);
const files = [
  "package.json", "extension.js", "model.js", "README.md",
  "media/viewer.js", "media/viewer.css",
];
for (const file of files) {
  const target = path.join(stage, "extension", file);
  fs.mkdirSync(path.dirname(target), { recursive: true });
  fs.copyFileSync(path.join(__dirname, file), target);
}
fs.writeFileSync(path.join(stage, "extension.vsixmanifest"), [
  '<?xml version="1.0" encoding="utf-8"?>',
  '<PackageManifest Version="2.0.0" xmlns="http://schemas.microsoft.com/developer/vsx-schema/2011" xmlns:d="http://schemas.microsoft.com/developer/vsx-schema-design/2011">',
  '<Metadata>',
  '<Identity Language="en-US" Id="' + xml(pkg.name) + '" Version="' +
    xml(pkg.version) + '" Publisher="' + xml(pkg.publisher) + '" />',
  '<DisplayName>' + xml(pkg.displayName) + '</DisplayName>',
  '<Description xml:space="preserve">' + xml(pkg.description) + '</Description>',
  '<Categories>Visualization</Categories><Tags>rrs,ensemble,plots</Tags>',
  '<Properties>',
  '<Property Id="Microsoft.VisualStudio.Code.Engine" Value="' + xml(pkg.engines.vscode) + '" />',
  '<Property Id="Microsoft.VisualStudio.Code.ExtensionKind" Value="workspace" />',
  '<Property Id="Microsoft.VisualStudio.Code.ExecutesCode" Value="true" />',
  '<Property Id="Microsoft.VisualStudio.Code.EnabledApiProposals" Value="" />',
  '</Properties></Metadata>',
  '<Installation><InstallationTarget Id="Microsoft.VisualStudio.Code" /></Installation>',
  '<Dependencies /><Assets>',
  '<Asset Type="Microsoft.VisualStudio.Code.Manifest" Path="extension/package.json" Addressable="true" />',
  '<Asset Type="Microsoft.VisualStudio.Services.Content.Details" Path="extension/README.md" Addressable="true" />',
  '</Assets></PackageManifest>',
].join("\n"));
fs.writeFileSync(path.join(stage, "[Content_Types].xml"), [
  '<?xml version="1.0" encoding="utf-8"?>',
  '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
  '<Default Extension="json" ContentType="application/json" />',
  '<Default Extension="js" ContentType="application/javascript" />',
  '<Default Extension="css" ContentType="text/css" />',
  '<Default Extension="md" ContentType="text/markdown" />',
  '<Default Extension="vsixmanifest" ContentType="text/xml" />',
  '</Types>',
].join("\n"));
fs.mkdirSync(outputDir, { recursive: true });
// zip builds a fresh archive inside the unique staging directory.
const stagedZip = path.join(stage, "viewer.vsix");
execFileSync("zip", ["-q", "-r", stagedZip, "extension",
  "extension.vsixmanifest", "[Content_Types].xml"], { cwd: stage });
execFileSync("unzip", ["-t", stagedZip], { stdio: "inherit" });
fs.copyFileSync(stagedZip, output);
// Remove only this build's uniquely created staging directory.
fs.rmSync(stage, { recursive: true });
console.log(output);
