# RRS Ensemble Viewer

A private, dependency-free VS Code extension for the existing ensemble PNGs.
It runs in the remote workspace extension host on Curry. No plots or NetCDF
data are bundled, uploaded, or modified. It does not run retrievals or shell
commands, fetch internet content, or collect telemetry.

## Open

In the Curry-connected VS Code window, open the Command Palette and run
**RRS: Open Ensemble Comparison**.

The viewer finds the repository in the open workspace. If it cannot find it,
select the **bottom_layer_XCO2_retrievals** directory when prompted.

## Layout and controls

| | Round 3 | Round 4 | Round 5 | Round 6 |
|---|---|---|---|---|
| Top | Aerosol + surface | Aerosol + surface | Aerosol + surface | Aerosol + surface |
| Bottom | CO2 + SIF | CO2 + SIF | CO2 + SIF | CO2 + SIF |

Each image retains all four surfaces: urban, rural, desert, and forest.

- Click inside the viewer, then **Left/Right** to cycle all eight panes together.
- There are 20 groups: five CO2 levels, two SIF cases, and two aerosol cases.
- Order follows truth state IDs: 001/021/041/061, 002/022/042/062, through
  020/040/060/080. Navigation wraps at either end.
- Use the dropdown to jump directly to a group.
- Click an image to enlarge it; click again or press **Escape** to restore
  all eight panes. Arrow navigation also works while enlarged.
- **Refresh** rechecks files after plots are regenerated.
- Unavailable images have placeholders; rounds are never shifted or substituted.

All eight images are committed together after loading. Rapid arrow presses
cannot display images from an older group after a newer group has loaded.
The selected group is remembered for this workspace.

For legibility, maximize the editor area or hide the sidebars. Use image
enlargement when eight figures are too small to inspect at once.

## Build and install

From this directory, run **npm test**, then **npm run package**. There are no
npm dependencies to install. Packaging requires the standard zip/unzip tools.
The archive is written to **dist/rrs-ensemble-viewer-0.1.0.vsix**.

Install from a VS Code terminal connected to Curry with:

    code --install-extension /absolute/path/to/rrs-ensemble-viewer-0.1.0.vsix

Alternatively use **Extensions: Install from VSIX...** in the Curry-connected
window. If the command is not listed immediately after installation, run
**Developer: Reload Window** and try again.

Extension ID: **sanghavi-local.rrs-ensemble-viewer**.
To remove it, uninstall **RRS Ensemble Viewer** from the Extensions view.

## Implementation references

Uses the official [VS Code webview API](https://code.visualstudio.com/api/extension-guides/webview)
with a restricted content-security policy and local resource roots, and
[workspace extension hosting](https://code.visualstudio.com/api/advanced-topics/remote-extensions)
for Remote SSH.
