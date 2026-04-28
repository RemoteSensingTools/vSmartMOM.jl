# vSmartMOM Module Overview

The vSmartMOM module wires parameters, model construction, and the solver into
end-to-end forward and linearized radiative-transfer runs.

For getting-started workflows see [Quick Start](../quickstart.md).
For Jacobians see [Compute Jacobians](../jacobians.md).
For equation-level solver mapping see [Core RT Theory](CoreRTTheory.md).

The core module is the reference-layer home for the `RTModel` hierarchy,
forward and linearized solver entry points, quadrature configuration, layer
operators, and surface coupling.

## Architecture

![Architecture diagram](vSmartMOMDiagram-vSmartMOM.drawio.png)
