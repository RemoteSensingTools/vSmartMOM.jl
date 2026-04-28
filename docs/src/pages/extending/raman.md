# Add a Raman Mode

**For:** method developers extending inelastic-scattering support.

**Next:** [Core RT Theory](../vSmartMOM/CoreRTTheory.md), [API Reference](../api_reference.md).

This page will become the Raman-extension guide. It will cover:

- the `AbstractRamanType` hierarchy;
- how `noRS`, `RRS`, `VS_0to1`, `VS_1to0`, and `_plus` variants dispatch;
- where Raman optical properties are computed;
- which CoreRT elemental, doubling, and interaction methods are involved;
- current limitations, including linearized inelastic support.

Do not add a SIF-specific page here until the SIF fixture/data policy is resolved.
