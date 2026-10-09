#!/usr/bin/env python3
"""Canonical plotting view of legacy and round-4/5/6 retrieval states.

Round 3 stores the two native SIF coordinates ``SIF760`` and ``mSIF`` in
every 30-element retrieval state.  Round 4 fixes the SIF radiance at 759 nm:
SIF-off states contain neither SIF coordinate, while SIF-on states contain
only ``mSIF``.  Plotting code still benefits from a common 30-coordinate view,
so this module reconstructs the fixed/derived coordinates with the exact
round-4 affine map recorded in each NetCDF file.
Round 5 fixes both the 759-nm radiance and slope; neither is an active
coordinate, and both contribute zero to tangent increments. Round 6 has the
same SIF boundary but restores the original UTLS aerosol prior.

Absolute states and tangent increments deliberately use separate methods.
The known 759-nm intercept belongs in an absolute state, but never in a
``G @ delta_y`` increment.
"""

from pathlib import Path

import numpy as np


ROUND4_STATE_MODEL = "round4_known_sif759"
ROUND5_STATE_MODEL = "round5_fixed_sif"
ROUND6_STATE_MODEL = "round6_fixed_sif"
SIF760_NAME = "SIF760"
MSIF_NAME = "mSIF"
ROUND4_BASE_STATE_COUNT = 28
ROUND4_CORE_STATE_COUNT = 30
WAVENUMBER_CONVERSION_NM_CM1 = 1.0e7
ROUND4_KNOWN_WAVELENGTH_NM = 759.0
CORE_SIF_REFERENCE_WAVELENGTH_NM = 760.0
ROUND4_SIF_ON_CASE = "angular_integral760_0p5"
ROUND4_KNOWN_LNU759 = 0.004818031987713776
ROUND4_TRUTH_SIF760 = 0.004596394756493938
ROUND4_TRUTH_MSIF = 1.2291230681458325e-5
ROUND4_TRUTH_ANGULAR_INTEGRAL760 = 0.5
# ``true_states_corrected_sif_v2.dat`` writes mSIF with 13 significant
# decimal digits (1.229123068146e-05).  This absolute tolerance admits that
# ASCII round trip while remaining orders of magnitude below any physically
# meaningful slope perturbation.  NetCDF provenance checks remain exact to
# their independently specified tolerances.
ROUND4_TRUTH_MSIF_ASCII_ATOL = 5.0e-18
ROUND4_BASE_NAMES = (
    "psurf",
    *("co2_vmr_layer%02d" % layer for layer in range(5, 17)),
    "ln_sulfate_aod760",
    "ln_organic_carbon_aod760",
    "ln_utls_sulfate_aod760",
    "ln_sulfate_z0",
    "ln_organic_carbon_z0",
    "ln_utls_sulfate_z0",
    "o2a_surface_P0",
    "o2a_surface_P1",
    "o2a_surface_P2",
    "weak_co2_surface_P0",
    "weak_co2_surface_P1",
    "weak_co2_surface_P2",
    "strong_co2_surface_P0",
    "strong_co2_surface_P1",
    "strong_co2_surface_P2",
)


def _required_attribute(dataset, name, source):
    if name not in dataset.ncattrs():
        raise ValueError("%s is missing required attribute %s" % (source, name))
    return dataset.getncattr(name)


def _close(actual, expected, label, source, atol=2.0e-12, rtol=2.0e-12):
    actual = float(actual)
    if not np.isfinite(actual) or not np.isclose(
            actual, float(expected), rtol=rtol, atol=atol):
        raise ValueError(
            "%s has %s=%.17g; expected %.17g" %
            (source, label, actual, float(expected))
        )
    return actual


class RetrievalPlotSchema:
    """Validated mapping from a saved active state to plotting coordinates."""

    def __init__(self, active_names, state_model, sif_mode,
                 known_wavelength_nm=None, known_lnu=None,
                 core_reference_wavelength_nm=CORE_SIF_REFERENCE_WAVELENGTH_NM,
                 delta_nu=None, source="retrieval", fixed_msif=None):
        self.active_names = tuple(active_names)
        self.state_model = str(state_model)
        self.sif_mode = str(sif_mode)
        self.known_wavelength_nm = known_wavelength_nm
        self.known_lnu = known_lnu
        self.core_reference_wavelength_nm = core_reference_wavelength_nm
        self.delta_nu = delta_nu
        self.fixed_msif = fixed_msif
        self.source = str(source)

        if self.has_fixed_sif:
            base = tuple(self.active_names[:ROUND4_BASE_STATE_COUNT])
            if base != ROUND4_BASE_NAMES:
                raise ValueError(
                    "%s has a non-canonical %s base-state order" %
                    (self.source, self.state_model)
                )
            self.canonical_names = base + (SIF760_NAME, MSIF_NAME)
        else:
            self.canonical_names = self.active_names

    @property
    def is_round4(self):
        return self.state_model == ROUND4_STATE_MODEL

    @property
    def is_round5(self):
        return self.state_model == ROUND5_STATE_MODEL

    @property
    def is_round6(self):
        return self.state_model == ROUND6_STATE_MODEL

    @property
    def is_fixed_sif_model(self):
        """Round 5/6 models with neither SIF coefficient in the active state."""
        return self.is_round5 or self.is_round6

    @property
    def has_fixed_sif(self):
        return self.is_round4 or self.is_fixed_sif_model

    @property
    def sif760_status(self):
        if self.is_fixed_sif_model:
            return "derived from fixed SIF759 and fixed mSIF" \
                if self.sif_mode == "on" else "fixed to zero"
        if not self.is_round4:
            return "retrieved"
        return "derived from known SIF759 and retrieved mSIF" \
            if self.sif_mode == "on" else "fixed to zero"

    @property
    def msif_status(self):
        if self.is_fixed_sif_model:
            return "fixed" if self.sif_mode == "on" else "fixed to zero"
        if not self.is_round4:
            return "retrieved"
        return "retrieved" if self.sif_mode == "on" else "fixed to zero"

    @property
    def identity(self):
        """Tuple suitable for checking paired corrected/uncorrected files."""
        return (
            self.state_model, self.sif_mode, self.active_names,
            self.known_wavelength_nm, self.known_lnu,
            self.core_reference_wavelength_nm, self.delta_nu,
            self.fixed_msif,
        )

    @classmethod
    def from_dataset(cls, dataset, source=None):
        source = str(source or getattr(dataset, "filepath", lambda: "retrieval")())
        active_names = tuple(str(
            _required_attribute(dataset, "parameter_names", source)
        ).split())
        if not active_names:
            raise ValueError("%s has an empty parameter_names attribute" % source)
        advertised = int(dataset.getncattr("state_dimension")) \
            if "state_dimension" in dataset.ncattrs() else len(active_names)
        if advertised != len(active_names):
            raise ValueError(
                "%s advertises state_dimension=%d but names %d parameters" %
                (source, advertised, len(active_names))
            )

        state_model = str(dataset.getncattr("retrieval_state_model")) \
            if "retrieval_state_model" in dataset.ncattrs() else "legacy"
        if state_model in (ROUND5_STATE_MODEL, ROUND6_STATE_MODEL):
            return cls._from_fixed_sif_dataset(dataset, active_names, source,
                                               state_model)
        if state_model != ROUND4_STATE_MODEL:
            missing = [
                name for name in (SIF760_NAME, MSIF_NAME)
                if name not in active_names
            ]
            if missing:
                raise ValueError(
                    "%s is not marked as round 4 and is missing %s" %
                    (source, ", ".join(missing))
                )
            return cls(active_names, state_model, "legacy", source=source)

        enabled = int(_required_attribute(
            dataset, "round4_known_sif_enabled", source
        ))
        if enabled != 1:
            raise ValueError(
                "%s marks round4_known_sif_enabled=%d; expected 1" %
                (source, enabled)
            )
        sif_mode = str(_required_attribute(dataset, "round4_sif_case", source))
        if sif_mode not in ("off", "on"):
            raise ValueError(
                "%s has round4_sif_case=%r; expected off or on" %
                (source, sif_mode)
            )
        truth_case = str(_required_attribute(dataset, "sif_case", source))
        expected_truth_case = "off" if sif_mode == "off" else ROUND4_SIF_ON_CASE
        if truth_case != expected_truth_case:
            raise ValueError(
                "%s has sif_case=%r for round4_sif_case=%r; expected %r" %
                (source, truth_case, sif_mode, expected_truth_case)
            )

        known_wavelength = _close(
            _required_attribute(
                dataset, "round4_known_sif_wavelength_nm", source
            ),
            ROUND4_KNOWN_WAVELENGTH_NM,
            "round4_known_sif_wavelength_nm", source,
        )
        core_reference = _close(
            _required_attribute(
                dataset, "round4_core_sif_reference_wavelength_nm", source
            ),
            CORE_SIF_REFERENCE_WAVELENGTH_NM,
            "round4_core_sif_reference_wavelength_nm", source,
        )
        known_lnu = float(_required_attribute(
            dataset, "round4_known_sif_Lnu_mW_m-2_sr-1_per_cm-1", source
        ))
        if not np.isfinite(known_lnu) or known_lnu < 0.0:
            raise ValueError("%s has an invalid known SIF759 radiance" % source)
        expected_delta = (
            WAVENUMBER_CONVERSION_NM_CM1 / core_reference -
            WAVENUMBER_CONVERSION_NM_CM1 / known_wavelength
        )
        delta_nu = _close(
            _required_attribute(
                dataset, "round4_delta_nu_760_minus_759_cm-1", source
            ),
            expected_delta, "round4_delta_nu_760_minus_759_cm-1", source,
        )

        if "active_core_parameter_index" not in dataset.variables:
            raise ValueError(
                "%s is missing active_core_parameter_index" % source
            )
        active_core = tuple(int(value) for value in np.asarray(
            dataset["active_core_parameter_index"][:], dtype=int
        ).reshape(-1))
        expected_core = tuple(range(1, ROUND4_BASE_STATE_COUNT + 1))
        if sif_mode == "on":
            expected_core += (ROUND4_CORE_STATE_COUNT,)
            expected_names = (
                len(active_names) == 29 and
                active_names[-1] == MSIF_NAME and
                SIF760_NAME not in active_names
            )
        else:
            expected_names = (
                len(active_names) == 28 and
                SIF760_NAME not in active_names and
                MSIF_NAME not in active_names
            )
        if active_core != expected_core:
            raise ValueError(
                "%s has active_core_parameter_index=%s; expected %s" %
                (source, active_core, expected_core)
            )
        if not expected_names:
            raise ValueError(
                "%s has parameter names inconsistent with round-4 SIF-%s" %
                (source, sif_mode)
            )
        if sif_mode == "off" and not np.isclose(known_lnu, 0.0, atol=0.0):
            raise ValueError("%s is SIF-off but known SIF759 is nonzero" % source)
        if sif_mode == "on":
            _close(
                known_lnu, ROUND4_KNOWN_LNU759,
                "round4_known_sif_Lnu_mW_m-2_sr-1_per_cm-1", source,
                atol=3.0e-18, rtol=0.0,
            )

        return cls(
            active_names, state_model, sif_mode,
            known_wavelength_nm=known_wavelength,
            known_lnu=known_lnu,
            core_reference_wavelength_nm=core_reference,
            delta_nu=delta_nu,
            source=source,
        )

    @classmethod
    def _from_fixed_sif_dataset(cls, dataset, active_names, source, state_model):
        """Validate the saved 28-to-30 fixed-SIF boundary, without relabeling it."""
        prefix = "round6" if state_model == ROUND6_STATE_MODEL else "round5"
        if active_names != ROUND4_BASE_NAMES:
            raise ValueError("%s has non-canonical %s parameter names" % (source, prefix))
        mode = str(_required_attribute(dataset, prefix + "_sif_case", source))
        if mode not in ("off", "on"):
            raise ValueError("%s has invalid %s_sif_case=%r" % (source, prefix, mode))
        expected_truth = "off" if mode == "off" else ROUND4_SIF_ON_CASE
        if str(_required_attribute(dataset, "sif_case", source)) != expected_truth:
            raise ValueError("%s has sif_case inconsistent with %s_sif_case" % (source, prefix))
        if state_model == ROUND6_STATE_MODEL:
            name = "round6_stratospheric_aerosol_sigma_scale"
            _close(_required_attribute(dataset, name, source), 1.0,
                   name, source, atol=0.0, rtol=0.0)
        for name, expected in (
                (prefix + "_active_state_dimension", ROUND4_BASE_STATE_COUNT),
                (prefix + "_core_state_dimension", ROUND4_CORE_STATE_COUNT)):
            _close(_required_attribute(dataset, name, source), expected,
                   name, source, atol=0.0, rtol=0.0)
        expected_core = tuple(range(1, ROUND4_BASE_STATE_COUNT + 1))
        if "active_core_parameter_index" not in dataset.variables:
            raise ValueError("%s is missing active_core_parameter_index" % source)
        actual_core = tuple(np.asarray(
            dataset["active_core_parameter_index"][:]
        ).reshape(-1))
        if actual_core != expected_core:
            raise ValueError("%s has invalid active_core_parameter_index" % source)
        for name, expected in (
                (prefix + "_active_to_core_parameter_index", expected_core),
                (prefix + "_active_to_full_parameter_index", (1,) + tuple(range(6, 33)))):
            actual = tuple(int(value) for value in str(
                _required_attribute(dataset, name, source)
            ).split())
            if actual != expected:
                raise ValueError("%s has invalid %s" % (source, name))
        wavelength = _close(
            _required_attribute(dataset, prefix + "_known_sif_wavelength_nm", source),
            ROUND4_KNOWN_WAVELENGTH_NM, prefix + "_known_sif_wavelength_nm", source,
        )
        anchor_name = prefix + "_known_sif_Lnu_mW_m-2_sr-1_per_cm-1"
        anchor = _close(
            _required_attribute(dataset, anchor_name, source),
            ROUND4_KNOWN_LNU759 if mode == "on" else 0.0,
            anchor_name, source, atol=3.0e-18 if mode == "on" else 0.0, rtol=0.0,
        )
        slope_name = prefix + "_fixed_mSIF_mW_m-2_sr-1_per_cm-2"
        slope = _close(
            _required_attribute(dataset, slope_name, source),
            ROUND4_TRUTH_MSIF if mode == "on" else 0.0,
            slope_name, source, atol=3.0e-20 if mode == "on" else 0.0, rtol=0.0,
        )
        delta_nu = (WAVENUMBER_CONVERSION_NM_CM1 / CORE_SIF_REFERENCE_WAVELENGTH_NM -
                    WAVENUMBER_CONVERSION_NM_CM1 / wavelength)
        core_name = prefix + "_core_SIF760_mW_m-2_sr-1_per_cm-1"
        _close(_required_attribute(dataset, core_name, source),
               anchor + delta_nu * slope, core_name, source,
               atol=3.0e-18 if mode == "on" else 0.0, rtol=0.0)
        return cls(active_names, state_model, mode,
                   known_wavelength_nm=wavelength, known_lnu=anchor,
                   delta_nu=delta_nu, fixed_msif=slope, source=source)

    def _validated_array(self, values, label):
        array = np.asarray(values, dtype=float)
        if array.ndim < 1 or array.shape[-1] != len(self.active_names):
            raise ValueError(
                "%s %s has trailing dimension %s; expected %d" %
                (self.source, label, array.shape, len(self.active_names))
            )
        if not np.all(np.isfinite(array)):
            raise ValueError("%s %s contains non-finite values" % (
                self.source, label,
            ))
        return array

    def expand_absolute_state(self, values):
        """Return absolute state(s) in the canonical 30-coordinate view."""
        active = self._validated_array(values, "state")
        if not self.has_fixed_sif:
            return np.array(active, copy=True)
        output = np.zeros(
            active.shape[:-1] + (ROUND4_CORE_STATE_COUNT,), dtype=float
        )
        output[..., :ROUND4_BASE_STATE_COUNT] = \
            active[..., :ROUND4_BASE_STATE_COUNT]
        if self.sif_mode == "on":
            msif = self.fixed_msif if self.is_fixed_sif_model else active[..., -1]
            output[..., -2] = self.known_lnu + self.delta_nu * msif
            output[..., -1] = msif
        return output

    def expand_tangent(self, values):
        """Return tangent increment(s), excluding the fixed 759-nm intercept."""
        active = self._validated_array(values, "tangent increment")
        if not self.has_fixed_sif:
            return np.array(active, copy=True)
        output = np.zeros(
            active.shape[:-1] + (ROUND4_CORE_STATE_COUNT,), dtype=float
        )
        output[..., :ROUND4_BASE_STATE_COUNT] = \
            active[..., :ROUND4_BASE_STATE_COUNT]
        if self.is_round4 and self.sif_mode == "on":
            delta_msif = active[..., -1]
            output[..., -2] = self.delta_nu * delta_msif
            output[..., -1] = delta_msif
        return output

    def native_dict(self, values, tangent=False):
        expanded = self.expand_tangent(values) if tangent \
            else self.expand_absolute_state(values)
        if expanded.ndim != 1:
            raise ValueError("native_dict requires one state vector")
        result = dict(zip(self.canonical_names, expanded))
        if self.has_fixed_sif:
            result["SIF759"] = 0.0 if tangent else self.known_lnu
        else:
            # Lnu759 = Lnu760 + mSIF*(nu759-nu760).
            result["SIF759"] = (
                result[SIF760_NAME] - self._legacy_delta_nu() *
                result[MSIF_NAME]
            )
        return result

    def truth_sif_coordinates(self, row, source="truth-table row"):
        """Return truth SIF coordinates in this retrieval's state space.

        The nonlinear truth template independently tabulates its value at
        760 nm.  Round 4, however, fixes the exact template value at 759 nm
        and retains only a linear spectral slope in the retrieval state.  Its
        comparison truth must therefore use that same affine state-space map.
        """
        self.validate_truth_row(row, source)
        try:
            native_sif760 = float(row[SIF760_NAME])
            native_msif = float(row[MSIF_NAME])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                "%s is missing finite SIF760/mSIF truth coordinates" % source
            ) from error
        if not np.isfinite(native_sif760) or not np.isfinite(native_msif):
            raise ValueError(
                "%s has non-finite SIF760/mSIF truth coordinates" % source
            )

        if not self.has_fixed_sif:
            return {
                "SIF759": native_sif760 - self._legacy_delta_nu() * native_msif,
                SIF760_NAME: native_sif760,
                MSIF_NAME: native_msif,
            }
        if self.sif_mode == "off":
            return {"SIF759": 0.0, SIF760_NAME: 0.0, MSIF_NAME: 0.0}

        _close(
            native_msif, ROUND4_TRUTH_MSIF, "mSIF", source,
            atol=ROUND4_TRUTH_MSIF_ASCII_ATOL, rtol=0.0,
        )
        if self.is_fixed_sif_model:
            native_msif = self.fixed_msif
        return {
            "SIF759": self.known_lnu,
            SIF760_NAME: self.known_lnu + self.delta_nu * native_msif,
            MSIF_NAME: native_msif,
        }

    def _legacy_delta_nu(self):
        return (
            WAVENUMBER_CONVERSION_NM_CM1 /
            CORE_SIF_REFERENCE_WAVELENGTH_NM -
            WAVENUMBER_CONVERSION_NM_CM1 /
            ROUND4_KNOWN_WAVELENGTH_NM
        )

    def description(self):
        if self.is_fixed_sif_model:
            return "round %d: SIF759 fixed at %.9g; mSIF fixed at %.9g; SIF760 derived" % (
                6 if self.is_round6 else 5, self.known_lnu, self.fixed_msif,
            )
        if not self.is_round4:
            return "legacy 30-coordinate state: SIF760 and mSIF retrieved"
        if self.sif_mode == "off":
            return "round 4: SIF759 and mSIF fixed to zero"
        return (
            "round 4: SIF759 fixed at %.9g; mSIF retrieved; SIF760 derived" %
            self.known_lnu
        )

    def validate_truth_row(self, row, source="truth-table row"):
        """Reject a truth row from a different SIF campaign convention.

        The round-4 SIF-on forward model is anchored to the exact nonlinear
        truth-template value at 759 nm.  Consequently we validate the campaign
        label here, but do not demand that a tangent inferred from the row's
        760-nm value reproduce that exact anchor.
        """
        if not self.has_fixed_sif:
            return
        truth_case = str(row.get("sif_case", ""))
        expected = "off" if self.sif_mode == "off" else ROUND4_SIF_ON_CASE
        if truth_case != expected:
            raise ValueError(
                "%s has sif_case=%r, but %s requires %r; select the exact "
                "fixed-SIF truth table instead of a legacy SIF table" %
                (source, truth_case, self.source, expected)
            )
        if self.sif_mode == "off":
            for name in (SIF760_NAME, MSIF_NAME):
                if name in row and not np.isclose(float(row[name]), 0.0,
                                                  rtol=0.0, atol=0.0):
                    raise ValueError(
                        "%s is SIF-off but %s is nonzero" % (source, name)
                    )
        else:
            required = (
                SIF760_NAME, MSIF_NAME, "sif_angular_integral760",
            )
            missing = [name for name in required if name not in row]
            if missing:
                raise ValueError(
                    "%s is missing corrected-v2 SIF field(s): %s" %
                    (source, ", ".join(missing))
                )
            _close(
                row[SIF760_NAME], ROUND4_TRUTH_SIF760, SIF760_NAME, source,
                atol=5.0e-15, rtol=0.0,
            )
            _close(
                row[MSIF_NAME], ROUND4_TRUTH_MSIF, MSIF_NAME, source,
                atol=ROUND4_TRUTH_MSIF_ASCII_ATOL, rtol=0.0,
            )
            _close(
                row["sif_angular_integral760"],
                ROUND4_TRUTH_ANGULAR_INTEGRAL760,
                "sif_angular_integral760", source,
                atol=0.0, rtol=0.0,
            )


def default_round4_local_root(campaign_root):
    """Return the established local round-4 no-SIF retrieval directory."""
    return Path(campaign_root) / "round4_known_sif759" / "retrievals_nosif"
