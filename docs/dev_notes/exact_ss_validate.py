"""
exact_ss_validate.py — Python translation of exact_ss_reference.jl for
validating the math in this environment (where Julia isn't available).

Mirrors the Julia file structure exactly so that bug-fixes here can be
ported back. Use as a cross-check, not as the primary reference.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss


# ============================================================================
# Contributor abstraction
# ============================================================================

class RayleighContributor:
    def __init__(self, tau):
        self.tau = np.asarray(tau, dtype=float)

    def phase(self, cos_theta):
        return 0.75 * (1.0 + cos_theta**2)

    def optical_depth(self, iz):
        return self.tau[iz]

    def ssa(self):
        return 1.0


class HGAerosolContributor:
    def __init__(self, g, omega, tau):
        self.g = g
        self.omega = omega
        self.tau = np.asarray(tau, dtype=float)

    def phase(self, cos_theta):
        g = self.g
        return (1.0 - g**2) / (1.0 + g**2 - 2*g * cos_theta)**1.5

    def optical_depth(self, iz):
        return self.tau[iz]

    def ssa(self):
        return self.omega


class AbsorptionContributor:
    def __init__(self, tau):
        self.tau = np.asarray(tau, dtype=float)

    def phase(self, cos_theta):
        return 0.0

    def optical_depth(self, iz):
        return self.tau[iz]

    def ssa(self):
        return 0.0


# ============================================================================
# Layer-effective optical properties
# ============================================================================

def layer_total_tau(contributors, iz):
    return sum(c.optical_depth(iz) for c in contributors)


def layer_effective_omega(contributors, iz):
    tau_total = layer_total_tau(contributors, iz)
    if tau_total == 0:
        return 0.0
    tau_scat = sum(c.optical_depth(iz) * c.ssa() for c in contributors)
    return tau_scat / tau_total


def layer_effective_phase(contributors, iz, cos_theta):
    weight_sum = 0.0
    weighted = 0.0
    for c in contributors:
        omega = c.ssa()
        if omega == 0:
            continue
        w = c.optical_depth(iz) * omega
        weight_sum += w
        weighted += w * c.phase(cos_theta)
    if weight_sum == 0:
        return 0.0
    return weighted / weight_sum


# ============================================================================
# Geometry
# ============================================================================

def scattering_angle_cosine(mu0, mu_v, dphi):
    return -mu0 * mu_v + np.sqrt(1 - mu0**2) * np.sqrt(1 - mu_v**2) * np.cos(dphi)


def cumulative_tau(contributors):
    n_layers = len(contributors[0].tau)
    tau_cum = np.zeros(n_layers + 1)
    for iz in range(n_layers):
        tau_cum[iz + 1] = tau_cum[iz] + layer_total_tau(contributors, iz)
    return tau_cum


# ============================================================================
# Quadrature on (0, 1]
# ============================================================================

def gauss_legendre_01(n):
    """Gauss-Legendre on (0,1]. Returns (nodes, weights)."""
    x_std, w_std = leggauss(n)  # on [-1, 1]
    x = 0.5 * (x_std + 1.0)
    w = 0.5 * w_std
    return x, w


# ============================================================================
# PATH 1
# ============================================================================

def path1_atmospheric_ss(contributors, mu0, mu_v, dphi, I0):
    assert mu0 > 0 and mu_v > 0
    cos_theta = scattering_angle_cosine(mu0, mu_v, dphi)
    tau_cum = cumulative_tau(contributors)
    n_layers = len(tau_cum) - 1
    a = 1.0/mu0 + 1.0/mu_v
    # Factor 1/mu_v from the RTE solution: integrating source along upward
    # slant path. dτ_slant = dτ_vertical / μ_v.
    prefactor = I0 / (4*np.pi * mu_v) / a

    L = 0.0
    for iz in range(n_layers):
        omega_eff = layer_effective_omega(contributors, iz)
        if omega_eff == 0:
            continue
        P_eff = layer_effective_phase(contributors, iz, cos_theta)
        layer_factor = np.exp(-tau_cum[iz] * a) - np.exp(-tau_cum[iz+1] * a)
        L += prefactor * omega_eff * P_eff * layer_factor
    return L


# ============================================================================
# PATH 2
# ============================================================================

def path2_surface_direct(tau_total, mu0, mu_v, albedo, I0):
    return (mu0 * I0 * albedo / np.pi) * np.exp(-tau_total/mu0) * np.exp(-tau_total/mu_v)


# ============================================================================
# PATHS 3 and 4: atm-surface coupled
# ============================================================================

def azimuthal_average_phase(contributor, mu_d, mu0, N_phi=64):
    """(1/2π) ∫₀^{2π} P(cos Θ(μ_d, μ₀, ϕ')) dϕ' via trapezoidal rule on periodic interval."""
    # cos Θ = μ₀ μ_d + sqrt(1-μ₀²) sqrt(1-μ_d²) cos(ϕ')
    a = mu0 * mu_d
    b = np.sqrt(1 - mu0**2) * np.sqrt(1 - mu_d**2)
    s = 0.0
    for k in range(N_phi):
        phi = 2*np.pi * k / N_phi
        cos_t = np.clip(a + b * np.cos(phi), -1.0, 1.0)
        s += contributor.phase(cos_t)
    return s / N_phi


def azimuthal_average_phase_rayleigh_analytic(mu_d, mu0):
    return 0.75 * (1.0 + (mu0 * mu_d)**2 + 0.5 * (1 - mu0**2) * (1 - mu_d**2))


def layer_effective_azimuthal_phase(contributors, iz, mu_d, mu0, N_phi=64):
    weight_sum = 0.0
    weighted = 0.0
    for c in contributors:
        omega = c.ssa()
        if omega == 0:
            continue
        w = c.optical_depth(iz) * omega
        weight_sum += w
        P_avg = azimuthal_average_phase(c, mu_d, mu0, N_phi=N_phi)
        weighted += w * P_avg
    if weight_sum == 0:
        return 0.0
    return weighted / weight_sum


def path3_layer_tau_integral(tau_top, tau_bot, tau_total, mu0, mu_d):
    """∫_{tau_top}^{tau_bot} exp(-tau/mu0) exp(-(tau_total - tau)/mu_d) dtau

    Numerically stable form: combine exponentials at each endpoint before
    subtracting, to avoid overflow when mu_d is small.
    """
    b = 1.0/mu0 - 1.0/mu_d
    # Endpoint values of the integrand (without the 1/b factor):
    f_top = np.exp(-tau_top/mu0 - (tau_total - tau_top)/mu_d)
    f_bot = np.exp(-tau_bot/mu0 - (tau_total - tau_bot)/mu_d)
    if abs(b) < 1e-10:
        # b ≈ 0: integrand is approximately constant, use trapezoidal
        return 0.5 * (f_top + f_bot) * (tau_bot - tau_top)
    return (f_top - f_bot) / b


def path3_atm_to_surface(contributors, mu0, mu_v, albedo, I0, N_quad=16, N_phi=64):
    tau_cum = cumulative_tau(contributors)
    tau_total = tau_cum[-1]
    n_layers = len(tau_cum) - 1
    mu_nodes, mu_weights = gauss_legendre_01(N_quad)

    F_surface = 0.0
    for iz in range(n_layers):
        omega_eff = layer_effective_omega(contributors, iz)
        if omega_eff == 0:
            continue
        tau_top = tau_cum[iz]
        tau_bot = tau_cum[iz+1]

        layer_sum = 0.0
        for k in range(N_quad):
            mu_d = mu_nodes[k]
            w_d = mu_weights[k]
            P_bar = layer_effective_azimuthal_phase(contributors, iz, mu_d, mu0, N_phi=N_phi)
            tau_int = path3_layer_tau_integral(tau_top, tau_bot, tau_total, mu0, mu_d)
            # Factor 0.5 = 2π/4π from azimuthal normalization.
            # No μ_d: the 1/μ_d (RTE attenuation along slant) and μ_d (radiance →
            # irradiance) cancel, leaving the integrand bare.
            layer_sum += w_d * P_bar * tau_int * 0.5
        F_surface += omega_eff * I0 * layer_sum

    L = (albedo / np.pi) * F_surface * np.exp(-tau_total/mu_v)
    return L


def path4_surface_to_atm(contributors, mu0, mu_v, albedo, I0, N_quad=16, N_phi=64):
    tau_cum = cumulative_tau(contributors)
    tau_total = tau_cum[-1]
    n_layers = len(tau_cum) - 1
    mu_nodes, mu_weights = gauss_legendre_01(N_quad)

    L = 0.0
    L_surface = (albedo / np.pi) * mu0 * I0 * np.exp(-tau_total/mu0)

    for iz in range(n_layers):
        omega_eff = layer_effective_omega(contributors, iz)
        if omega_eff == 0:
            continue
        tau_top = tau_cum[iz]
        tau_bot = tau_cum[iz+1]

        layer_sum = 0.0
        for k in range(N_quad):
            mu_u = mu_nodes[k]
            w_u = mu_weights[k]
            P_bar = layer_effective_azimuthal_phase(contributors, iz, mu_u, mu_v, N_phi=N_phi)
            b = 1.0/mu_v - 1.0/mu_u
            f_top = np.exp(-tau_top/mu_v - (tau_total - tau_top)/mu_u)
            f_bot = np.exp(-tau_bot/mu_v - (tau_total - tau_bot)/mu_u)
            if abs(b) < 1e-10:
                tau_int = 0.5 * (f_top + f_bot) * (tau_bot - tau_top)
            else:
                tau_int = (f_top - f_bot) / b
            # No μ_u: the source-into-view integral over incoming directions
            # carries no μ_u factor (radiance source from radiance input).
            layer_sum += w_u * P_bar * tau_int * 0.5
        # 1/μ_v from RTE solution along upward view path.
        L += L_surface * omega_eff * layer_sum / mu_v

    return L


def exact_ss_toa(contributors, mu0, mu_v, dphi, albedo, I0, N_quad=16, N_phi=64):
    tau_cum = cumulative_tau(contributors)
    tau_total = tau_cum[-1]
    p1 = path1_atmospheric_ss(contributors, mu0, mu_v, dphi, I0)
    p2 = path2_surface_direct(tau_total, mu0, mu_v, albedo, I0)
    p3 = path3_atm_to_surface(contributors, mu0, mu_v, albedo, I0, N_quad, N_phi)
    p4 = path4_surface_to_atm(contributors, mu0, mu_v, albedo, I0, N_quad, N_phi)
    return {'path1': p1, 'path2': p2, 'path3': p3, 'path4': p4,
            'total': p1 + p2 + p3 + p4}


# ============================================================================
# Tests
# ============================================================================

def assert_close(actual, expected, rtol=1e-12, atol=0.0, msg=""):
    if abs(expected) < atol + 1e-300:
        ok = abs(actual) <= atol
    else:
        ok = abs(actual - expected) <= rtol * abs(expected) + atol
    if not ok:
        raise AssertionError(
            f"{msg}: actual={actual}, expected={expected}, "
            f"diff={actual - expected}, rtol={rtol}, atol={atol}"
        )


def test_geometry():
    assert_close(scattering_angle_cosine(0.5, 0.5, 0.0), 0.5,       msg="μ₀=μ_v=0.5, Δϕ=0")
    assert_close(scattering_angle_cosine(0.5, 0.5, np.pi), -1.0,    msg="μ₀=μ_v=0.5, Δϕ=π")
    assert_close(scattering_angle_cosine(1.0, 1.0, 0.0), -1.0,      msg="zenith")
    assert_close(scattering_angle_cosine(1.0, 1.0, np.pi), -1.0,    msg="zenith Δϕ=π")
    print("  geometry: OK")


def test_phase_normalization():
    # Rayleigh
    n_q = 100000
    mus = np.linspace(-1, 1, n_q, endpoint=False) + 1.0/n_q
    ray = RayleighContributor([0.5])
    norm = np.mean([ray.phase(m) for m in mus])
    assert_close(norm, 1.0, rtol=1e-3, msg="Rayleigh normalization")

    for g in [0.0, 0.3, 0.6, -0.4]:
        hg = HGAerosolContributor(g, 1.0, [0.1])
        norm = np.mean([hg.phase(m) for m in mus])
        assert_close(norm, 1.0, rtol=1e-3, msg=f"HG g={g} normalization")

    hg_iso = HGAerosolContributor(0.0, 1.0, [0.1])
    for c in [-1.0, -0.3, 0.0, 0.5, 1.0]:
        assert_close(hg_iso.phase(c), 1.0, msg=f"HG g=0 isotropic at cosΘ={c}")
    print("  phase normalization: OK")


def test_rayleigh_azimuthal_average():
    ray = RayleighContributor([0.5])
    for mu_d in [0.1, 0.4, 0.8]:
        for mu0 in [0.2, 0.5, 0.9]:
            num = azimuthal_average_phase(ray, mu_d, mu0, N_phi=128)
            ana = azimuthal_average_phase_rayleigh_analytic(mu_d, mu0)
            assert_close(num, ana, rtol=1e-10,
                         msg=f"Rayleigh azimuthal avg μ_d={mu_d}, μ₀={mu0}")
    print("  Rayleigh azimuthal avg analytic vs numerical: OK")


def test_gauss_legendre():
    for n in [4, 8, 16]:
        x, w = gauss_legendre_01(n)
        assert_close(np.sum(w), 1.0, rtol=1e-12, msg=f"GL n={n} ∫1=1")
        assert_close(np.sum(w*x), 0.5, rtol=1e-12, msg=f"GL n={n} ∫x=0.5")
        for k in range(1, 2*n):
            assert_close(np.sum(w * x**k), 1.0/(k+1), rtol=1e-10,
                         msg=f"GL n={n} ∫x^{k}")
    print("  Gauss-Legendre: OK")


def test_path1_zenith_rayleigh():
    # μ₀=μ_v=1, Δϕ=0, P_Ray(-1)=1.5, single layer τ=0.1.
    # Expected: (1/(4π·μ_v)) · ϖ · P · [1 - exp(-τ·a)]/a
    #         = (1/4π · 1) · 1 · 1.5 · (1 - exp(-0.2)) / 2
    # μ_v=1 so this passes either way; not a sufficient test for the 1/μ_v factor.
    ray = RayleighContributor([0.1])
    r = path1_atmospheric_ss([ray], 1.0, 1.0, 0.0, 1.0)
    expected = (1.5 / (4*np.pi)) * (1.0 - np.exp(-0.2)) / 2.0
    assert_close(r, expected, rtol=1e-12, msg="path1 zenith Rayleigh")
    print(f"  path1 zenith Rayleigh: {r:.6e} (expected {expected:.6e})  OK")


def test_path1_offzenith_isotropic():
    """Off-zenith isotropic scattering, the test GPT recommended that does
    catch the 1/μ_v factor since μ_v ≠ 1."""
    # Isotropic: HG with g=0, ϖ=1, P=1.
    # μ₀=0.5, μ_v=0.25, τ=0.1.
    # cos Θ = -0.5·0.25 + sqrt(0.75)·sqrt(0.9375)·1 = -0.125 + 0.838... = 0.713...
    # P=1 (isotropic), so the cos Θ doesn't matter.
    # Expected: (1/(4π·μ_v)) · ϖ · P · [1 - exp(-τ·a)] / a
    iso = HGAerosolContributor(0.0, 1.0, [0.1])
    mu0, mu_v = 0.5, 0.25
    tau = 0.1
    a = 1.0/mu0 + 1.0/mu_v
    expected = (1.0 / (4*np.pi * mu_v)) * (1.0 - np.exp(-tau * a)) / a
    r = path1_atmospheric_ss([iso], mu0, mu_v, 0.0, 1.0)
    assert_close(r, expected, rtol=1e-12,
                 msg="path1 off-zenith isotropic (1/μ_v factor test)")
    print(f"  path1 off-zenith isotropic: {r:.6e} (expected {expected:.6e})  OK")


def test_path1_layer_subdivision():
    r1 = RayleighContributor([0.10])
    r2 = RayleighContributor([0.04, 0.06])
    for mu0 in [0.3, 0.7]:
        for mu_v in [0.4, 0.9]:
            for dphi in [0.0, np.pi/2, np.pi]:
                v1 = path1_atmospheric_ss([r1], mu0, mu_v, dphi, 1.0)
                v2 = path1_atmospheric_ss([r2], mu0, mu_v, dphi, 1.0)
                assert_close(v1, v2, rtol=1e-12,
                             msg=f"layer subdiv μ₀={mu0}, μ_v={mu_v}, Δϕ={dphi}")
    print("  path1 layer subdivision invariance: OK")


def test_mixed_contributors():
    ray = RayleighContributor([0.05])
    hg  = HGAerosolContributor(0.6, 0.9, [0.1])
    cos_t = -0.7
    omega_eff = layer_effective_omega([ray, hg], 0)
    P_eff = layer_effective_phase([ray, hg], 0, cos_t)

    omega_exp = (0.05*1.0 + 0.1*0.9) / 0.15
    P_ray_at = 0.75 * (1 + cos_t**2)
    g = 0.6
    P_hg_at = (1 - g**2) / (1 + g**2 - 2*g*cos_t)**1.5
    P_exp = (0.05*1.0*P_ray_at + 0.1*0.9*P_hg_at) / (0.05*1.0 + 0.1*0.9)

    assert_close(omega_eff, omega_exp, msg="ϖ_eff Rayleigh+HG")
    assert_close(P_eff, P_exp, msg="P_eff Rayleigh+HG")
    print("  mixed Rayleigh+HG contributor consistency: OK")


def test_path2_no_atmosphere():
    empty = RayleighContributor([0.0])
    r = exact_ss_toa([empty], 0.5, 0.7, 0.3, 0.4, 1.0)
    expected = 0.5 * 1.0 * 0.4 / np.pi
    assert_close(r['path2'], expected, msg="path2 no atm")
    assert_close(r['path1'], 0.0, atol=1e-15, msg="path1 should be 0")
    assert_close(r['path3'], 0.0, atol=1e-15, msg="path3 should be 0")
    assert_close(r['path4'], 0.0, atol=1e-15, msg="path4 should be 0")
    print(f"  path2 no atmosphere: {r['path2']:.6e}  OK")


def test_black_surface():
    ray = RayleighContributor([0.2])
    r = exact_ss_toa([ray], 0.6, 0.4, np.pi/3, 0.0, 1.0)
    assert r['path2'] == 0.0, f"path2 should be 0, got {r['path2']}"
    assert r['path3'] == 0.0, f"path3 should be 0, got {r['path3']}"
    assert r['path4'] == 0.0, f"path4 should be 0, got {r['path4']}"
    assert_close(r['total'], r['path1'], msg="total = path1 black surf")
    print("  black surface: OK")


def test_absorption_only_atm():
    abs_c = AbsorptionContributor([0.3])
    r = exact_ss_toa([abs_c], 0.5, 0.7, 0.0, 0.5, 1.0)
    assert_close(r['path1'], 0.0, atol=1e-15, msg="path1 abs only")
    assert_close(r['path3'], 0.0, atol=1e-15, msg="path3 abs only")
    assert_close(r['path4'], 0.0, atol=1e-15, msg="path4 abs only")
    expected_p2 = 0.5 * 1.0 * 0.5 / np.pi * np.exp(-0.3/0.5) * np.exp(-0.3/0.7)
    assert_close(r['path2'], expected_p2, msg="path2 abs only")
    print(f"  absorption-only atmosphere: path2={r['path2']:.6e}  OK")


def test_path3_path4_brdf_reciprocity():
    """With corrected formulas, the sum (path3 + path4) should satisfy BRDF
    reciprocity. By the structure of the corrections, path3(μ₀, μ_v)/μ₀ should
    equal path4(μ_v, μ₀)/μ_v — i.e., swapping sun and view directions and
    swapping the role of "atmosphere first" vs "surface first" maps the
    contributions onto each other.
    """
    ray = RayleighContributor([0.15])
    hg  = HGAerosolContributor(0.4, 0.95, [0.05])
    mu0, mu_v = 0.6, 0.4
    albedo = 0.3
    r_orig = exact_ss_toa([ray, hg], mu0, mu_v, 0.0, albedo, 1.0,
                          N_quad=64, N_phi=128)
    r_swap = exact_ss_toa([ray, hg], mu_v, mu0, 0.0, albedo, 1.0,
                          N_quad=64, N_phi=128)
    print(f"  path3+path4 reciprocity test:")
    print(f"    orig: path3={r_orig['path3']:.6e}, path4={r_orig['path4']:.6e}")
    print(f"    swap: path3={r_swap['path3']:.6e}, path4={r_swap['path4']:.6e}")
    print(f"    orig.path3/μ₀ = {r_orig['path3']/mu0:.6e}, swap.path4/μ_v = {r_swap['path4']/mu_v:.6e}")
    print(f"    orig.path4/μ₀ = {r_orig['path4']/mu0:.6e}, swap.path3/μ_v = {r_swap['path3']/mu_v:.6e}")
    print(f"    total/μ₀ = {r_orig['total']/mu0:.6e}, total_swap/μ_v = {r_swap['total']/mu_v:.6e}")
    # Full BRDF reciprocity: total/μ₀ symmetric under (μ₀, μ_v) swap
    assert_close(r_orig['total']/mu0, r_swap['total']/mu_v, rtol=1e-4,
                 msg="full BRDF reciprocity")
    # Path-by-path reciprocity (paths 3 and 4 swap roles under (μ₀, μ_v) swap)
    assert_close(r_orig['path3']/mu0, r_swap['path4']/mu_v, rtol=1e-4,
                 msg="path3(μ₀,μ_v)/μ₀ ≈ path4(μ_v,μ₀)/μ_v")
    assert_close(r_orig['path4']/mu0, r_swap['path3']/mu_v, rtol=1e-4,
                 msg="path4(μ₀,μ_v)/μ₀ ≈ path3(μ_v,μ₀)/μ_v")
    print("  paths 3+4 BRDF reciprocity: OK")


def test_path1_brdf_reciprocity():
    """With the corrected 1/μ_v prefactor, path 1 satisfies BRDF reciprocity:
        path1(μ₀, μ_v, Δϕ) / μ₀ = path1(μ_v, μ₀, Δϕ) / μ_v
    Equivalently: path1 · (1/μ₀ - 1/μ_v) is symmetric under swap
    multiplied by appropriate factors.
    """
    ray = RayleighContributor([0.15])
    hg  = HGAerosolContributor(0.4, 0.95, [0.05])
    for (mu0, mu_v) in [(0.6, 0.4), (0.3, 0.8), (0.1, 0.9)]:
        for dphi in [0.0, np.pi/3, np.pi]:
            p_orig = path1_atmospheric_ss([ray, hg], mu0, mu_v, dphi, 1.0)
            p_swap = path1_atmospheric_ss([ray, hg], mu_v, mu0, dphi, 1.0)
            # BRDF reciprocity: f_r(μ₀, μ_v) = f_r(μ_v, μ₀)
            # where f_r ∝ I/μ₀, so I(μ₀,μ_v)/μ₀ = I(μ_v,μ₀)/μ_v
            assert_close(p_orig/mu0, p_swap/mu_v, rtol=1e-12,
                         msg=f"path1 reciprocity (μ₀={mu0}, μ_v={mu_v}, Δϕ={dphi})")
    print("  path1 BRDF reciprocity: OK")


def test_path2_normalized_symmetry():
    """path2(μ₀, μ_v) / (μ₀ μ_v) is symmetric — the prefactor μ₀ in path2's
    definition is the only asymmetric piece. Verifies the closed-form structure."""
    tau_total = 0.2
    albedo = 0.3
    for (mu0, mu_v) in [(0.6, 0.4), (0.3, 0.8), (0.5, 0.5)]:
        p_orig = path2_surface_direct(tau_total, mu0, mu_v, albedo, 1.0)
        p_swap = path2_surface_direct(tau_total, mu_v, mu0, albedo, 1.0)
        # path2(μ₀, μ_v) = μ₀ · g(τ_total/μ₀, τ_total/μ_v) by the formula structure
        # Test: path2(μ₀,μ_v)/μ₀ = path2(μ_v,μ₀)/μ_v
        assert_close(p_orig/mu0, p_swap/mu_v, rtol=1e-12,
                     msg=f"path2 normalized (μ₀={mu0}, μ_v={mu_v})")
    print("  path2 normalized symmetry: OK")


def test_quadrature_convergence():
    ray = RayleighContributor([0.2])
    hg  = HGAerosolContributor(0.5, 0.9, [0.1])
    results = []
    for N_q in [8, 16, 32, 64]:
        r = exact_ss_toa([ray, hg], 0.5, 0.6, np.pi/4, 0.3, 1.0,
                         N_quad=N_q, N_phi=128)
        results.append(r['total'])
    print(f"  convergence (N_quad=8,16,32,64): {results}")
    diffs = [abs(results[i+1] - results[i]) for i in range(3)]
    print(f"    successive differences: {diffs}")
    # Differences should generally decrease (with possible plateauing once
    # quadrature error falls below other numerical noise)
    assert diffs[1] <= diffs[0] * 1.5, f"diff[1]={diffs[1]} should be < diff[0]={diffs[0]}"
    assert diffs[2] <= diffs[1] * 1.5, f"diff[2]={diffs[2]} should be < diff[1]={diffs[1]}"
    r_ref = exact_ss_toa([ray, hg], 0.5, 0.6, np.pi/4, 0.3, 1.0,
                         N_quad=128, N_phi=256)
    rel_err = abs(results[-1] - r_ref['total']) / abs(r_ref['total'])
    print(f"    rel error vs N_quad=128: {rel_err:.2e}")
    assert rel_err < 1e-4, f"convergence error {rel_err} should be < 1e-4"
    print("  quadrature convergence: OK")


def test_thin_limit_linearity():
    """In the optically thin limit (τ<<1), SS contributions are linear in τ.
    Use truly small base τ so exponentials are near 1."""
    base_tau = 0.001
    ray_unit = RayleighContributor([base_tau])
    r_unit = exact_ss_toa([ray_unit], 0.6, 0.5, 0.0, 0.0, 1.0)  # black surface
    base = r_unit['total']
    print(f"  thin-limit linearity check (base τ={base_tau}):")
    for scale in [0.5, 0.25, 0.1]:
        ray_scaled = RayleighContributor([base_tau * scale])
        r = exact_ss_toa([ray_scaled], 0.6, 0.5, 0.0, 0.0, 1.0)
        ratio = r['total'] / base
        print(f"    scale={scale}: ratio={ratio:.6f} (expected ≈ {scale})")
        assert_close(ratio, scale, rtol=0.01,
                     msg=f"linearity at scale={scale}")
    print("  thin-limit linearity: OK")


def test_no_atmosphere_total():
    empty = RayleighContributor([0.0])
    mu0, mu_v, albedo = 0.5, 0.3, 0.6
    r = exact_ss_toa([empty], mu0, mu_v, 0.4, albedo, 1.0)
    expected = mu0 * 1.0 * albedo / np.pi
    assert_close(r['total'], expected, msg="no-atm total")
    print(f"  no atmosphere reflectance limit: {r['total']:.6e}  OK")


def test_energy_budget_rayleigh_thin():
    """Hemispherical reflectance for thin Rayleigh + black surface should be small and ≈linear in τ."""
    mu0 = 0.7
    I0 = 1.0
    n_q = 16
    mu_nodes, mu_weights = gauss_legendre_01(n_q)
    n_p = 16

    def hemispheric_R(tau):
        ray = RayleighContributor([tau])
        R = 0.0
        for k in range(n_q):
            mu_v = mu_nodes[k]
            I_avg = 0.0
            for j in range(n_p):
                dphi = 2*np.pi * j / n_p
                r = exact_ss_toa([ray], mu0, mu_v, dphi, 0.0, I0,
                                 N_quad=8, N_phi=32)
                I_avg += r['total']
            I_avg /= n_p
            R += mu_weights[k] * I_avg * mu_v * 2*np.pi
        return R / (mu0 * I0)

    R_full = hemispheric_R(0.01)
    R_half = hemispheric_R(0.005)

    print(f"  hemispheric R (Rayleigh black surface):")
    print(f"    τ=0.010: R={R_full:.6f}")
    print(f"    τ=0.005: R={R_half:.6f}")
    print(f"    ratio:  {R_half/R_full:.4f} (expected ≈ 0.5 in thin limit)")
    assert 0 < R_full < 1, f"R_full = {R_full} out of bounds"
    assert_close(R_half / R_full, 0.5, rtol=0.02,
                 msg="hemispheric R linearity")
    print("  energy budget thin-Rayleigh: OK")


def main():
    print("Running exact-SS reference validation tests...")
    print()
    test_geometry()
    test_phase_normalization()
    test_rayleigh_azimuthal_average()
    test_gauss_legendre()
    test_path1_zenith_rayleigh()
    test_path1_offzenith_isotropic()
    test_path1_layer_subdivision()
    test_mixed_contributors()
    test_path2_no_atmosphere()
    test_black_surface()
    test_absorption_only_atm()
    test_no_atmosphere_total()
    test_thin_limit_linearity()
    test_quadrature_convergence()
    test_path1_brdf_reciprocity()
    test_path2_normalized_symmetry()
    test_path3_path4_brdf_reciprocity()
    test_energy_budget_rayleigh_thin()
    print()
    print("All validation tests passed.")


if __name__ == "__main__":
    main()
