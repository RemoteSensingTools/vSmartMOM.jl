# # Surface BRDF Models: Theory and Usage
#
# vSmartMOM supports a family of Bidirectional Reflectance Distribution Function
# (BRDF) models for the lower boundary condition in radiative transfer.  This
# tutorial covers the mathematical background, code usage, and guidance for every
# available surface type -- from the simplest isotropic Lambertian to the fully
# polarized Cox-Munk ocean surface.
#
# **Contents**
# 1. Overview and model comparison
# 2. Lambertian surfaces
# 3. RPV (Rahman-Pinty-Verstraete) model
# 4. Ross-Li kernel-based model
# 5. Cox-Munk polarized ocean surface
# 6. Canopy-coupled surfaces
# 7. Polarization and Mueller matrices
# 8. Surface Jacobians (linearized RT)
# 9. Choosing the right model
# 10. Running a full example
# 11. YAML configuration
# 12. References
#
# ---

using vSmartMOM
using vSmartMOM.CoreRT
using CairoMakie

# ## 1) Overview
#
# Every surface model is a subtype of `AbstractSurfaceType`.  The RT solver
# dispatches through `create_surface_layer!`, which computes the
# Fourier-decomposed surface reflectance matrix $R_\text{surf}$ for each
# azimuthal moment $m$.
#
# | Model | Constructor | Params | Polarized? | Jacobians? | Use case |
# |:------|:-----------|:------:|:----------:|:----------:|:---------|
# | Lambertian (scalar) | `LambertianSurfaceScalar(α)` | 1 | no | yes | Simple land, snow, dark ocean |
# | Lambertian (spectral) | `LambertianSurfaceSpectrum([α₁,…])` | N | no | no | Wavelength-varying surfaces |
# | RPV | `rpvSurfaceScalar(ρ₀, ρ_c, k, Θ)` | 4 | no | no | Vegetation, soils |
# | Ross-Li | `RossLiSurfaceScalar(fvol, fgeo, fiso)` | 3 | no | no | MODIS-style vegetation |
# | Cox-Munk | `CoxMunkSurface(wind_speed=U,…)` | 1-5 | **yes** | yes | Ocean, sun glint |
# | Canopy | `CanopySurface(soil=…, LAI=…,…)` | many | via soil | partial | Vegetation canopy + soil |
#
# Two additional Lambertian variants (`LambertianSurfaceLegendre`,
# `LambertianSurfaceSpline`) are available internally for polynomial or
# spline-based spectral albedo but are not exported.
#
# ---

# ## 2) Lambertian Surface
#
# The simplest model: reflectance is independent of viewing and illumination
# angles.  A perfect Lambertian surface satisfies
#
# $$\rho(\mu_i, \mu_r, \phi) \;=\; \frac{\alpha}{\pi}$$
#
# where $\alpha \in [0,1]$ is the hemispherical albedo.  Because the BRDF has
# no azimuthal dependence, only the $m = 0$ Fourier moment is nonzero.  The
# code applies a factor of 2 to $m = 0$ for consistency with the Fourier
# series normalization used in the RT solver.
#
# ### Scalar albedo
#
# Use `LambertianSurfaceScalar` for a single albedo across all wavelengths in
# the band:

lam_scalar = LambertianSurfaceScalar(0.15)

# ### Spectral albedo
#
# Use `LambertianSurfaceSpectrum` when the albedo varies per spectral point
# (one value per grid point in the band):

lam_spectrum = LambertianSurfaceSpectrum([0.12, 0.13, 0.14, 0.15])

# ### Legendre and spline variants
#
# Not exported, but accessible from `CoreRT`:
#
# - `CoreRT.LambertianSurfaceLegendre([c₀, c₁, …])`: polynomial spectral
#   albedo $\alpha(\nu) = c_0 P_0(\hat\nu) + c_1 P_1(\hat\nu) + \cdots$,
#   where $\hat\nu$ is the spectral grid rescaled to $[-1,1]$.
# - `CoreRT.LambertianSurfaceSpline(interpolator, wlGrid)`: spline-interpolated
#   albedo from a user-supplied `Interpolations.jl` object.
#
# ---

# ## 3) RPV (Rahman-Pinty-Verstraete) Model
#
# The RPV model (Rahman, Pinty & Verstraete, 1993) provides a four-parameter
# empirical BRDF widely used for vegetation and soils:
#
# $$\rho \;=\; \rho_0 \;\cdot\; M(\mu_i, \mu_r, k) \;\cdot\; F(\Theta, \cos g) \;\cdot\; H(\rho_c, G)$$
#
# where the three factors are:
#
# **Minnaert function** (angular shape):
# $$M(\mu_i, \mu_r, k) \;=\; \frac{(\mu_i\, \mu_r)^{k-1}}{(\mu_i + \mu_r)^{1-k}}$$
#
# **Hot-spot function** (Henyey-Greenstein-like):
# $$F(\Theta, \cos g) \;=\; \frac{1 - \Theta^2}{\bigl(1 + \Theta^2 + 2\Theta \cos g\bigr)^{3/2}}$$
#
# **Geometric bowl function**:
# $$H(\rho_c, G) \;=\; 1 + \frac{1 - \rho_c}{1 + G}$$
#
# with the scattering angle cosine
# $\cos g = -\mu_i \mu_r + \sin\theta_i \sin\theta_r \cos\phi$
# and the geometric factor
# $G = [\tan^2\theta_i + \tan^2\theta_r + 2\tan\theta_i \tan\theta_r \cos\phi]^{1/2}$.
#
# | Parameter | Symbol | Typical range | Physical meaning |
# |:----------|:------:|:-------------:|:-----------------|
# | Overall amplitude | $\rho_0$ | 0 - 0.5 | Isotropic reflectance level |
# | Minnaert exponent | $k$ | 0.3 - 1.5 | $k < 1$: bowl shape, $k > 1$: bell shape |
# | Hot-spot asymmetry | $\Theta$ | $-1$ to $1$ | $\Theta < 0$: backscatter, $\Theta > 0$: forward |
# | Geometric amplitude | $\rho_c$ | $-1$ to $1$ | 1 = no geometric term |
#
# RPV is scalar-only: it returns zero for Stokes components beyond $I$.

rpv = rpvSurfaceScalar(0.15, 0.1, 0.7, -0.3)  ## ρ₀, ρ_c, k, Θ

# ---

# ## 4) Ross-Li Model
#
# The Ross-Li model (Lucht, Schaaf & Strahler, 2000) decomposes the BRDF into
# a linear combination of semi-physical kernels:
#
# $$\rho \;=\; f_\text{iso}\, K_\text{iso} \;+\; f_\text{vol}\, K_\text{vol}(\theta_i, \theta_r, \phi) \;+\; f_\text{geo}\, K_\text{geo}(\theta_i, \theta_r, \phi)$$
#
# **Isotropic kernel**: $K_\text{iso} = 1$.
#
# **Ross Thick volumetric kernel** (dense canopy scattering):
# $$K_\text{vol} = \frac{(\pi/2 - \xi)\cos\xi + \sin\xi}{\cos\theta_i + \cos\theta_r} - \frac{\pi}{4}$$
# where $\xi = \arccos\!\bigl(\cos\theta_i\cos\theta_r + \sin\theta_i\sin\theta_r\cos\phi\bigr)$
# is the scattering angle.
#
# **Li Sparse geometric kernel** (shadowing by sparsely distributed objects):
# $$K_\text{geo} = O(\theta_i', \theta_r', \phi) - \sec\theta_i' - \sec\theta_r' + \tfrac{1}{2}(1 + \cos\xi')\sec\theta_i'\sec\theta_r'$$
# where $O$ is the overlap function and primed angles are scaled by the
# crown shape ratio ($h/b = 2$, $b/r = 1$ RAMI defaults).
#
# | Parameter | Symbol | Typical range | Physical meaning |
# |:----------|:------:|:-------------:|:-----------------|
# | Isotropic | $f_\text{iso}$ | 0 - 0.3 | Baseline reflectance |
# | Volumetric | $f_\text{vol}$ | 0 - 0.1 | Dense-canopy scattering |
# | Geometric | $f_\text{geo}$ | 0 - 0.1 | Mutual shadowing |
#
# Ross-Li is scalar-only: it returns zero for Stokes components beyond $I$.

rossli = RossLiSurfaceScalar(0.05, 0.03, 0.1)  ## fvol, fgeo, fiso

# ---

# ## 5) Cox-Munk Polarized Ocean Surface
#
# `CoxMunkSurface` is the first **fully polarized** surface model in vSmartMOM.
# It models the ocean as an ensemble of specular wave facets whose tilt
# distribution follows a wind-dependent Gaussian.  Fresnel reflection from
# each facet produces a full Mueller-matrix BRDF that couples all four Stokes
# components.  This is critical for remote sensing of ocean scenes (OCO-2/3,
# PACE, etc.) where sun glint polarization is a dominant signal.
#
# ### 5a) Slope distribution (Cox & Munk, 1954)
#
# The ocean surface is approximated as a collection of small flat facets
# ("capillary waves") whose slopes $(z_x, z_y)$ follow an isotropic Gaussian:
#
# $$P(z_x, z_y) \;=\; \frac{1}{2\pi\sigma^2}\,\exp\!\left(-\frac{z_x^2 + z_y^2}{2\sigma^2}\right)$$
#
# where the slope variance depends linearly on 10-m wind speed $U$ (m/s):
#
# $$\sigma^2 \;=\; 0.003 + 0.00512\, U$$
#
# Low wind ($U \lesssim 2$ m/s) gives a nearly mirror-like surface with a
# narrow specular glint; high wind ($U \gtrsim 10$ m/s) spreads the glint
# over a wide angular range.
#
# ### 5b) Fresnel reflection
#
# Each facet reflects specularly according to Fresnel's equations.  For a
# complex relative refractive index $n = n_r + i\,n_i$ (water/air) and local
# incidence angle $\theta_\text{local}$:
#
# $$r_s = \frac{\cos\theta_i - n\cos\theta_t}{\cos\theta_i + n\cos\theta_t}, \qquad r_p = \frac{n\cos\theta_i - \cos\theta_t}{n\cos\theta_i + \cos\theta_t}$$
#
# where $\cos\theta_t = \sqrt{1 - \sin^2\theta_i / n^2}$ (Snell's law).
#
# The Fresnel **Mueller matrix** for reflection is:
#
# $$\mathbf{M}_F = \begin{pmatrix} \tfrac{|r_s|^2+|r_p|^2}{2} & \tfrac{|r_s|^2-|r_p|^2}{2} & 0 & 0 \\ \tfrac{|r_s|^2-|r_p|^2}{2} & \tfrac{|r_s|^2+|r_p|^2}{2} & 0 & 0 \\ 0 & 0 & \operatorname{Re}(r_s r_p^*) & \operatorname{Im}(r_s r_p^*) \\ 0 & 0 & -\operatorname{Im}(r_s r_p^*) & \operatorname{Re}(r_s r_p^*) \end{pmatrix}$$
#
# The off-diagonal (1,2) and (2,1) elements couple the $I$ and $Q$ Stokes
# parameters -- this is the polarizing effect of the ocean surface.
#
# ### 5c) BRDF Mueller matrix
#
# For a given geometry $(\mu_i, \mu_r, \Delta\phi)$ the full BRDF is:
#
# $$\boldsymbol{\rho} \;=\; \frac{P(z_x, z_y)\;\cdot\; S(\mu_i, \mu_r)}{4\,\mu_i\,\mu_r\,\cos^4\!\beta}\;\; \mathbf{L}(\alpha_2)\;\mathbf{M}_F(\cos\theta_\text{local})\;\mathbf{L}(-\alpha_1)$$
#
# where:
# - $\beta$ is the facet tilt angle ($\cos\beta = n_z$, the $z$-component of the facet normal)
# - $\theta_\text{local}$ is the local incidence angle on the tilted facet
# - $\alpha_1$, $\alpha_2$ are rotation angles between the scattering plane
#   and the facet incidence/reflection planes
# - $\mathbf{L}(\phi)$ is the Stokes rotation matrix:
#
# $$\mathbf{L}(\phi) = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos 2\phi & \sin 2\phi & 0 \\ 0 & -\sin 2\phi & \cos 2\phi & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$
#
# - $S$ is the Smith (1967) bistatic shadow masking factor
#
# ### 5d) Shadow masking (Smith, 1967)
#
# At low grazing angles, parts of the surface are shadowed by neighboring wave
# crests.  The Smith shadowing function corrects for this:
#
# $$S(\mu_i, \mu_r) \;=\; \frac{1}{1 + \Lambda(\mu_i) + \Lambda(\mu_r)}$$
#
# where $\Lambda(\mu) = \tfrac{1}{2}\!\left[\tfrac{\exp(-\nu^2)}{\sqrt{2\pi}\,\nu} - \operatorname{erfc}(\nu)\right]$
# with $\nu = \cot\theta / (\sqrt{2}\,\sigma)$ and $\sigma = \sqrt{\sigma^2}$.
#
# ### 5e) Whitecap contribution (Monahan & O'Muircheartaigh, 1980)
#
# At higher wind speeds, breaking waves form whitecaps that scatter light
# diffusely.  The fractional whitecap coverage is:
#
# $$f_\text{wc} \;=\; 2.95 \times 10^{-6} \cdot U^{3.52}$$
#
# Whitecaps are modeled as unpolarized Lambertian reflectors with albedo
# $\alpha_\text{wc}$ (default 0.22, Koepke 1984).  The total BRDF is:
#
# $$\boldsymbol{\rho}_\text{total} \;=\; (1 - f_\text{wc})\,\boldsymbol{\rho}_\text{glint} \;+\; f_\text{wc}\,\boldsymbol{\rho}_\text{wc}$$
#
# where $\boldsymbol{\rho}_\text{wc}$ has only a $(1,1)$ element equal to
# $\alpha_\text{wc}/\pi$.
#
# ### 5f) Water refractive index
#
# When `n_water = nothing` (the default), the code uses a built-in lookup
# table from Segelstein (1981) covering 200 nm -- 2600 nm.  The function
# `water_refractive_index(λ_nm)` returns the complex refractive index at any
# wavelength via log-linear interpolation:

n_water_550 = water_refractive_index(550.0)
println("Water RI at 550 nm:  n = ", round(real(n_water_550), digits=4),
        ", k = ", round(imag(n_water_550), sigdigits=3))

# You can also provide a fixed value or a per-wavelength vector via the
# `n_water` keyword of `CoxMunkSurface`.
#
# ### 5g) TMS single-scattering correction
#
# The specular sun-glint peak can be extremely narrow (especially at low wind
# speed), requiring many Fourier moments to resolve.  Rather than increasing
# `max_m`, the code applies a **Truncated Multiple Scattering (TMS)**
# correction after the Fourier loop: it adds the difference between the exact
# single-scattering surface contribution and the truncated Fourier
# reconstruction at each viewing geometry.  This is handled automatically by
# `apply_ss_correction!` inside `rt_run`.
#
# ### 5h) Construction

ocean = CoxMunkSurface(wind_speed = 5.0)  ## uses built-in water RI, whitecaps on

ocean_custom = CoxMunkSurface(
    wind_speed      = 7.5,
    n_water         = complex(1.34, 1e-7),  ## override refractive index
    whitecap_albedo = 0.22,
    include_whitecaps = true,
    shadowing       = true,
)

# ---

# ## 6) Canopy-Coupled Surface
#
# `CanopySurface` represents a vegetation canopy backed by a soil BRDF.  It
# internally solves canopy sub-layers via the adding-doubling method before
# presenting an effective reflectance to the atmospheric RT.  This is covered
# in detail in the [Canopy Tutorial](Tutorial_Canopy.md); here is a brief
# example:

canopy = CanopySurface(
    soil = LambertianSurfaceScalar(0.1),
    LAI  = 3.0,
    n_layers = 1,
    leaf_reflectance   = 0.45,
    leaf_transmittance = 0.05,
)

# `CanopySurface_from_prospect(leaf, wl_grid; ...)` provides a convenience
# constructor that computes spectral leaf optics from the PROSPECT model.
# See the Canopy Tutorial for multi-layer setups, spectral leaf optics, and
# within-canopy atmospheric absorption.
#
# ---

# ## 7) Polarization and Mueller Matrices
#
# The Lambertian, RPV, and Ross-Li models are **scalar-only**: they populate
# only the $(1,1)$ Stokes block of the reflectance matrix and return zero for
# all off-diagonal blocks.  This is sufficient for intensity-only ($I$) RT.
#
# The Cox-Munk surface is the first model to fill the **full Mueller matrix**
# at the surface boundary.  The Fresnel reflection off tilted wave facets
# couples the $I$-$Q$ and $U$-$V$ Stokes components.  Concretely:
#
# - For `Stokes_I()`: the reflectance matrix is $[N_\mu \times N_\mu]$.
# - For `Stokes_IQUV()`: it becomes $[4 N_\mu \times 4 N_\mu]$ with
#   all Stokes cross-coupling blocks populated.
#
# The Fourier decomposition uses polarization-aware azimuthal kernels:
# $\cos(m\Delta\phi)$ for the $I$-$Q$ / $I$-$Q$ blocks and
# $\sin(m\Delta\phi)$ for the cross blocks, matching the postprocessing
# reconstruction convention.
#
# ---

# ## 8) Surface Jacobians (Linearized RT)
#
# Several surface types support analytical Jacobians for use with
# `rt_run_lin`:
#
# ### Lambertian (analytical)
# The derivative of the surface reflectance matrix w.r.t. albedo $\alpha$ is
# trivially:
# $$\frac{\partial R_\text{surf}}{\partial \alpha} \;=\; 2\,\mathbf{E}\,\mathbf{w}^T$$
# for $m = 0$, where $\mathbf{E}$ is the identity-like Stokes selector and
# $\mathbf{w}$ are the quadrature weights.
#
# ### Cox-Munk (analytical)
# The derivative of the surface reflectance w.r.t. wind speed $U$ is computed
# analytically via the chain rule:
# $$\frac{\partial \boldsymbol{\rho}}{\partial U} \;=\; (1 - f_\text{wc})\!\left[\frac{\partial P}{\partial \sigma^2}\,S + P\,\frac{\partial S}{\partial \sigma^2}\right]\!\frac{d\sigma^2}{dU}\;\frac{\mathbf{L}_2\,\mathbf{M}_F\,\mathbf{L}_1}{4\mu_i\mu_r\cos^4\!\beta} \;+\; \frac{df_\text{wc}}{dU}\left(\boldsymbol{\rho}_\text{wc} - \boldsymbol{\rho}_\text{glint}\right)$$
#
# where $d\sigma^2/dU = 0.00512$.  The geometry, Fresnel coefficients, and
# rotation matrices are all independent of $U$ and shared with the forward
# evaluation.
#
# ### Canopy
# Jacobians with respect to the soil albedo use finite-difference
# perturbation of the internal adding-doubling solve.
#
# ### RPV, Ross-Li
# No built-in linearization.  Use external finite differences or ForwardDiff.
#
# ---

# ## 9) Choosing the Right Model
#
# | Scenario | Recommended | Reason |
# |:---------|:------------|:-------|
# | Quick test, isotropic surface | `LambertianSurfaceScalar` | Fastest, simplest |
# | Ocean, sun glint, polarimetry | `CoxMunkSurface` | Full polarization, wind-speed Jacobians |
# | Vegetation (MODIS/RAMI) | `RossLiSurfaceScalar` | Standard kernel-based, 3 params |
# | Vegetation (empirical) | `rpvSurfaceScalar` | 4-param empirical with hotspot |
# | Physical canopy model | `CanopySurface` | Leaf-level scattering, spectral optics |
# | Spectral albedo from data | `LambertianSurfaceSpectrum` | Per-wavelength albedo vector |
# | Polarized atmosphere, scalar surface | `LambertianSurfaceScalar` | Only $(1,1)$ block filled |
#
# ---

# ## 10) Running a Full Example
#
# Load a parameter file, override the surface, and run the RT:

yaml_path = joinpath(pkgdir(vSmartMOM),
                     "test", "test_parameters", "PureRayleighParameters.yaml")
params = read_parameters(yaml_path)
params.brdf[1] = CoxMunkSurface(wind_speed = 5.0)
model = model_from_parameters(params)
R, T = rt_run(model)
println("R shape: ", size(R))
println("R(nadir, I): ", R[1, 1, 1])

# Compare with a Lambertian surface at the same geometry:
params2 = read_parameters(yaml_path)
params2.brdf[1] = LambertianSurfaceScalar(0.06)
model2 = model_from_parameters(params2)
R2, T2 = rt_run(model2)
println("R(nadir, I) Lambertian 0.06: ", R2[1, 1, 1])

# Compare the reflectance spectra from the two surface types:

fig = Figure(size=(700, 450))
ax = Axis(fig[1,1],
    xlabel = "Spectral index",
    ylabel = "TOA Reflectance (Stokes I)")
lines!(ax, R[1, 1, :],  label="Cox-Munk (U=5 m/s)")
lines!(ax, R2[1, 1, :], label="Lambertian (α=0.06)")
axislegend(ax, position=:rt)
fig

# ---

# ## 11) YAML Configuration
#
# In a YAML parameter file, the surface is specified as a constructor string
# under `radiative_transfer.surface` (one entry per spectral band):
#
# ```yaml
# radiative_transfer:
#   surface:
#     - LambertianSurfaceScalar(0.15)
#     - rpvSurfaceScalar(0.12, 0.08, 0.75, -0.25)
#     - RossLiSurfaceScalar(0.04, 0.02, 0.08)
#     - CoxMunkSurface(wind_speed=5.0)
# ```
#
# For `CoxMunkSurface`, optional keyword arguments can be specified inline:
# ```yaml
#     - CoxMunkSurface(wind_speed=7.5, include_whitecaps=true, shadowing=true)
# ```
#
# For canopy surfaces, use the dedicated `canopy` YAML section (see the
# Canopy Tutorial), which wraps the `surface` entry as the soil BRDF.
#
# ---

# ## 12) References
#
# - Cox, C. & Munk, W. (1954). Measurement of the roughness of the sea surface from photographs of the sun's glitter. *JOSA* 44(11), 838-850. [doi:10.1364/JOSA.44.000838](https://doi.org/10.1364/JOSA.44.000838)
# - Koepke, P. (1984). Effective reflectance of oceanic whitecaps. *Applied Optics* 23(11), 1816-1824. [doi:10.1364/AO.23.001816](https://doi.org/10.1364/AO.23.001816)
# - Lucht, W., Schaaf, C.B. & Strahler, A.H. (2000). An algorithm for the retrieval of albedo from space using semiempirical BRDF models. *IEEE Trans. Geosci. Remote Sens.* 38(2), 977-998. [doi:10.1109/36.841980](https://doi.org/10.1109/36.841980)
# - Mishchenko, M.I. & Travis, L.D. (1997). Satellite retrieval of aerosol properties over the ocean using polarization as well as intensity of reflected sunlight. *JGR* 102(D14), 16989-17013. [doi:10.1029/97JD01084](https://doi.org/10.1029/97JD01084)
# - Monahan, E.C. & O'Muircheartaigh, I.G. (1980). Optimal power-law description of oceanic whitecap coverage dependence on wind speed. *J. Phys. Oceanogr.* 10, 2094-2099.
# - Rahman, H., Pinty, B. & Verstraete, M.M. (1993). Coupled surface-atmosphere reflectance (CSAR) model: 2. Semiempirical surface model usable with NOAA advanced very high resolution radiometer data. *JGR* 98(D11), 20791-20801. [doi:10.1029/93JD00564](https://doi.org/10.1029/93JD00564)
# - Segelstein, D.J. (1981). The complex refractive index of water. M.S. Thesis, University of Missouri-Kansas City.
# - Smith, B.G. (1967). Geometrical shadowing of a random rough surface. *IEEE Trans. Antennas Propag.* 15(5), 668-671. [doi:10.1109/TAP.1967.1138991](https://doi.org/10.1109/TAP.1967.1138991)
# - Zhai, P.-W., Hu, Y., Chowdhary, J., Trepte, C.R., Lucker, P.L. & Josset, D.B. (2010). A vector radiative transfer model for coupled atmosphere and ocean systems with a rough interface. *JQSRT* 111(7-8), 1025-1040. [doi:10.1016/j.jqsrt.2009.12.005](https://doi.org/10.1016/j.jqsrt.2009.12.005)
