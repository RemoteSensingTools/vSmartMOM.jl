# Ordering in the supplied Fortran/C++ LBL implementation

Inspected 2026-09-06, read-only, at
`/home/sanghavi/data/vSmartMOM1.6_lbl_OCO2`. Source hashes are recorded in
`evidence/fortran_source_manifest.json`. This is a static examination of the
provided snapshot, not a claim that its present source matches a historical
executable or the exact implementation timed in either Jacobian paper.

## Which files form the selected program

The uppercase `Makefile` is empty. The lowercase `makefile`, lines 35–60,
selects `misr_image_4_new3.f90`, `rt_atmo_new2_homog.f90`,
`rt_elemental.f90`, `rt_doubling.f90`, `rt_interaction.f90`, and
`rt_atmo_surf3.f90`, alongside `atmo_prop_new6.f90` and
`phase_function_new.f90`. It links a C++ retrieval/optics driver including
`aerosol.cpp`. Several similarly named older alternatives are present;
in particular, `rt_atmo_surf3_lbl.f90` is not the selected surface routine.
The CMake file concerns the separate levmar library.

`parameter_block.f90`, lines 15–29, sets `SFI=.false.`, four Stokes components,
four aerosol single-scattering variables, four microphysical parameters and
two vertical-profile parameters per aerosol. These settings matter when
comparing its ordering and cost with modern Julia SFI runs.

## Wavelengths are sequential, with substantial reuse

The active driver has this order:

```text
band                                      misr_image_4_new3.f90:1111
  prepare aerosol optics / derivative coefficients for this band
  detect duplicate spectral bin profiles and common atmospheric prefixes
  Fourier moment                          :1491
    compute phase matrices                :1500
    atmospheric layer                     :1527
      distinct absorption bin             :1539
        ATMO_PROP                         :1555
        RT_ATMO_homog                     :1571
          ELEMENTAL                       rt_atmo_new2_homog.f90:466
          DOUBLING                        :485
        cache doubled operators and tangents
    spectral case jctr                    misr_image_4_new3.f90:1820
      load cached atmospheric prefix if available
      layer                               :1867
        select cached bin; transform local gas derivatives
        INTERACTION                       :2102
      RT_SURF / RT_ATMO_SURF               :2201 / :2205
      accumulate this Fourier contribution
```

The third/fourth dimensions of the cached arrays do not have the modern Julia
meaning of one wavelength-batched kernel. Each elemental/doubling call processes
one layer/bin matrix. The driver then assembles spectral cases in a serial loop.
However, it does **not** redo elemental and doubling work at every wavelength:
`inumbins` and `ilambdabin` let several wavelengths reuse a doubled layer.

`duplicate_wl`, `duplicate_z` and `write_z` also detect repeated complete bin
profiles or prefixes (driver lines 1386–1430). `read_data` at 1846 resumes from
a cached composite and `write_data` at 2128 saves useful prefixes. Thus a runtime
comparison must account for how many *distinct bins and prefixes* are solved,
not just the number of wavelengths. Bin discretization and any associated
spectral approximation are separate from the Jacobian factorization question.

## The Jacobian contraction order is hybrid

### Aerosols: physical parameters enter before doubling

`parameter_block.f90:22` describes the four core blocks as τ, ω, P and ftrunc.
But these are intermediate elemental partials, not the final doubling basis.

1. `rt_elemental.f90:331–355` constructs elemental derivatives with respect to
   those single-scattering blocks, including the local phase derivative.
2. At 358–379 it contracts them into each of the four microphysical directions,
   with `mp_drtv_tau`, `mp_drtv_w`, `cdrtv_mp_mp/pp`, and `mp_drtv_f`.
3. At 382–399 it constructs the vertical-profile tangents.
4. At 590–616 it packs AOT, microphysical and profile derivatives into
   `cdrtv_{r,t}_elt_*`. The AOT derivative is converted to reference-wavelength
   AOT with `k_ratio`.
5. `rt_doubling.f90:270–297` loops over `tot_aer_var` and propagates each of
   those matrices through every doubling step.

The driver defines `tot_aer_var = NAer + tot_mp + tot_vp` at 476–481.
With the checked constants this is **7 NAer**, not four universal core
columns. Aerosol derivatives are already expanded before the call to DOUBLING.
The forward inverse and several forward products are shared across those
columns (`rt_doubling.f90:223–233`).

### Rayleigh/absorption: local optical directions are contracted after doubling

This is the concrete delayed-contraction precedent in the supplied source.
`tau_comp_ray` includes molecular scattering **and gas absorption**; despite
the field name it is not solely Rayleigh scattering optical depth. The driver
forms `tau_comp_ray = tau_scat + tau_abs` and
`w0_ray = tau_scat / tau_comp_ray` at 1236–1252.

`atmo_prop_new6.f90:234–245` differentiates the mixture with respect to these
two local optical quantities. `rt_elemental.f90:296–304` forms the resulting
two operator tangents, and `rt_doubling.f90:235–263` doubles them. They are
cached per layer/bin (`misr_image_4_new3.f90:1704–1717`).

Only when loading a doubled layer for adding does the driver apply

```math
 L_{p_0}=L_{\tau_r}\,\tau_{r,p_0}+L_{\omega_r}\,\omega_{r,p_0},\qquad
 L_{v_{CO_2}}=L_{\tau_r}\,\tau_{r,v_{CO_2}}
             +L_{\omega_r}\,\omega_{r,v_{CO_2}}.
```

The first-layer mapping is at 1933–1956; subsequent layers use 1968–1994.
The scalar coefficient arrays are constructed at 1258–1265. The two resulting
columns represent pressure and a **single CO₂ VMR**, not an independent
multi-gas profile state with a column for every layer. Adding propagates those
two columns (`rt_interaction.f90:493–527`) and the `7 NAer` aerosol columns
(531 onward).

For one aerosol this snapshot therefore doubles seven aerosol directions plus
two molecular optical directions at m≤2. It is not a comparable 60-column
profile retrieval. The two molecular tangent families are also gated at m≤2
in several places. That historical gate must not be interpreted as proof that
gas absorption has zero effect on higher-order aerosol Fourier radiances.

## Truncation and source handling

The Fortran boundary explicitly receives `in_drtv_trunc_f` and Greek tangent
arrays, placing them in `mp_drtv_f` and `mp_drtv_l*` at driver lines 1339–1354.
The elemental chain rule includes the f derivative, as shown above.

The selected C++ optics code includes δ-M and fitted-truncation alternatives.
The active normalized-Greek equations in `delta_fit_truncation2`
(`aerosol.cpp:2148–2161`) have the expected form
`d(g_fit/c0) = (dg_fit − g_truncated dc0)/c0`, with `c0=1−f`.
This is the same normalization principle checked independently in Julia.

However, the supplied C++ fitted-truncation implementation has apparent defects:

- At `aerosol.cpp:1940`, the degree loop assigns `gsl_vector_get(x,0)` to
  **every** `drtv_c_beta[i][l]`, instead of selecting component l.
- Within that function the only assignment to the output `drtv_truncf` is in
  the commented-out old linearization block. The active replacement computes
  `drtv_c_beta` but does not write `drtv_truncf = -dc0` to the caller's vector.

These are static source findings, not reproduced failures from a rebuilt
executable. They make this snapshot useful for architecture but unsuitable as
an unquestioned numerical reference. The original directory was not modified.

With `SFI=false`, this build obtains solar radiance through the embedded solar
operator representation. It does not establish correctness or performance of
the modern Julia external-solar source tangent path. Dormant SFI code does move
the above-layer forward beam attenuation to the cached-layer loading step
(driver 1891–1904); a complete corresponding source-Jacobian implementation
must still be independently validated. The Julia delayed-contraction experiment
already checks the differentiated attenuation term explicitly.

## Consequences for the Julia work

The user's direction is to take useful ideas while retaining everything the
modern code already does better. The current Julia solver is the baseline;
this source inspection does not change its equations or numerical behavior.

| Decision | Reason |
|---|---|
| Carry forward local optical tangents with delayed retrieval contraction | The molecular path is a concrete precedent; the modern 96-assertion experiment independently verifies the general doubling identity and source attenuation. |
| Keep GPU wavelength batching, fused pivoted solves and reusable scratch | These are measured improvements in the modern implementation; sequential old Fortran loops are not the target architecture. |
| Keep finite-thickness elemental formulas, current polarization conventions and validated embedded/external-solar source paths | Historical settings and dormant source code do not supersede the current math and tests. |
| Defer spectral-bin/prefix reuse as a separate workload-dependent experiment | Binning changes the representation of spectral optical states and could add overhead or approximation; it is not needed to implement the local Jacobian basis. |
| Leave historical derivative defects and fixed two-column gas layout behind | A modern multi-gas profile retrieval requires its actual parameter layout and finite-difference validation. |

The source supports generalizing **local optical propagation followed by
retrieval contraction**, because it already does that for molecular optics.
It does not show the proposed complete aerosol phase-basis method already
implemented: aerosol microphysical/profile columns are expanded early.
The next Julia design can combine a local molecular basis with a factored
truncated aerosol phase basis, preserving wavelength batching on the GPU.

Keep three improvements separate when measuring gains:

- compact local tangent bases versus retrieval-column propagation;
- wavelength batching and kernel fusion;
- reuse across repeated spectral optical states / layer prefixes.

Finally, the Fortran forward cost structure differs even apart from these
choices: `.x.` uses intrinsic `MATMUL` (`matrix_multiply.f90:43`), and the
selected `INVERT_MATRIX2` (`invert_matrix.f90:236–327`) performs handwritten
augmented-matrix elimination rather than modern batched pivoted LU.
No runtime ratio was measured for the supplied code in this examination.
Its quadrature builder also explicitly deduplicates camera zenith cosines
(`spl_quad_pts.f90:102–115`) before adding nodes (163–164), reinforcing that
301 viewing directions do not automatically imply 301 distinct angular nodes.
