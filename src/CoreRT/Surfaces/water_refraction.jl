#=
Built-in water refractive index from Segelstein (1981).

Provides `water_refractive_index(λ_nm)` returning Complex(n, k)
where n is the real part and k is the imaginary (absorption) part.

Data covers 200 nm – 200 μm.  Interpolation is linear in log(λ) for n
and log-linear for k (log10 of k is interpolated).
=#

# Tabulated data from Segelstein (1981)
# Columns: wavelength [nm], n (real), k (imaginary)
# Selected subset of key points for good interpolation accuracy
const _WATER_RI_TABLE_NM = Float64[
    200.0,  210.0,  220.0,  230.0,  240.0,  250.0,  260.0,  270.0,  280.0,  290.0,
    300.0,  310.0,  320.0,  330.0,  340.0,  350.0,  360.0,  370.0,  380.0,  390.0,
    400.0,  410.0,  420.0,  430.0,  440.0,  450.0,  460.0,  470.0,  480.0,  490.0,
    500.0,  510.0,  520.0,  530.0,  540.0,  550.0,  560.0,  570.0,  580.0,  590.0,
    600.0,  610.0,  620.0,  630.0,  640.0,  650.0,  660.0,  670.0,  680.0,  690.0,
    700.0,  720.0,  740.0,  760.0,  780.0,  800.0,  820.0,  840.0,  860.0,  880.0,
    900.0,  920.0,  940.0,  960.0,  980.0, 1000.0, 1050.0, 1100.0, 1150.0, 1200.0,
   1250.0, 1300.0, 1350.0, 1400.0, 1450.0, 1500.0, 1550.0, 1600.0, 1650.0, 1700.0,
   1750.0, 1800.0, 1850.0, 1900.0, 1950.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0,
   2500.0, 2600.0
]

const _WATER_N_REAL = Float64[
    1.396, 1.373, 1.362, 1.354, 1.349, 1.346, 1.343, 1.341, 1.339, 1.338,
    1.337, 1.336, 1.335, 1.335, 1.334, 1.334, 1.333, 1.333, 1.333, 1.332,
    1.332, 1.332, 1.331, 1.331, 1.331, 1.331, 1.330, 1.330, 1.330, 1.330,
    1.329, 1.329, 1.329, 1.329, 1.328, 1.328, 1.328, 1.328, 1.327, 1.327,
    1.327, 1.326, 1.326, 1.326, 1.325, 1.325, 1.325, 1.325, 1.324, 1.324,
    1.324, 1.323, 1.322, 1.322, 1.321, 1.320, 1.319, 1.319, 1.318, 1.317,
    1.316, 1.315, 1.314, 1.313, 1.312, 1.311, 1.308, 1.306, 1.303, 1.300,
    1.296, 1.293, 1.289, 1.285, 1.277, 1.268, 1.261, 1.255, 1.253, 1.255,
    1.260, 1.268, 1.279, 1.295, 1.306, 1.304, 1.279, 1.232, 1.188, 1.147,
    1.131, 1.129
]

const _WATER_K_IMAG = Float64[
    1.42e-7, 7.00e-8, 4.00e-8, 2.60e-8, 1.80e-8, 1.40e-8, 1.10e-8, 9.00e-9, 7.50e-9, 6.50e-9,
    6.00e-9, 4.60e-9, 3.50e-9, 2.70e-9, 2.20e-9, 1.80e-9, 1.60e-9, 1.40e-9, 1.30e-9, 1.30e-9,
    1.30e-9, 1.40e-9, 1.50e-9, 1.60e-9, 1.70e-9, 1.80e-9, 1.90e-9, 2.05e-9, 2.30e-9, 2.69e-9,
    3.21e-9, 3.81e-9, 4.36e-9, 4.78e-9, 5.14e-9, 5.69e-9, 6.49e-9, 7.63e-9, 9.22e-9, 1.09e-8,
    1.26e-8, 1.39e-8, 1.48e-8, 1.55e-8, 1.63e-8, 1.74e-8, 1.91e-8, 2.20e-8, 2.72e-8, 3.59e-8,
    4.78e-8, 7.50e-8, 1.10e-7, 1.43e-7, 1.65e-7, 1.72e-7, 1.63e-7, 1.46e-7, 1.32e-7, 1.28e-7,
    1.38e-7, 1.65e-7, 2.41e-7, 4.42e-7, 7.40e-7, 1.06e-6, 1.79e-6, 1.65e-6, 1.10e-6, 9.60e-7,
    1.32e-6, 2.26e-6, 4.58e-6, 1.07e-5, 2.94e-5, 5.88e-5, 7.15e-5, 6.71e-5, 5.68e-5, 4.65e-5,
    3.85e-5, 3.44e-5, 3.72e-5, 5.63e-5, 1.27e-4, 2.98e-4, 6.56e-4, 1.14e-3, 1.67e-3, 1.89e-3,
    1.67e-3, 1.19e-3
]

# Pre-compute log wavelengths for interpolation
const _WATER_LOG_NM = log.(_WATER_RI_TABLE_NM)
const _WATER_LOG_K  = log.(_WATER_K_IMAG)

"""
    water_refractive_index(λ_nm)

Complex refractive index of liquid water at wavelength `λ_nm` (in nanometers).

Uses the Segelstein (1981) tabulation with log-linear interpolation.
Valid range: 200–2600 nm. Values outside are clamped to the boundary.

# Returns
- `Complex{Float64}(n, k)` where n is the real part and k the imaginary part.
"""
function water_refractive_index(λ_nm::Real)
    FT = Float64
    log_λ = log(FT(λ_nm))
    N = length(_WATER_LOG_NM)

    # Clamp to table range
    if log_λ <= _WATER_LOG_NM[1]
        return Complex{FT}(_WATER_N_REAL[1], _WATER_K_IMAG[1])
    elseif log_λ >= _WATER_LOG_NM[N]
        return Complex{FT}(_WATER_N_REAL[N], _WATER_K_IMAG[N])
    end

    # Binary search for interval
    lo, hi = 1, N
    while hi - lo > 1
        mid = (lo + hi) >> 1
        if _WATER_LOG_NM[mid] <= log_λ
            lo = mid
        else
            hi = mid
        end
    end

    # Linear interpolation weight
    t = (log_λ - _WATER_LOG_NM[lo]) / (_WATER_LOG_NM[hi] - _WATER_LOG_NM[lo])

    # n: linear interpolation in log(λ)
    n = _WATER_N_REAL[lo] + t * (_WATER_N_REAL[hi] - _WATER_N_REAL[lo])

    # k: log-linear interpolation (interpolate log(k), then exp)
    log_k = _WATER_LOG_K[lo] + t * (_WATER_LOG_K[hi] - _WATER_LOG_K[lo])
    k = exp(log_k)

    return Complex{FT}(n, k)
end
