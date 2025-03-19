

"""
    $(FUNCTIONNAME)(RS_type, λ, grid_in)

Returns the Raman SSA per layer at reference wavelength `λ` from nearby source wavelengths for RRS and a (currently) single incident wavelength for VRS/RVRS
(N₂,O₂ atmosphere, i.e. terrestrial)

Input: 
    - `RS_type` Raman scattering type (RRS/RVRS/VRS)
    - `λ` wavelength in `[μm]`
    - `grid_in` wavenumber grid with equidistant gridpoints
"""
function getRamanSSProp(RS_type::noRS, depol,  λ, grid_in)

    return nothing
end

function getRamanSSProp!(RS_type::Union{VS_0to1,VS_1to0}, depol, λ, grid_in)
    @unpack n2,o2 =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    λ_scatt = 2.e7/(grid_in[1]+grid_in[end])     
    #determine Rayleigh scattering cross-section at mean scattered wavelength
    atmo_σ_Rayl_scatt = compute_optical_Rayl!(λ_scatt, n2, o2)
    # determine Rayleigh scattering cross-section at incident wavelength λ 
    atmo_σ_Rayl = compute_optical_Rayl!(λ, n2, o2)
    # determine RRS cross-sections to λ₀ from nSpecRaman wavelengths around λ₀  
    index_VRSgrid_out, atmo_σ_VRS, index_RVRSgrid_out, atmo_σ_RVRS = 
        compute_optical_RS!(RS_type, grid_in, λ, n2, o2)
    # declare ϖ_Raman to be a grid of length raman grid
    RS_type.ϖ_VRS = atmo_σ_VRS/atmo_σ_Rayl
    RS_type.i_VRS = index_VRSgrid_out
    RS_type.ϖ_RVRS = atmo_σ_RVRS/atmo_σ_Rayl
    RS_type.i_RVRS = index_RVRSgrid_out
    RS_type.k_Rayl_scatt = atmo_σ_Rayl_scatt/atmo_σ_Rayl
    RS_type.n_Raman = 1;
    return nothing
end
#function getRamanLayerSSA(RS_type::VRS_1to0, T, λ, grid_in)
#    @unpack n2,o2 =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    # determine Rayleigh scattering cross-section at single monochromatic wavelength λ of the spectral band (assumed constant throughout the band)
#    atmo_σ_Rayl = compute_optical_Rayl!(λ, n2, o2)
    # determine RRS cross-sections to λ₀ from nSpecRaman wavelengths around λ₀  
#    index_VRSgrid_out, atmo_σ_VRS_1to0, index_RVRSgrid_out, atmo_σ_RVRS_1to0 = 
#        compute_optical_RS!(RS_type, grid_in, λ, n2, o2)
    # declare ϖ_Raman to be a grid of length raman grid
#    ϖ_VRS = atmo_σ_VRS_1to0/atmo_σ_Rayl
#    i_VRS = index_VRSgrid_out
#    ϖ_RVRS = atmo_σ_RVRS_1to0/atmo_σ_Rayl
#    i_RVRS = index_RVRSgrid_out
#    return ϖ_RVRS, i_RVRS, ϖ_VRS, i_VRS
#end



function getRamanSSProp!(RS_type::RRS, depol, λ, grid_in) 
    @unpack n2,o2 =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    # determine Rayleigh scattering cross-section at central wavelength λ of the spectral band (assumed constant throughout the band)
    atmo_σ_Rayl = compute_optical_Rayl(λ, n2, o2)
    RS_type.greek_raman = get_greek_raman(RS_type, n2, o2)
    #get_greek_raman!(RS_type, n2, o2)
    RS_type.ϖ_Cabannes .= compute_ϖ_Cabannes(RS_type, λ)
    # @show RS_type.ϖ_Cabannes
    # determine RRS cross-sections to λ₀ from nSpecRaman wavelengths around λ₀  
    index_raman_grid, atmo_σ_RRS = compute_optical_RS!(RS_type, grid_in, λ, n2, o2)
    # declare ϖ_Raman to be a grid of length raman grid
    #RS_type.ϖ_λ₁λ₀ = atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl * (1-RS_type.ϖ_Cabannes[1])/sum(atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl) #the grid gets inverted because the central wavelength is now seen as the recipient of RRS from neighboring source wavelengths
    RS_type.ϖ_λ₁λ₀ = (atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl) #the grid gets inverted because the central wavelength is now seen as the recipient of RRS from neighboring source wavelengths
    #@show RS_type.ϖ_λ₁λ₀
    #@show sum(RS_type.ϖ_λ₁λ₀)
    #@show RS_type.ϖ_Cabannes
    RS_type.i_λ₁λ₀ = index_raman_grid[end:-1:1]
    RS_type.n_Raman = length(RS_type.ϖ_λ₁λ₀)
    return nothing
end

function getRamanSSProp!(RS_type::RRS, λ, grid_in) 
    @unpack n2,o2 =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    # determine Rayleigh scattering cross-section at central wavelength λ of the spectral band (assumed constant throughout the band)
    atmo_σ_Rayl = compute_optical_Rayl(λ, n2, o2)
    RS_type.greek_raman = get_greek_raman(RS_type, n2, o2)
    #get_greek_raman!(RS_type, n2, o2)
    RS_type.ϖ_Cabannes .= compute_ϖ_Cabannes(RS_type, λ)
    # @show RS_type.ϖ_Cabannes
    # determine RRS cross-sections to λ₀ from nSpecRaman wavelengths around λ₀  
    index_raman_grid, atmo_σ_RRS = compute_optical_RS!(RS_type, grid_in, λ, n2, o2)
    # declare ϖ_Raman to be a grid of length raman grid
    #RS_type.ϖ_λ₁λ₀ = atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl * (1-RS_type.ϖ_Cabannes[1])/sum(atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl) #the grid gets inverted because the central wavelength is now seen as the recipient of RRS from neighboring source wavelengths
    RS_type.ϖ_λ₁λ₀ = (atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl) #the grid gets inverted because the central wavelength is now seen as the recipient of RRS from neighboring source wavelengths
    #@show RS_type.ϖ_λ₁λ₀
    #@show sum(RS_type.ϖ_λ₁λ₀)
    #@show RS_type.ϖ_Cabannes
    RS_type.i_λ₁λ₀ = index_raman_grid #index_raman_grid[end:-1:1]
    RS_type.n_Raman = length(RS_type.ϖ_λ₁λ₀)
    return nothing
end 

function getRamanSSProp!(RS_type::RRS_plus, depol) 
    @unpack n2, o2, 
            iBand, grid_in,
            greek_raman,
            ϖ_Cabannes, fscattRayl,
            ϖ_λ₁λ₀, i_λ₁λ₀,
            n_Raman =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    # determine Rayleigh scattering cross-section at central wavelength λ of the spectral band (assumed constant throughout the band)
    greek_raman = get_greek_raman(RS_type, n2, o2)
    t_w = Vector{Vector{FT}}(undef,0)
    t_i = Vector{Vector{Int}}(undef,0)
    nBands = length(iBand)
    ϖ_Cabannes = zeros(FT, nBands) 
    fscattRayl = zeros(FT, nBands)
    n_Raman = 0
    atmo_σ_Rayl = compute_optical_Rayl(λ, n2, o2)
    ϖ_Cabannes[iB] = compute_ϖ_Cabannes(RS_type, λ)
    for iB = 1:nBands
        _grid_in = grid_in[iB]
        λ = nm_per_m/(0.5*(_grid_in[1]+_grid_in[end]))
        #@show ϖ_Cabannes
        # determine RRS cross-sections to λ₀ from nSpecRaman wavelengths around λ₀  
        index_raman_grid, atmo_σ_RRS = compute_optical_RS!(RS_type, grid_in, λ, n2, o2)
        # declare ϖ_λ₁λ₀ to be a grid of length N_raman 
        t_ϖ_λ₁λ₀ = atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl; #the grid gets inverted because the central wavelength is now seen as the recipient of RRS from neighboring source wavelengths
        t_i_λ₁λ₀ = index_raman_grid; #index_raman_grid[end:-1:1];
        push!(t_w, t_ϖ_λ₁λ₀);
        push!(t_i, t_i_λ₁λ₀);
        t_n_Raman = length(t_ϖ_λ₁λ₀);
        n_Raman = (t_n_Raman>n_Raman) ? t_n_Raman : n_Raman
    end
    i_λ₁λ₀ = zeros(Int, nBands, n_Raman);
    ϖ_λ₁λ₀ = zeros( FT, nBands, n_Raman);
    for iB = 1:nBands
        for Δn = 1:length(t_i[iB])
            i_λ₁λ₀[iB,Δn] = t_i[iB][Δn];
            ϖ_λ₁λ₀[iB,Δn] = t_w[iB][Δn];
        end
    end
    @pack! RS_type = greek_raman, ϖ_Cabannes, fscattRayl, ϖ_λ₁λ₀, i_λ₁λ₀, n_Raman 
    return nothing;
end

function getRamanSSProp!(
            RS_type::VS_0to1_plus, depol, λ_inc)

    @unpack n2,o2,
            iBand, grid_in, bandSpecLim, 
            greek_raman, greek_raman_VS_n2, greek_raman_VS_o2,
            ϖ_Cabannes, fscattRayl,
            ϖ_λ₁λ₀, i_λ₁λ₀,
            ϖ_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_n2,
            ϖ_λ₁λ₀_VS_o2, i_λ₁λ₀_VS_o2,
            i_λ₁λ₀_all =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    iBand = []
    grid_in = []
    bandSpecLim = []

    nm_per_m = 1.e7;
    greek_raman = get_greek_raman(RS_type, n2, o2)
    greek_raman_VS_n2 = get_greek_raman_VS(RS_type, n2)
    greek_raman_VS_o2 = get_greek_raman_VS(RS_type, o2)
    #λ_scatt = 2.e7/(grid_in[1]+grid_in[end])     
    # determine Rayleigh scattering cross-section at incident wavelength λ 
    atmo_σ_Rayl = compute_optical_Rayl(λ_inc, n2, o2)
    # determine VRS_0to1 bands for λ₀ due to n2 and o2  
    nBand = 3; #for n2 and o2 each
    iBand = [1,2,3];
    ν̄ = nm_per_m/λ_inc
    push!(grid_in, ν̄:0.3:ν̄)
    molec = [n2,o2]
    #========================================================================#
    for mol in molec

        y1 = mol.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2;
        y2 = mol.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2;

        band_min = minimum([minimum(y1[y1.!=0]), minimum(y2[y2.!=0])]);
        band_min += nm_per_m/λ_inc;
        band_min -= 2.;
    
        band_max = maximum([maximum(y1[y1.!=0]), maximum(y2[y2.!=0])]);
        band_max += nm_per_m/λ_inc;
        band_max += 2.;

        push!(grid_in, band_min:0.3:band_max)

    end
    #========================================================================#
    FT = eltype(λ_inc);
    t_w_VS = Vector{Vector{FT}}(undef,0)
    t_i_VS = Vector{Vector{Int}}(undef,0)
    t_w_RVS = Vector{Vector{FT}}(undef,0)
    t_i_RVS = Vector{Vector{Int}}(undef,0)

    ϖ_Cabannes = zeros(FT, nBand) 
    fscattRayl = zeros(FT, nBand)
    #k_Rayl_scatt = zeros(FT, nBands)

    for iB = 1:nBand
        _grid_in = grid_in[iB]
        λ = nm_per_m/(0.5*(_grid_in[1]+_grid_in[end]))
        
        if iB==1
            #@show InelasticScattering.compute_ϖ_Cabannes(RS_type, depol, λ_inc, n2, o2)
            ϖ_Cabannes[iB] = 1.    #@show ϖ_Cabannes
        else
            ϖ_Cabannes[iB] = InelasticScattering.compute_ϖ_Cabannes(RS_type, λ_inc)
            # 1.
        end
        if iB==1
            t_ϖ_VRS = [0.]
            t_i_VRS = [1]
            t_ϖ_RVRS = [0.]
            t_i_RVRS = [1]
        else 
            index_VRSgrid_out, atmo_σ_VRS, index_RVRSgrid_out, atmo_σ_RVRS = 
                InelasticScattering.compute_optical_RS!(RS_type, 
                                            _grid_in, λ_inc, n2, o2)

            t_ϖ_VRS = atmo_σ_VRS/atmo_σ_Rayl
            t_i_VRS = index_VRSgrid_out
            t_ϖ_RVRS = atmo_σ_RVRS/atmo_σ_Rayl
            t_i_RVRS = index_RVRSgrid_out
        end
        push!(t_w_VS, t_ϖ_VRS);
        push!(t_i_VS, t_i_VRS);

        push!(t_w_RVS, t_ϖ_RVRS);
        push!(t_i_RVS, t_i_RVRS);

        #k_Rayl_scatt[iB] = compute_optical_Rayl(λ₀, n2, o2)/atmo_σ_Rayl
        
    end
    n_Raman = 1;
    bandSpecLim = [] # (1:τ_abs[iB])#zeros(Int64, iBand, 2) #Suniti: how to do this?
    #Suniti: make bandSpecLim a part of RS_type (including noRS) so that it can be passed into rt_kernel and elemental/doubling/interaction and postprocessing_vza without major syntax changes
    nSpec = 0;
    for iB in iBand
        nSpec0 = nSpec+1;
        nSpec += size(grid_in[iB], 1); # Number of spectral points
        push!(bandSpecLim,nSpec0:nSpec);                
    end

    i_λ₁λ₀ = zeros(Int, length(t_i_RVS[2])+length(t_i_RVS[3]));
    ϖ_λ₁λ₀ = zeros( FT, length(t_i_RVS[2])+length(t_i_RVS[3]));
    i_λ₁λ₀_VS_n2 = zeros(Int, length(t_i_VS[2])+length(t_i_VS[3]));
    ϖ_λ₁λ₀_VS_n2 = zeros( FT, length(t_i_VS[2])+length(t_i_VS[3]));
    i_λ₁λ₀_VS_o2 = zeros(Int, length(t_i_VS[2])+length(t_i_VS[3]));
    ϖ_λ₁λ₀_VS_o2 = zeros( FT, length(t_i_VS[2])+length(t_i_VS[3]));

    for Δn = 1:length(t_i_RVS[2])
        i_λ₁λ₀[Δn] = bandSpecLim[2][1] - 1 + t_i_RVS[2][Δn];
        ϖ_λ₁λ₀[Δn] = t_w_RVS[2][Δn];
    end
    for Δn = (length(t_i_RVS[2])+1):(length(t_i_RVS[2])+length(t_i_RVS[3]))
        i_λ₁λ₀[Δn] = bandSpecLim[3][1] - 1 + t_i_RVS[3][Δn-length(t_i_RVS[2])];
        ϖ_λ₁λ₀[Δn] = t_w_RVS[3][Δn-length(t_i_RVS[2])];
    end

    for Δn = 1:length(t_i_VS[2])
        i_λ₁λ₀_VS_n2[Δn] = bandSpecLim[2][1] - 1 + t_i_VS[2][Δn];
        ϖ_λ₁λ₀_VS_n2[Δn] = t_w_VS[2][Δn];
    end

    for Δn = (length(t_i_VS[2])+1):(length(t_i_VS[2])+length(t_i_VS[3]))
        i_λ₁λ₀_VS_o2[Δn] = bandSpecLim[3][1] - 1 + t_i_VS[3][Δn-length(t_i_VS[2])];
        ϖ_λ₁λ₀_VS_o2[Δn] = t_w_VS[3][Δn-length(t_i_VS[2])];
    end
    i_λ₁λ₀_all = unique(cat(i_λ₁λ₀, i_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_o2, dims = (1)))

    @pack! RS_type =
            iBand, grid_in, bandSpecLim,  
            i_λ₁λ₀, ϖ_λ₁λ₀, 
            i_λ₁λ₀_VS_n2, ϖ_λ₁λ₀_VS_n2,
            i_λ₁λ₀_VS_o2, ϖ_λ₁λ₀_VS_o2,
            i_λ₁λ₀_all,
            ϖ_Cabannes, fscattRayl,
            greek_raman,
            greek_raman_VS_n2,
            greek_raman_VS_o2
    return nothing
end

function getRamanSSProp!(
    RS_type::VS_1to0_plus, depol, λ_inc)

    @unpack n2,o2,
        iBand, grid_in, bandSpecLim, 
        greek_raman, greek_raman_VS_n2, greek_raman_VS_o2,
        ϖ_Cabannes, fscattRayl,
        ϖ_λ₁λ₀, i_λ₁λ₀,
        ϖ_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_n2,
        ϖ_λ₁λ₀_VS_o2, i_λ₁λ₀_VS_o2,
        i_λ₁λ₀_all =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    iBand = []
    grid_in = []
    bandSpecLim = []

    nm_per_m = 1.e7;
    greek_raman = get_greek_raman(RS_type, n2, o2)
    greek_raman_VS_n2 = get_greek_raman_VS(RS_type, n2)
    greek_raman_VS_o2 = get_greek_raman_VS(RS_type, o2)
    #λ_scatt = 2.e7/(grid_in[1]+grid_in[end])     
    # determine Rayleigh scattering cross-section at incident wavelength λ 
    atmo_σ_Rayl = compute_optical_Rayl(λ_inc, n2, o2)
    # determine VRS_0to1 bands for λ₀ due to n2 and o2  
    nBand = 3; #for n2 and o2 each
    iBand = [1,2,3];
    ν̄ = nm_per_m/λ_inc
    push!(grid_in, ν̄:0.3:ν̄)
    molec = [n2,o2]
    #========================================================================#
    for mol in molec

    y1 = mol.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2;
    y2 = mol.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2;

    band_min = minimum([minimum(y1[y1.!=0]), minimum(y2[y2.!=0])]);
    band_min += nm_per_m/λ_inc;
    band_min -= 2.;

    band_max = maximum([maximum(y1[y1.!=0]), maximum(y2[y2.!=0])]);
    band_max += nm_per_m/λ_inc;
    band_max += 2.;

    push!(grid_in, band_min:0.3:band_max)

    end
    #========================================================================#
    FT = eltype(λ_inc);
    t_w_VS = Vector{Vector{FT}}(undef,0)
    t_i_VS = Vector{Vector{Int}}(undef,0)
    t_w_RVS = Vector{Vector{FT}}(undef,0)
    t_i_RVS = Vector{Vector{Int}}(undef,0)

    ϖ_Cabannes = zeros(FT, nBand) 
    fscattRayl = zeros(FT, nBand)
    #k_Rayl_scatt = zeros(FT, nBands)

    for iB = 1:nBand
    _grid_in = grid_in[iB]
    λ = nm_per_m/(0.5*(_grid_in[1]+_grid_in[end]))

    if iB==1
        ϖ_Cabannes[iB] = InelasticScattering.compute_ϖ_Cabannes(RS_type, λ_inc)
        #@show ϖ_Cabannes
    else
        ϖ_Cabannes[iB] = 1.
    end
    if iB==1
        t_ϖ_VRS = [0.]
        t_i_VRS = [1]
        t_ϖ_RVRS = [0.]
        t_i_RVRS = [1]
    else 
        index_VRSgrid_out, atmo_σ_VRS, index_RVRSgrid_out, atmo_σ_RVRS = 
            InelasticScattering.compute_optical_RS!(RS_type, 
                                        _grid_in, λ_inc, n2, o2)

        t_ϖ_VRS = atmo_σ_VRS/atmo_σ_Rayl
        t_i_VRS = index_VRSgrid_out
        t_ϖ_RVRS = atmo_σ_RVRS/atmo_σ_Rayl
        t_i_RVRS = index_RVRSgrid_out
    end
    push!(t_w_VS, t_ϖ_VRS);
    push!(t_i_VS, t_i_VRS);

    push!(t_w_RVS, t_ϖ_RVRS);
    push!(t_i_RVS, t_i_RVRS);

    #k_Rayl_scatt[iB] = compute_optical_Rayl(λ₀, n2, o2)/atmo_σ_Rayl

    end
    n_Raman = 1;
    bandSpecLim = [] # (1:τ_abs[iB])#zeros(Int64, iBand, 2) #Suniti: how to do this?
    #Suniti: make bandSpecLim a part of RS_type (including noRS) so that it can be passed into rt_kernel and elemental/doubling/interaction and postprocessing_vza without major syntax changes
    nSpec = 0;
    for iB in iBand
    nSpec0 = nSpec+1;
    nSpec += size(grid_in[iB], 1); # Number of spectral points
    push!(bandSpecLim,nSpec0:nSpec);                
    end

    i_λ₁λ₀ = zeros(Int, length(t_i_RVS[2])+length(t_i_RVS[3]));
    ϖ_λ₁λ₀ = zeros( FT, length(t_i_RVS[2])+length(t_i_RVS[3]));
    i_λ₁λ₀_VS_n2 = zeros(Int, length(t_i_VS[2])+length(t_i_VS[3]));
    ϖ_λ₁λ₀_VS_n2 = zeros( FT, length(t_i_VS[2])+length(t_i_VS[3]));
    i_λ₁λ₀_VS_o2 = zeros(Int, length(t_i_VS[2])+length(t_i_VS[3]));
    ϖ_λ₁λ₀_VS_o2 = zeros( FT, length(t_i_VS[2])+length(t_i_VS[3]));

    for Δn = 1:length(t_i_RVS[2])
    i_λ₁λ₀[Δn] = bandSpecLim[2][1] - 1 + t_i_RVS[2][Δn];
    ϖ_λ₁λ₀[Δn] = t_w_RVS[2][Δn];
    end
    for Δn = (length(t_i_RVS[2])+1):(length(t_i_RVS[2])+length(t_i_RVS[3]))
    i_λ₁λ₀[Δn] = bandSpecLim[3][1] - 1 + t_i_RVS[3][Δn-length(t_i_RVS[2])];
    ϖ_λ₁λ₀[Δn] = t_w_RVS[3][Δn-length(t_i_RVS[2])];
    end

    for Δn = 1:length(t_i_VS[2])
    i_λ₁λ₀_VS_n2[Δn] = bandSpecLim[2][1] - 1 + t_i_VS[2][Δn];
    ϖ_λ₁λ₀_VS_n2[Δn] = t_w_VS[2][Δn];
    end

    for Δn = (length(t_i_VS[2])+1):(length(t_i_VS[2])+length(t_i_VS[3]))
    i_λ₁λ₀_VS_o2[Δn] = bandSpecLim[3][1] - 1 + t_i_VS[3][Δn-length(t_i_VS[2])];
    ϖ_λ₁λ₀_VS_o2[Δn] = t_w_VS[3][Δn-length(t_i_VS[2])];
    end
    i_λ₁λ₀_all = unique(cat(i_λ₁λ₀, i_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_o2, dims = (1)))

    @pack! RS_type =
        iBand, grid_in, bandSpecLim,  
        i_λ₁λ₀, ϖ_λ₁λ₀, 
        i_λ₁λ₀_VS_n2, ϖ_λ₁λ₀_VS_n2,
        i_λ₁λ₀_VS_o2, ϖ_λ₁λ₀_VS_o2,
        i_λ₁λ₀_all,
        ϖ_Cabannes, fscattRayl,
        greek_raman,
        greek_raman_VS_n2,
        greek_raman_VS_o2
    return nothing
end

function getRamanSSProp!(
    RS_type::Union{VS_0to1_plus, VS_1to0_plus}, depol, λ_inc,  target_grid)

    @unpack n2,o2,
        iBand, grid_in, bandSpecLim, 
        greek_raman, greek_raman_VS_n2, greek_raman_VS_o2,
        ϖ_Cabannes, fscattRayl,
        ϖ_λ₁λ₀, i_λ₁λ₀,
        ϖ_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_n2,
        ϖ_λ₁λ₀_VS_o2, i_λ₁λ₀_VS_o2,
        i_λ₁λ₀_all =  RS_type
    #n2, o2 = getRamanAtmoConstants(1.e7/λ, T)
    iBand = []
    grid_in = []
    bandSpecLim = []

    nm_per_m = 1.e7;
    greek_raman = get_greek_raman(RS_type, n2, o2)
    greek_raman_VS_n2 = get_greek_raman_VS(RS_type, n2)
    greek_raman_VS_o2 = get_greek_raman_VS(RS_type, o2)
    #λ_scatt = 2.e7/(grid_in[1]+grid_in[end])     
    # determine Rayleigh scattering cross-section at incident wavelength λ 
    atmo_σ_Rayl = compute_optical_Rayl(λ_inc, n2, o2)
    # determine VRS_0to1 bands for λ₀ due to n2 and o2  
    nBand = 2; #for n2 and o2 each
    iBand = [1,2];
    ν̄ = nm_per_m/λ_inc
    push!(grid_in, ν̄:0.3:ν̄)
    push!(grid_in, target_grid)
    molec = [n2,o2]
    #========================================================================#
    #=
    for mol in molec

        y1 = mol.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2;
        y2 = mol.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2;

        band_min = minimum([minimum(y1[y1.!=0]), minimum(y2[y2.!=0])]);
        band_min += nm_per_m/λ_inc;
        band_min -= 2.;

        band_max = maximum([maximum(y1[y1.!=0]), maximum(y2[y2.!=0])]);
        band_max += nm_per_m/λ_inc;
        band_max += 2.;

        push!(grid_in, band_min:0.3:band_max)

    end
    =#
    #========================================================================#
    FT = eltype(λ_inc);
    t_w_VS_n2 = Vector{Vector{FT}}(undef,0)
    t_i_VS_n2 = Vector{Vector{Int}}(undef,0)
    t_w_VS_o2 = Vector{Vector{FT}}(undef,0)
    t_i_VS_o2 = Vector{Vector{Int}}(undef,0)
    t_w_RVS = Vector{Vector{FT}}(undef,0)
    t_i_RVS = Vector{Vector{Int}}(undef,0)

    ϖ_Cabannes = zeros(FT, nBand) 
    fscattRayl = zeros(FT, nBand)
    #k_Rayl_scatt = zeros(FT, nBands)

    for iB = 1:nBand
        _grid_in = grid_in[iB]
        #λ = nm_per_m/(0.5*(_grid_in[1]+_grid_in[end]))

        if iB==1
            #@show InelasticScattering.compute_ϖ_Cabannes(RS_type, λ_inc)
            ϖ_Cabannes[iB] = InelasticScattering.compute_ϖ_Cabannes(RS_type, λ_inc)
            #@show ϖ_Cabannes
        else
            ϖ_Cabannes[iB] = 1.
        end
        if iB==1
            t_ϖ_VRS_n2 = [0.]
            t_i_VRS_n2 = [1]
            t_ϖ_VRS_o2 = [0.]
            t_i_VRS_o2 = [1]
            t_ϖ_RVRS = [0.]
            t_i_RVRS = [1]
        else 
            index_RVRSgrid_out, atmo_σ_RVRS = 
                InelasticScattering.compute_optical_RVRS!(RS_type, 
                                            _grid_in, λ_inc, n2, o2)
            index_VRSgrid_n2_out, atmo_σ_VRS_n2 = 
                InelasticScattering.compute_optical_VRS!(RS_type, 
                                        _grid_in, λ_inc, n2) 
            index_VRSgrid_o2_out, atmo_σ_VRS_o2 = 
                InelasticScattering.compute_optical_VRS!(RS_type, 
                                        _grid_in, λ_inc, o2)                                                            

            t_ϖ_VRS_n2 = atmo_σ_VRS_n2/atmo_σ_Rayl
            t_i_VRS_n2 = index_VRSgrid_n2_out
            t_ϖ_VRS_o2 = atmo_σ_VRS_o2/atmo_σ_Rayl
            t_i_VRS_o2 = index_VRSgrid_o2_out
            t_ϖ_RVRS = atmo_σ_RVRS/atmo_σ_Rayl
            t_i_RVRS = index_RVRSgrid_out
        end
        push!(t_w_VS_n2, t_ϖ_VRS_n2);
        push!(t_i_VS_n2, t_i_VRS_n2);
        push!(t_w_VS_o2, t_ϖ_VRS_o2);
        push!(t_i_VS_o2, t_i_VRS_o2);
        push!(t_w_RVS, t_ϖ_RVRS);
        push!(t_i_RVS, t_i_RVRS);

        #k_Rayl_scatt[iB] = compute_optical_Rayl(λ₀, n2, o2)/atmo_σ_Rayl

    end
    n_Raman = 1;
    bandSpecLim = [] # (1:τ_abs[iB])#zeros(Int64, iBand, 2) #Suniti: how to do this?
    #Suniti: make bandSpecLim a part of RS_type (including noRS) so that it can be passed into rt_kernel and elemental/doubling/interaction and postprocessing_vza without major syntax changes
    nSpec = 0;
    for iB in iBand
        nSpec0 = nSpec+1;
        nSpec += size(grid_in[iB], 1); # Number of spectral points
        push!(bandSpecLim,nSpec0:nSpec);                
    end

    i_λ₁λ₀ = zeros(Int, length(t_i_RVS[2]));
    ϖ_λ₁λ₀ = zeros( FT, length(t_i_RVS[2]));
    i_λ₁λ₀_VS_n2 = zeros(Int, length(t_i_VS_n2[2]));
    ϖ_λ₁λ₀_VS_n2 = zeros( FT, length(t_i_VS_n2[2]));
    i_λ₁λ₀_VS_o2 = zeros(Int, length(t_i_VS_o2[2]));
    ϖ_λ₁λ₀_VS_o2 = zeros( FT, length(t_i_VS_o2[2]));

    for Δn = 1:length(t_i_RVS[2])
        i_λ₁λ₀[Δn] = bandSpecLim[2][1] - 1 + t_i_RVS[2][Δn];
        ϖ_λ₁λ₀[Δn] = t_w_RVS[2][Δn];
    end

    for Δn = 1:length(t_i_VS_n2[2])
        i_λ₁λ₀_VS_n2[Δn] = bandSpecLim[2][1] - 1 + t_i_VS_n2[2][Δn];
        ϖ_λ₁λ₀_VS_n2[Δn] = t_w_VS_n2[2][Δn];
    end

    for Δn = 1:length(t_i_VS_o2[2])
        i_λ₁λ₀_VS_o2[Δn] = bandSpecLim[2][1] - 1 + t_i_VS_o2[2][Δn];
        ϖ_λ₁λ₀_VS_o2[Δn] = t_w_VS_o2[2][Δn];
    end
    i_λ₁λ₀_all = unique(cat(i_λ₁λ₀, i_λ₁λ₀_VS_n2, i_λ₁λ₀_VS_o2, dims = (1)))

    @pack! RS_type =
        iBand, grid_in, bandSpecLim,  
        i_λ₁λ₀, ϖ_λ₁λ₀, 
        i_λ₁λ₀_VS_n2, ϖ_λ₁λ₀_VS_n2,
        i_λ₁λ₀_VS_o2, ϖ_λ₁λ₀_VS_o2,
        i_λ₁λ₀_all,
        ϖ_Cabannes, fscattRayl,
        greek_raman,
        greek_raman_VS_n2,
        greek_raman_VS_o2
    return nothing
end

