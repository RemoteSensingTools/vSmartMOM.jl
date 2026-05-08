#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#
### New update:
function rt_kernel_multisensor!(RS_type::noRS{FT},
                    sensor_levels,
                    pol_type, SFI,
                    added_layer,
                    composite_layer,
                    computed_layer_properties::CoreScatteringOpticalProperties,
                    scattering_interface,
                    τ_sum,
                    m, quad_points,
                    I_static,
                    architecture,
                    qp_μN, iz, arr_type;
                    dτ_max_threshold::Union{Nothing,Real} = nothing,
                    dτ_min_floor::Union{Nothing,Real} = nothing) where {FT}

    (; μ₀) = quad_points
    (; F₀) = RS_type
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
    # Centralised dτ/ndoubl helper (see noRS rt_kernel! variant).
    dτ, ndoubl = get_dtau_ndoubl(computed_layer_properties, quad_points;
                                 dτ_max_threshold = dτ_max_threshold,
                                 dτ_min_floor = dτ_min_floor)
    scatter = true # edit later
    arr_type = array_type(architecture)
    expk = arr_type(exp.(-dτ /μ₀))
    

    # If there is scattering, perform the elemental and doubling steps
    if scatter
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, dτ, F₀,
                                        computed_layer_properties, 
                                        m, ndoubl, scatter, quad_points,  
                                        added_layer,  architecture)
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, 
                                        expk, ndoubl, 
                                        added_layer, 
                                        I_static, architecture)
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        temp = collect(exp.(-τ_λ./qp_μN'))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
        end
    end
    #M1 = collect(added_layer.t⁺⁺)
    #M2 = collect(added_layer.r⁺⁻)
    #M3 = collect(added_layer.j₀⁻)
    #M4 = collect(added_layer.j₀⁺)
    #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1]
    
    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the bottom composite layer
    if (iz == 1)
        for ims=1:length(sensor_levels)
            if(sensor_levels[ims]==0)
                #@show sensor_levels[ims], iz, (iz==1)
                composite_layer.botT⁺⁺[ims][:], composite_layer.botT⁻⁻[ims][:] = 
                    collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                composite_layer.botR⁻⁺[ims][:], composite_layer.botR⁺⁻[ims][:] = 
                    collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                composite_layer.botJ₀⁺[ims][:], composite_layer.botJ₀⁻[ims][:] = 
                    collect(added_layer.j₀⁺), collect(added_layer.j₀⁻)
            else
                #@show sensor_levels[ims], iz, (iz==1)
                composite_layer.topT⁺⁺[ims][:], composite_layer.topT⁻⁻[ims][:] = 
                    collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                composite_layer.topR⁻⁺[ims][:], composite_layer.topR⁺⁻[ims][:] = 
                    collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                composite_layer.topJ₀⁺[ims][:], composite_layer.topJ₀⁻[ims][:] = 
                    collect(added_layer.j₀⁺), collect(added_layer.j₀⁻)
            end
        end
    # If this is not the TOA, perform the interaction step
    else
        #@timeit "interaction_multisensor" interaction_multisensor!(RS_type, sensor_levels, scattering_interface, SFI, composite_layer, added_layer, I_static)
        for ims=1:length(sensor_levels)
            if sensor_levels[ims]==0
                #@show sensor_levels[ims], iz, (iz!=1)
                @timeit "interaction_multisensor" interaction_bot!(ims, 
                RS_type, 
                scattering_interface, 
                SFI, 
                composite_layer, 
                added_layer, 
                I_static,
                arr_type)
            else
                if sensor_levels[ims]==(iz-1) #include ims==Nz with ims==0
                    #@show sensor_levels[ims], iz, (iz==sensor_levels[ims]+1)
                    composite_layer.botT⁺⁺[ims][:], composite_layer.botT⁻⁻[ims][:] = 
                        collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                    composite_layer.botR⁻⁺[ims][:], composite_layer.botR⁺⁻[ims][:] = 
                        collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                    composite_layer.botJ₀⁺[ims][:], composite_layer.botJ₀⁻[ims][:] = 
                        collect(added_layer.j₀⁺), collect(added_layer.j₀⁻ )
                elseif sensor_levels[ims]<(iz-1) 
                    #@show sensor_levels[ims], iz, (iz>sensor_levels[ims]+1)
                    @timeit "interaction_multisensor" interaction_bot!(ims, 
                                                                    RS_type, 
                                                                    scattering_interface, 
                                                                    SFI, 
                                                                    composite_layer, 
                                                                    added_layer, 
                                                                    I_static,
                                                                    arr_type)
                elseif sensor_levels[ims]>=iz 
                    #@show sensor_levels[ims], iz, (iz<=sensor_levels[ims])
                    @timeit "interaction_multisensor" interaction_top!(ims, 
                                                                    RS_type, 
                                                                    scattering_interface, 
                                                                    SFI, 
                                                                    composite_layer, 
                                                                    added_layer, 
                                                                    I_static,
                                                                    arr_type)    
                end
                if iz==2
                    M1 = (composite_layer.botT⁺⁺[1]);
                    M2 = (composite_layer.botR⁺⁻[1]);
                    M3 = (composite_layer.botT⁻⁻[1]);
                    M4 = (composite_layer.botR⁻⁺[1]);
                    M5 = (composite_layer.botJ₀⁻[1]);
                    M6 = (composite_layer.botJ₀⁺[1]);
                    #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1], M5[1,1,1], M6[1,1,1]
                end
            end
        end
    end
end

function rt_kernel_multisensor!(RS_type::Union{RRS{FT}, VS_0to1_plus{FT}, VS_1to0_plus{FT}},
                                sensor_levels,
                                pol_type,
                                SFI,
                                added_layer,
                                composite_layer,
                                computed_layer_properties::CoreScatteringOpticalProperties,
                                scattering_interface,
                                τ_sum,
                                m,
                                quad_points,
                                I_static,
                                architecture,
                                qp_μN,
                                iz,
                                arr_type;
                                dτ_max_threshold::Union{Nothing,Real} = nothing,
                                dτ_min_floor::Union{Nothing,Real} = nothing)  where {FT}
    (; μ₀) = quad_points
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
    # Centralised dτ/ndoubl helper (see noRS rt_kernel! variant).
    dτ, ndoubl = get_dtau_ndoubl(computed_layer_properties, quad_points;
                                 dτ_max_threshold = dτ_max_threshold,
                                 dτ_min_floor = dτ_min_floor)
    scatter = true # edit later
    expk = arr_type(exp.(-dτ /μ₀))

    (; Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, F₀) = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, dτ, ϖ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                F₀,
                                                m, ndoubl, scatter, 
                                                quad_points,  added_layer,  
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ, F₀, computed_layer_properties, m, ndoubl, scatter, quad_points,  added_layer,  architecture)
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        added_layer.ier⁻⁺[:] .= 0;
        added_layer.ier⁺⁻[:] .= 0;
        added_layer.ieJ₀⁻[:] .= 0;
        added_layer.iet⁻⁻[:] .= 0;
        added_layer.iet⁺⁺[:] .= 0;
        added_layer.ieJ₀⁺[:] .= 0;
        temp = collect(exp.(-τ_λ./qp_μN'))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
        end
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        for ims=1:length(sensor_levels)
            if(sensor_levels[ims]==0)
                # bottom composite layer for TOA/BOA sensors
                composite_layer.botT⁺⁺[ims][:], composite_layer.botT⁻⁻[ims][:] = 
                    collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                composite_layer.botR⁻⁺[ims][:], composite_layer.botR⁺⁻[ims][:] = 
                    collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                composite_layer.botJ₀⁺[ims][:], composite_layer.botJ₀⁻[ims][:] = 
                    collect(added_layer.j₀⁺), collect(added_layer.j₀⁻ )
                composite_layer.botieT⁺⁺[ims][:], composite_layer.botieT⁻⁻[ims][:] = 
                    collect(added_layer.iet⁺⁺), collect(added_layer.iet⁻⁻)
                composite_layer.botieR⁻⁺[ims][:], composite_layer.botieR⁺⁻[ims][:] = 
                    collect(added_layer.ier⁻⁺), collect(added_layer.ier⁺⁻)
                composite_layer.botieJ₀⁺[ims][:], composite_layer.botieJ₀⁻[ims][:] = 
                    collect(added_layer.ieJ₀⁺), collect(added_layer.ieJ₀⁻ )
            else
                composite_layer.topT⁺⁺[ims][:], composite_layer.topT⁻⁻[ims][:] = 
                    collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                composite_layer.topR⁻⁺[ims][:], composite_layer.topR⁺⁻[ims][:] = 
                    collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                composite_layer.topJ₀⁺[ims][:], composite_layer.topJ₀⁻[ims][:] = 
                    collect(added_layer.j₀⁺), collect(added_layer.j₀⁻ )
                composite_layer.topieT⁺⁺[ims][:], composite_layer.topieT⁻⁻[ims][:] = 
                    collect(added_layer.iet⁺⁺), collect(added_layer.iet⁻⁻)
                composite_layer.topieR⁻⁺[ims][:], composite_layer.topieR⁺⁻[ims][:] = 
                    collect(added_layer.ier⁻⁺), collect(added_layer.ier⁺⁻)
                composite_layer.topieJ₀⁺[ims][:], composite_layer.topieJ₀⁻[ims][:] = 
                    collect(added_layer.ieJ₀⁺), collect(added_layer.ieJ₀⁻ )
            end
        end
    # If this is not the TOA, perform the interaction step
    else
        #@timeit "interaction_multisensor" interaction_multisensor!(RS_type, sensor_levels, scattering_interface, SFI, composite_layer, added_layer, I_static)
        for ims=1:length(sensor_levels)
            if sensor_levels[ims]==(iz-1) #include ims==Nz with ims==0
                composite_layer.botT⁺⁺[ims][:], composite_layer.botT⁻⁻[ims][:] = 
                    collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                composite_layer.botR⁻⁺[ims][:], composite_layer.botR⁺⁻[ims][:] = 
                    collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                composite_layer.botJ₀⁺[ims][:], composite_layer.botJ₀⁻[ims][:] = 
                    collect(added_layer.j₀⁺), collect(added_layer.j₀⁻ )
                composite_layer.botieT⁺⁺[ims][:], composite_layer.botieT⁻⁻[ims][:] = 
                    collect(added_layer.iet⁺⁺), collect(added_layer.iet⁻⁻)
                composite_layer.botieR⁻⁺[ims][:], composite_layer.botieR⁺⁻[ims][:] = 
                    collect(added_layer.ier⁻⁺), collect(added_layer.ier⁺⁻)
                composite_layer.botieJ₀⁺[ims][:], composite_layer.botieJ₀⁻[ims][:] = 
                    collect(added_layer.ieJ₀⁺), collect(added_layer.ieJ₀⁻ )
            elseif sensor_levels[ims]<(iz-1) 
                @timeit "interaction_multisensor" interaction_bot!(ims, 
                                                                RS_type, 
                                                                scattering_interface, 
                                                                SFI, 
                                                                composite_layer, 
                                                                added_layer, 
                                                                I_static,
                                                                arr_type)
            elseif sensor_levels[ims]>=iz 
                @timeit "interaction_multisensor" interaction_top!(ims, 
                                                                RS_type, 
                                                                scattering_interface, 
                                                                SFI, 
                                                                composite_layer, 
                                                                added_layer, 
                                                                I_static,
                                                                arr_type)                
            end
        end
    end
end
#=
function rt_kernel_multisensor!(
            RS_type::Union{VS_0to1_plus{FT}, VS_1to0_plus{FT}},
            sensor_levels,              
            pol_type, 
            SFI, 
            added_layer, 
            composite_layer, 
            computed_layer_properties::CoreScatteringOpticalProperties, 
            scattering_interface, 
            τ_sum,
            m, 
            quad_points, 
            I_static, 
            architecture, 
            qp_μN, 
            iz, 
            arr_type)  where {FT}
    @unpack qp_μ, μ₀ = quad_points
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    # SUNITI, check? Also, better to write function here
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    scatter = true # edit later
    #arr_type = array_type(architecture)
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    expk = arr_type(exp.(-dτ /μ₀))

    @unpack Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, dτ, ϖ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                m, ndoubl, scatter, 
                                                quad_points,  added_layer,  
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ, computed_layer_properties, m, ndoubl, scatter, quad_points,  added_layer,  architecture)
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        added_layer.ier⁻⁺[:] .= 0;
        added_layer.ier⁺⁻[:] .= 0;
        added_layer.ieJ₀⁻[:] .= 0;
        added_layer.iet⁻⁻[:] .= 0;
        added_layer.iet⁺⁺[:] .= 0;
        added_layer.ieJ₀⁺[:] .= 0;
        temp = collect(exp.(-τ./qp_μN'))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
        end
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        for ims=1:length(sensor_levels)
            if(ims==0) # bottom composite layer for satellite/groundbased sensors
                composite_layer.botT⁺⁺[ims][:], composite_layer.botT⁻⁻[ims][:] = 
                    collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                composite_layer.botR⁻⁺[ims][:], composite_layer.botR⁺⁻[ims][:] = 
                    collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                composite_layer.botJ₀⁺[ims][:], composite_layer.botJ₀⁻[ims][:] = 
                    collect(added_layer.j₀⁺), collect(added_layer.j₀⁻ )
                composite_layer.botieT⁺⁺[ims][:], composite_layer.botieT⁻⁻[ims][:] = 
                    collect(added_layer.iet⁺⁺), collect(added_layer.iet⁻⁻)
                composite_layer.botieR⁻⁺[ims][:], composite_layer.botieR⁺⁻[ims][:] = 
                    collect(added_layer.ier⁻⁺), collect(added_layer.ier⁺⁻)
                composite_layer.botieJ₀⁺[ims][:], composite_layer.botieJ₀⁻[ims][:] = 
                    collect(added_layer.ieJ₀⁺), collect(added_layer.ieJ₀⁻ )
            else
                composite_layer.topT⁺⁺[ims][:], composite_layer.topT⁻⁻[ims][:] = 
                    collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                composite_layer.topR⁻⁺[ims][:], composite_layer.topR⁺⁻[ims][:] = 
                    collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                composite_layer.topJ₀⁺[ims][:], composite_layer.topJ₀⁻[ims][:] = 
                    collect(added_layer.j₀⁺), collect(added_layer.j₀⁻ )
                composite_layer.topieT⁺⁺[ims][:], composite_layer.topieT⁻⁻[ims][:] = 
                    collect(added_layer.iet⁺⁺), collect(added_layer.iet⁻⁻)
                composite_layer.topieR⁻⁺[ims][:], composite_layer.topieR⁺⁻[ims][:] = 
                    collect(added_layer.ier⁻⁺), collect(added_layer.ier⁺⁻)
                composite_layer.topieJ₀⁺[ims][:], composite_layer.topieJ₀⁻[ims][:] = 
                    collect(added_layer.ieJ₀⁺), collect(added_layer.ieJ₀⁻ )
            end
        end
        
    else # If this is not the TOA, perform the interaction step
        for ims=1:length(sensor_levels)
            if sensor_levels[ims]==(iz+1) #include ims==Nz with ims==0
                composite_layer.botT⁺⁺[ims][:], composite_layer.botT⁻⁻[ims][:] = 
                    collect(added_layer.t⁺⁺), collect(added_layer.t⁻⁻)
                composite_layer.botR⁻⁺[ims][:], composite_layer.botR⁺⁻[ims][:] = 
                    collect(added_layer.r⁻⁺), collect(added_layer.r⁺⁻)
                composite_layer.botJ₀⁺[ims][:], composite_layer.botJ₀⁻[ims][:] = 
                    collect(added_layer.j₀⁺), collect(added_layer.j₀⁻ )
                composite_layer.botieT⁺⁺[ims][:], composite_layer.botieT⁻⁻[ims][:] = 
                    collect(added_layer.iet⁺⁺), collect(added_layer.iet⁻⁻)
                composite_layer.botieR⁻⁺[ims][:], composite_layer.botieR⁺⁻[ims][:] = 
                    collect(added_layer.ier⁻⁺), collect(added_layer.ier⁺⁻)
                composite_layer.botieJ₀⁺[ims][:], composite_layer.botieJ₀⁻[ims][:] = 
                    collect(added_layer.ieJ₀⁺), collect(added_layer.ieJ₀⁻ )
            elseif sensor_levels[ims]>(iz+1) 
                @timeit "interaction_multisensor" interaction_bot!(ims, 
                                                                RS_type, 
                                                                scattering_interface, 
                                                                SFI, 
                                                                composite_layer, 
                                                                added_layer, 
                                                                I_static,
                                                                arr_type)
            elseif sensor_levels[ims]<=iz 
                @timeit "interaction_multisensor" interaction_top!(ims, 
                                                                RS_type, 
                                                                scattering_interface, 
                                                                SFI, 
                                                                composite_layer, 
                                                                added_layer, 
                                                                I_static,
                                                                arr_type)                
            end
        end
    end
end
=#
