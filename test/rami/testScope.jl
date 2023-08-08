using CanopyLayers

FT = Float64;
nLayer = 20
LAI = 3.0

can = create_canopy_rt(FT, nLayer=nLayer, LAI=LAI);
wls = create_wave_length(FT, collect(456:20:533));

# 2. create rt_dims from wls and can
rt_dim = create_rt_dims(can, wls);

# 3. create can_rad, can_opt, and etc from rt_dim and wls
can_rad = create_canopy_rads(FT, rt_dim);
can_opt = create_canopy_opticals(FT, rt_dim);
in_rad  = create_incoming_radiation(wls);
in_rad.E_diffuse .= 0.06 * in_rad.E_direct
soil    = SoilOpticals{FT}(wls);
angles  = SolarAngles{FT}();
angles.sza = 30.0
rt_con  = create_rt_cache(FT, rt_dim);

# Create an array of standard leaves
leaves = [create_leaf_bios(FT, rt_dim) for i in 1:nLayer];
for i in 1:nLayer
    fluspect!(leaves[i], wls);
end
soil.ρ_SW .=0.04439
for leaf in leaves
    leaf.τ_SW .= 0.08408;
    leaf.ρ_SW .= 0.1323;
end


vza_start = 1.0  # rami_measures[1]["vza_start"]["value"]
vza_end   = 75.0 # rami_measures[1]["vza_end"]["value"]
vza_step  = 2.0  # rami_measures[1]["vza_step"]["value"]
vzas = collect(vza_start:vza_step:vza_end)
vza_  = [reverse(vzas); -vzas; reverse(vzas); -vzas]
vaz_  = [repeat([0.0], length(vzas)); repeat([0.0], length(vzas));repeat([90.0], length(vzas));repeat([-90.0], length(vzas)) ]


R = similar(vza_)

for j in eachindex(vza_)
    
        angles.vza = vza_[j]
        angles.vaa = vaz_[j]
        angles.raa = vaz_[j]


        canopy_geometry!(can, angles, can_opt, rt_con);
        # 2. Update scattering coefficients (required)
        canopy_matrices!(leaves, can_opt);
        # 3. Simulate short wave simulation (required)
        short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
        R[j] = can_rad.alb_obs[1]

end

R2 = similar(R)
can.LIDFa = 1.0
can.LIDFb = 0.0
can = Canopy4RT{FT}(LAI    = LAI   ,nLayer = nLayer, LIDFa = 1.0, LIDFb = 0.0)
# create_canopy_rt(FT, nLayer=nLayer, LAI=LAI, LIDFa = 1.0, LIDFb = 0.0);
for j in eachindex(vza_)
    
    angles.vza = vza_[j]
    angles.vaa = vaz_[j]
    angles.raa = vaz_[j]


    canopy_geometry!(can, angles, can_opt, rt_con);
    # 2. Update scattering coefficients (required)
    canopy_matrices!(leaves, can_opt);
    # 3. Simulate short wave simulation (required)
    short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
    R2[j] = can_rad.alb_obs[1]

end

R3 = similar(R)
can = Canopy4RT{FT}(LAI    = LAI   ,nLayer = nLayer, LIDFa = -1.0, LIDFb = 0.0)
# create_canopy_rt(FT, nLayer=nLayer, LAI=LAI, LIDFa = 1.0, LIDFb = 0.0);
for j in eachindex(vza_)
    
    angles.vza = vza_[j]
    angles.vaa = vaz_[j]
    angles.raa = vaz_[j]

    canopy_geometry!(can, angles, can_opt, rt_con);
    # 2. Update scattering coefficients (required)
    canopy_matrices!(leaves, can_opt);
    # 3. Simulate short wave simulation (required)
    short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
    R3[j] = can_rad.alb_obs[1]

end

plot(hdrf35,label="VS Erectophile")
plot!(hdrf25,label="VS Planophile")
plot!(hdrf45,label="VS Uniform")

plot!(R3,label="SCOPE Erectophile")
plot!(R2,label="SCOPE Planophile")
plot!(R,label="SCOPE Uniform")