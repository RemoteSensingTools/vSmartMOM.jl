using Plots 
using LegendrePolynomials
using DelimitedFiles
#using PlotlyJS

fit_window = [758.0, 759.2]
#sza_str   = "50"
#psurf_str = "1000"
#alb_str   = "0p2"
FT = Float32
a=[]
for i=0:20
    push!(a, acosd(i/20))  
end
sza=reverse(Int.(ceil.(a[8:21])))
ρ = zeros(FT,15)
ρ_str = []
for iρ = 1:15
    ρ[iρ] = (iρ-1)*0.05
    push!(ρ_str, replace(string(round(ρ[iρ], digits=2)),"."=>"p"))
end
psurf=[1000 750 500]
xSIF = zeros(length(ρ_str), length(sza))
xSIF_noRS = zeros(length(ρ_str), length(sza))
isurf = 3
iρ = 1
iA = 1
for isurf = 1:3 # 1:1 #2:2 # 
    for iρ = 1:15 #3 #1:15
        for iA = 1:14
            sza_str  = string(sza[iA])
            alb_str  = ρ_str[iρ]
            psurf_str= string(psurf[isurf])

            fname0   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1   = "/home/sanghavi/RamanSIFgrid/raylSIF_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"

            fname0_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_nors_ABO2.dat"
            fname1_0   = "/home/sanghavi/RamanSIFgrid/rayl_sza"*sza_str*"_alb"*alb_str*"_psurf"*psurf_str*"hpa_rrs_ABO2.dat"


            specNoRS = readdlm(fname0)
            specRRS  = readdlm(fname1)

            specNoRS0 = readdlm(fname0_0)
            specRRS0  = readdlm(fname1_0)

            wl = 1e7./specNoRS[:,1]
            radRRS_I = specRRS[:,2] + specRRS[:,5] 
            radRRS_I_el = specRRS[:,2] 
            radRRS_I_ie = specRRS[:,5] 
            radnoRS_I = specNoRS[:,2] 

            radRRS_I0 = specRRS0[:,2] + specRRS0[:,5] 
            radRRS_I_el0 = specRRS0[:,2] 
            radRRS_I_ie0 = specRRS0[:,5] 
            radnoRS_I0 = specNoRS0[:,2] 

            F0 = specNoRS[:,5] 
            # Find indices for wavelength window:
            ind = findall(x->x>fit_window[1] && x<fit_window[2], wl);

            # Fit SIF
            # Grid for Legendre Polynomials:
            iLeg = range(-1,1,length(ind))
            # Take noRS run as solar reference (can in principle also use raw solar spectrum, shouldn't matter)
            Io = specNoRS[ind,5] #radnoRS_I[ind]
            # Define polynomial terms
            poly = Pl.(iLeg,(0:3)');
            # Multiply polynomial with solar spectrum
            K_ = Io .* poly 
            # Add a constant SIF term to the Jacobian (can add Shape later)
            K = [K_ ones(length(ind))];
            K1 = [K_ ones(length(ind))];
            # Fit with Least squares:
            x = K \ radRRS_I[ind];

            # COmpute fitted spectrum:
            radRRS_I_fit = K * x;
            # Fitted SIF:
            SIF_fit = x[end]
            relative_SIF = SIF_fit / maximum(radRRS_I[ind])

            # Fit SIF without Raman
            # Fit with Least squares:
            x_noRS = K \ radnoRS_I[ind];

            # COmpute fitted spectrum:
            radnoRS_I_fit0 = K * x_noRS;
            # Fitted SIF:
            SIF_fit_noRS = x_noRS[end]
            relative_SIF_noRS = SIF_fit_noRS / maximum(radnoRS_I[ind])

            xSIF[iρ, iA] = relative_SIF # SIF_fit
            xSIF_noRS[iρ, iA] = relative_SIF_noRS #SIF_fit_noRS
#        end
#    end
#end

#heatmap(x, x, f, c = :thermal)
heatmap(ρ, sza, xSIF', c = :thermal)

#N = 100000
plot(scattergl(
    x=ρ[1:13], y=sza, mode="markers",
    marker=attr(color=xSIF, colorscale="Viridis", line_width=1)
))

l = @layout [a ; b]
p1 = plot(wl[ind], radRRS_I[ind], label="Measured")
plot!(wl[ind], radRRS_I_fit, label="Fitted")
p2 = plot(wl[ind], (radRRS_I[ind] .- radRRS_I_fit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)



# Fit without SIF:
#x2 = K_ \ radRRS_I[ind];
#radStupidFit = K_ * x2;
#l = @layout [a ; b]
# Fit with Least squares:



p1 = plot(wl[ind], radRRS_I0[ind], label="Measured")
plot!(wl[ind], radRRS_I_fit0, label="Fit without SIF")
p2 = plot(wl[ind], (radRRS_I0[ind] .- radRRS_I_fit0)./radRRS_I0[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)

#plot!(wl[ind], radStupidFit, label="Fit without SIF")
#p2 = plot(wl[ind], (radRRS_I[ind] .- radStupidFit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
#plot(p1, p2,  layout = l)
 

# Fit without Raman:
#x2 = K_ \ radRRS_I[ind];
#radStupidFit = K_ * x2;
#l = @layout [a ; b]
# Fit with Least squares:
#x_noRS = K \ radnoRS_I[ind];

# COmpute fitted spectrum:
#radnoRS_I_fit = K * x_noRS;


p1 = plot(wl[ind], radnoRS_I[ind], label="Measured")
plot!(wl[ind], radnoRS_I_fit, label="Fit without SIF")
p2 = plot(wl[ind], (radnoRS_I[ind] .- radnoRS_I_fit)./radnoRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)

#plot!(wl[ind], radStupidFit, label="Fit without SIF")
#p2 = plot(wl[ind], (radRRS_I[ind] .- radStupidFit)./radRRS_I[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
#plot(p1, p2,  layout = l)
 
#x_noRS0 = K \ radnoRS_I0[ind];

# COmpute fitted spectrum:
#radnoRS_I_fit0 = K * x_noRS0;


p1 = plot(wl[ind], radnoRS_I0[ind], label="Measured")
plot!(wl[ind], radnoRS_I_fit0, label="Fit without SIF")
p2 = plot(wl[ind], (radnoRS_I0[ind] .- radnoRS_I_fit0)./radnoRS_I0[ind]*1000, label="Residuals", xlabel="Wavelength (nm)", ylabel="Residuals  in permille")
plot(p1, p2,  layout = l)