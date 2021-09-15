#=
 
This file contains complex error functions that are to be used in Voigt 
line-shape calculations. 

The `w(z)` functions near the bottom are those directly called by line_shape in 
compute_absorption_cross_section.jl. All other functions are helpers for 
calculating rational approximations. 
 
=#

"""
    Humlicek (1982) rational approximation:  region I.
"""
function humlicek1(z::Complex)
    FT = eltype(real(z));
    #t = z.imag-1j*z.real;   w = t * recSqrtPi / (0.5 + t*t)
    w = 1im*FT(1/sqrt(pi)) *z / (z*z-FT(0.5))
end

"""
    Humlicek (1982) rational approximation:  region II.
"""
function humlicek2(z::Complex)
    FT = eltype(real(z));
    # this implementation corresponds to the fortran code
    t = imag(z)-1im*real(z)
    u = t * t
    w = (t * (FT(1.410474) + u*FT(1/sqrt(pi))))/ (3/4 + (u *(3+u)))
end

"""
    Humlicek (1982) rational approximation:  region II.
"""
function _humlicek2(z::Complex)
    FT = eltype(real(z));
    # this implementation corresponds to Eq. (12) of the manuscript
    zz = z*z
    w = 1im* (z * (zz*FT(1/sqrt(pi))-FT(1.410474)))/ (3/4 + zz*(zz-3))
end

"""
    Humlicek (1982) rational approximation:  region III.
"""
function humlicek3(z::Complex)
    FT = eltype(real(z));
    t = imag(z)-1im*real(z)
    w =  (FT(16.4955) + t * (FT(20.20933) + t * (FT(11.96482) + t * (FT(3.778987) + FT(0.5642236)*t)))) /
         (FT(16.4955) + t * (FT(38.82363) + t * (FT(39.27121) + t * (FT(21.69274) + t * (FT(6.699398) + t)))))
end

"""
    Humlicek (1982) rational approximation:  region IV.
"""
function humlicek4(z::Complex)
    FT = eltype(real(z));
    t = imag(z)-1im*real(z)
    u = t*t
    nom = t*(FT(36183.31)-u*(FT(3321.99)-u*(FT(1540.787)-u*(FT(219.031)-u*(FT(35.7668)-u*(FT(1.320522)-u*FT(.56419)))))))
    den = FT(32066.6)-u*(FT(24322.8)-u*(FT(9022.23)-u*(FT(2186.18)-u*(FT(364.219)-u*(FT(61.5704)-u*(FT(1.84144)-u))))))
    w  = exp(u) - nom/den
end

"""
    Humlicek (1982) complex probability function:  w4 rational approximations.
"""
function humlicek(z::Complex)
    s = abs(real(z)) + imag(z)
    # Choices depending on regions:
    if s>15.0
        w = humlicek1(z)
    elseif s>5.5
        w = humlicek2(z)
    else
        if imag(z)>=0.195*abs(real(z))-0.176
            w = humlicek3(z)
        else
            w = humlicek4(z)
        end
    end
    return w
end


"""
    Humlicek (1979) complex probability function:  rational approximation for y>0.85 OR |x|<18.1*y+1.65.
"""
function cpf12a(z::Complex)
    FT = eltype(real(z));
    cr  = FT(1.5);
    crr = FT(2.25);
    #a = similar(z)
    ct1 = FT(.3142403762544);  ct2 = FT(.9477883912402);  ct3 = FT(1.5976826351526);
    ct4 = FT(2.2795070805011); ct5 = FT(3.0206370251209); ct6 = FT(3.88972489786978);
    
    ca1 = FT(-1.393236997981977);  ca2 = FT(-0.2311524061886763);  ca3 = FT(+0.1553514656420944);
    ca4 = FT(-0.006218366236965554); ca5 = FT(9.190829861057117e-5); ca6 = FT(+6.275259577e-7);
    
    cb1 = FT(1.011728045548831);  cb2 = FT(-0.7519714696746353);  cb3 = FT(0.01255772699323164);
    cb4 = FT(0.01002200814515897); cb5 = FT(-2.420681348155727e-4); cb6 = FT(5.008480613664576e-7);

    x = real(z)
    y = imag(z)
    ry = cr + y
    ryy = ry * ry
    wk =  (ca1*(x-ct1) + cb1*ry) / ((x-ct1)^2 + ryy) - (ca1*(x+ct1) - cb1*ry) / ((x+ct1)^2 + ryy) +
          (ca2*(x-ct2) + cb2*ry) / ((x-ct2)^2 + ryy) - (ca2*(x+ct2) - cb2*ry) / ((x+ct2)^2 + ryy) +
          (ca3*(x-ct3) + cb3*ry) / ((x-ct3)^2 + ryy) - (ca3*(x+ct3) - cb3*ry) / ((x+ct3)^2 + ryy) +
          (ca4*(x-ct4) + cb4*ry) / ((x-ct4)^2 + ryy) - (ca4*(x+ct4) - cb4*ry) / ((x+ct4)^2 + ryy) +
          (ca5*(x-ct5) + cb5*ry) / ((x-ct5)^2 + ryy) - (ca5*(x+ct5) - cb5*ry) / ((x+ct5)^2 + ryy) +
          (ca6*(x-ct6) + cb6*ry) / ((x-ct6)^2 + ryy) - (ca6*(x+ct6) - cb6*ry) / ((x+ct6)^2 + ryy);
    wl =  (cb1*(x-ct1) - ca1*ry) / ((x-ct1)^2 + ryy) + (cb1*(x+ct1) + ca1*ry) / ((x+ct1)^2 + ryy) +
          (cb2*(x-ct2) - ca2*ry) / ((x-ct2)^2 + ryy) + (cb2*(x+ct2) + ca2*ry) / ((x+ct2)^2 + ryy) +
          (cb3*(x-ct3) - ca3*ry) / ((x-ct3)^2 + ryy) + (cb3*(x+ct3) + ca3*ry) / ((x+ct3)^2 + ryy) +
          (cb4*(x-ct4) - ca4*ry) / ((x-ct4)^2 + ryy) + (cb4*(x+ct4) + ca4*ry) / ((x+ct4)^2 + ryy) +
          (cb5*(x-ct5) - ca5*ry) / ((x-ct5)^2 + ryy) + (cb5*(x+ct5) + ca5*ry) / ((x+ct5)^2 + ryy) +
          (cb6*(x-ct6) - ca6*ry) / ((x-ct6)^2 + ryy) + (cb6*(x+ct6) + ca6*ry) / ((x+ct6)^2 + ryy);
    a = wk + im * wl
    #println(a)
    return a#a   # wk, wl
end

"""
    Humlicek (1979) complex probability function:  rational approximation for y<0.85 AND |x|>18.1*y+1.65.
"""
function cpf12b(z::Complex)
    FT = eltype(real(z));
    ct1 = FT(.3142403762544);  ct2 = FT(.9477883912402);  ct3 = FT(1.5976826351526);
    ct4 = FT(2.2795070805011); ct5 = FT(3.0206370251209); ct6 = FT(3.88972489786978);
    
    ca1 = FT(-1.393236997981977);  ca2 = FT(-0.2311524061886763);  ca3 = FT(+0.1553514656420944);
    ca4 = FT(-0.006218366236965554); ca5 = FT(9.190829861057117e-5); ca6 = FT(+6.275259577e-7);
    
    cb1 = FT(1.011728045548831);  cb2 = FT(-0.7519714696746353);  cb3 = FT(0.01255772699323164);
    cb4 = FT(0.01002200814515897); cb5 = FT(-2.420681348155727e-4); cb6 = FT(5.008480613664576e-7);
    cr  = FT(1.5);
    crr = FT(2.25);

    x = real(z);
    y = imag(z);
    ry   = cr+y      # yy0   = y+1.5
    y2r  = y +2cr  # y2y0  = y+3.0
    rry  = cr*ry     # y0yy0 = 1.5*(y+1.5)
    ryry = ry^2     # yy0sq = (y+1.5)^2
    wk =  ( cb1*((x-ct1)^2-rry) - ca1*(x-ct1)*y2r ) / (((x-ct1)^2+crr)*((x-ct1)^2+ryry)) +
          ( cb1*((x+ct1)^2-rry) + ca1*(x+ct1)*y2r ) / (((x+ct1)^2+crr)*((x+ct1)^2+ryry)) +
          ( cb2*((x-ct2)^2-rry) - ca2*(x-ct2)*y2r ) / (((x-ct2)^2+crr)*((x-ct2)^2+ryry)) +
          ( cb2*((x+ct2)^2-rry) + ca2*(x+ct2)*y2r ) / (((x+ct2)^2+crr)*((x+ct2)^2+ryry)) +
          ( cb3*((x-ct3)^2-rry) - ca3*(x-ct3)*y2r ) / (((x-ct3)^2+crr)*((x-ct3)^2+ryry)) +
          ( cb3*((x+ct3)^2-rry) + ca3*(x+ct3)*y2r ) / (((x+ct3)^2+crr)*((x+ct3)^2+ryry)) +
          ( cb4*((x-ct4)^2-rry) - ca4*(x-ct4)*y2r ) / (((x-ct4)^2+crr)*((x-ct4)^2+ryry)) +
          ( cb4*((x+ct4)^2-rry) + ca4*(x+ct4)*y2r ) / (((x+ct4)^2+crr)*((x+ct4)^2+ryry)) +
          ( cb5*((x-ct5)^2-rry) - ca5*(x-ct5)*y2r ) / (((x-ct5)^2+crr)*((x-ct5)^2+ryry)) +
          ( cb5*((x+ct5)^2-rry) + ca5*(x+ct5)*y2r ) / (((x+ct5)^2+crr)*((x+ct5)^2+ryry)) +
          ( cb6*((x-ct6)^2-rry) - ca6*(x-ct6)*y2r ) / (((x-ct6)^2+crr)*((x-ct6)^2+ryry)) +
          ( cb6*((x+ct6)^2-rry) + ca6*(x+ct6)*y2r ) / (((x+ct6)^2+crr)*((x+ct6)^2+ryry))
    wl =  (cb1*(x-ct1) - ca1*ry) / ((x-ct1)^2 + ryry) + (cb1*(x+ct1) + ca1*ry) / ((x+ct1)^2 + ryry) +
          (cb2*(x-ct2) - ca2*ry) / ((x-ct2)^2 + ryry) + (cb2*(x+ct2) + ca2*ry) / ((x+ct2)^2 + ryry) +
          (cb3*(x-ct3) - ca3*ry) / ((x-ct3)^2 + ryry) + (cb3*(x+ct3) + ca3*ry) / ((x+ct3)^2 + ryry) +
          (cb4*(x-ct4) - ca4*ry) / ((x-ct4)^2 + ryry) + (cb4*(x+ct4) + ca4*ry) / ((x+ct4)^2 + ryry) +
          (cb5*(x-ct5) - ca5*ry) / ((x-ct5)^2 + ryry) + (cb5*(x+ct5) + ca5*ry) / ((x+ct5)^2 + ryry) +
          (cb6*(x-ct6) - ca6*ry) / ((x-ct6)^2 + ryry) + (cb6*(x+ct6) + ca6*ry) / ((x+ct6)^2 + ryry)
    a = exp(-x * x) + y * wk + im * wl   # wk, wl
    return a
end

"""
    Complex error function --- J.A.C. Weideman: SIAM J. Num. Anal. 31 1497-1518 (1994); equation (38.I) and table I.
"""
function weideman32a(z::Complex)
    FT = eltype(real(z));
    L32 = FT(sqrt(32/sqrt(2)));
    a1=  FT(2.5722534081245696e+00);  a2 = FT(2.2635372999002676e+00);  a3 = FT(1.8256696296324824e+00);  a4 = FT(1.3455441692345453e+00);
    a5=  FT(9.0192548936480144e-01);  a6 = FT(5.4601397206393498e-01);  a7 = FT(2.9544451071508926e-01);  a8 = FT(1.4060716226893769e-01);
    a9=  FT(5.7304403529837900e-02);  a10= FT(1.9006155784845689e-02);  a11= FT(4.5195411053501429e-03);  a12= FT(3.9259136070122748e-04);
    a13= FT(-2.4532980269928922e-04); a14= FT(-1.3075449254548613e-04); a15= FT(-2.1409619200870880e-05); a16= FT(6.8210319440412389e-06);
    a17= FT(4.4015317319048931e-06);  a18= FT(4.2558331390536872e-07);  a19= FT(-4.1840763666294341e-07); a20= FT(-1.4813078891201116e-07);
    a21= FT(2.2930439569075392e-08);  a22= FT(2.3797557105844622e-08);  a23= FT(8.1248960947953431e-10);  a24= FT(-3.2080150458594088e-09);
    a25= FT(-5.2310170266050247e-10); a26= FT(4.1537465934749353e-10);  a27= FT(1.1658312885903929e-10);  a28= FT(-5.5441820344468828e-11);
    a29= FT(-2.1542618451370239e-11); a30= FT(8.0314997274316680e-12);  a31= FT(3.7424975634801558e-12);  a32= FT(-1.3031797863050087e-12);

    iz  = 1im*real(z) - imag(z)
    lpiz    = L32 + iz  # wL - y + x*complex(0.,1.)
    lmiz    = L32 - iz  # wL + y - x*complex(0.,1.)
    recLmiZ = 1 / lmiz
    Z       = lpiz * recLmiZ
    polynom = a1+(a2+(a3+(a4+(a5+(a6+(a7+(a8+(a9+(a10+(a11+(a12+(a13+(a14+(a15+(a16+(a17+(a18+(a19+(a20+(a21+(a22+(a23+(a24+(a25+(a26+(a27+(a28+(a29+(a30+(a31+a32*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z
    w       = (FT(1/sqrt(pi))  +  2 * polynom * recLmiZ)  *  recLmiZ
    return w
end

"""
    Humlicek (1979) complex probability function w(z): a single rational approximation.
"""
function w(::CPF12ErrorFunction, z::Complex)
    FT = eltype(real(z));
    # Needs to be double checked, implemented slightly differently from Schreier
    if (abs(real(z))<FT(18.1)*imag(z)+FT(1.65)) || (imag(z)>FT(0.85))
        return cpf12a(z)
    else
        return cpf12b(z)
    end
end

"""
    Complex error function w(z) using Humlicek's and Weideman's rational approximation:
        |x|+y>15:  Humlicek (JQSRT, 1982) rational approximation for region I;
        else:      J.A.C. Weideman (SIAM-NA 1994); equation (38.I) and table I.
"""
function w(::HumlicekWeidemann32VoigtErrorFunction, z::Complex)
    
    FT = eltype(real(z));
    if abs(real(z))+imag(z)>15
        w = 1im*FT(1/sqrt(pi)) * z / (z*z-FT(0.5))  # Humlicek (1982) approx 1 for s>15
    else
        w = weideman32a(z)
    end
    return w
end

"""
 Complex error function w(z) using Humlicek's and Weideman's rational approximation:
        |x|+y>8:  Humlicek (JQSRT, 1982) rational approximation for region II;
        else:      J.A.C. Weideman (SIAM-NA 1994); equation (38.I) and table I.
"""
function w(::HumlicekWeidemann32SDErrorFunction, z::Complex)
    FT = eltype(real(z));
    if abs(real(z))+imag(z)>=8
        w = humlicek2(z)  # Humlicek (1982) approx 1 for s>8
    else
        w = weideman32a(z)
    end
    return w
end

"""
    Complex error function w(z) using Humlicek's and Weideman's rational approximation:
        |x|+y>15:  Humlicek (JQSRT, 1982) rational approximation for region I;
        else:      using erfc(-iz) from Special Functions
"""
function w(::ErfcHumliErrorFunctionVoigt, z::Complex)
    FT = eltype(real(z));
    if abs(real(z))+imag(z)>15
        w = 1im*FT(1/sqrt(pi)) * z / (z*z-FT(0.5))  # Humlicek (1982) approx 1 for s>15
    else
        w = erfcx(-1im*z)
    end
    return w
end

"""
 Complex error function w(z) using Humlicek's and Weideman's rational approximation:
        |x|+y>8:  Humlicek (JQSRT, 1982) rational approximation for region II;
        else:     using erfc(-iz) from Special Functions
"""
function w(::ErfcHumliErrorFunctionSD, z::Complex)
    FT = eltype(real(z));
    if abs(real(z))+imag(z)>=8
        w = humlicek2(z)  # Humlicek (1982) approx 1 for s>8
    else
        w = erfcx(-1im*z)
    end
    return w
end

"""
 Complex error function w(z) using erfc(-iz) from Special Functions
"""
function w(::ErfcErrorFunction, z::Complex)
    erfcx(-im*z)
end
