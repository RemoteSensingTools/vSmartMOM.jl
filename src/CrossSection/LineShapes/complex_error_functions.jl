


"
    Humlicek (1982) rational approximation:  region I. 
"
function humlicek1(z)
    FT = eltype(real(z));
    #t = z.imag-1j*z.real;   w = t * recSqrtPi / (0.5 + t*t)
    w = 1im*FT(1/sqrt(pi)) *z / (z*z-FT(0.5))
end


"
    Humlicek (1982) rational approximation:  region II.
"
function humlicek2(z)
    FT = eltype(real(z));
    # this implementation corresponds to the fortran code
    t = imag(z)-1im*real(z)
    u = t * t
    w = (t * (FT(1.410474) + u*FT(1/sqrt(pi))))/ (3/4 + (u *(3+u)))
end

"
    Humlicek (1982) rational approximation:  region II.
"
function Humlicek2(z)
    FT = eltype(real(z));
    # this implementation corresponds to Eq. (12) of the manuscript
    zz = z*z
    w = 1im* (z * (zz*FT(1/sqrt(pi))-FT(1.410474)))/ (3/4 + zz*(zz-3))
end

"
    Humlicek (1982) rational approximation:  region III. 
"
function humlicek3(z)
    FT = eltype(real(z));
    t = imag(z)-1im*real(z)
    w =  (FT(16.4955) + t * (FT(20.20933) + t * (FT(11.96482) + t * (FT(3.778987) + FT(0.5642236)*t)))) /
         (FT(16.4955) + t * (FT(38.82363) + t * (FT(39.27121) + t * (FT(21.69274) + t * (FT(6.699398) + t)))))
end

"
    Humlicek (1982) rational approximation:  region IV. 
"
function humlicek4(z)
    FT = eltype(real(z));
    t = imag(z)-1im*real(z)
    u = t*t
    nom = t*(FT(36183.31)-u*(FT(3321.99)-u*(FT(1540.787)-u*(FT(219.031)-u*(FT(35.7668)-u*(FT(1.320522)-u*FT(.56419)))))))
    den = FT(32066.6)-u*(FT(24322.8)-u*(FT(9022.23)-u*(FT(2186.18)-u*(FT(364.219)-u*(FT(61.5704)-u*(FT(1.84144)-u))))))
    w  = exp(u) - nom/den
end

" 
    Humlicek (1982) complex probability function:  w4 rational approximations.
"
function humlicek(z)
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


####################################################################################################################################



cr=1.5; crr=2.25  #  y0 and y0^2 of the original
ct = np.array([0., .3142403762544,  .9477883912402, 1.5976826351526, 2.2795070805011, 3.0206370251209, 3.88972489786978])
ca = np.array([0., -1.393236997981977,    -0.2311524061886763,   +0.1553514656420944,
                   -0.006218366236965554, -9.190829861057117e-5, +6.275259577e-7])
cb = np.array([0.,  1.011728045548831,    -0.7519714696746353, 0.01255772699323164,
                    0.01002200814515897, -2.420681348155727e-4,  5.008480613664576e-7])

"
    Humlicek (1979) complex probability function:  rational approximation for y>0.85 OR |x|<18.1*y+1.65. 
"
function cpf12a(z)
    FT = eltype(real(z));
    ct  = Array{FT,1}([.3142403762544,  .9477883912402, 1.5976826351526, 2.2795070805011, 3.0206370251209, 3.88972489786978])
    ca  = Array{FT,1}([-1.393236997981977,  -0.2311524061886763,   +0.1553514656420944, -0.006218366236965554, -9.190829861057117e-5, +6.275259577e-7])
    cb  = Array{FT,1}([1.011728045548831,    -0.7519714696746353, 0.01255772699323164, 0.01002200814515897, -2.420681348155727e-4,  5.008480613664576e-7])
    cr  = FT(1.5);
    crr = FT(2.25);
    
    x = real(z)
    y = imag(z)
    ry = cr+y
    #ryry = ry^2
    wk =  (ca[1]*(x-ct[1]) + cb[1]*ry) / ((x-ct[1])^2 + ry^2) - (ca[1]*(x+ct[1]) - cb[1]*ry) / ((x+ct[1])^2 + ry^2) +
          (ca[2]*(x-ct[2]) + cb[2]*ry) / ((x-ct[2])^2 + ry^2) - (ca[2]*(x+ct[2]) - cb[2]*ry) / ((x+ct[2])^2 + ry^2) +
          (ca[3]*(x-ct[3]) + cb[3]*ry) / ((x-ct[3])^2 + ry^2) - (ca[3]*(x+ct[3]) - cb[3]*ry) / ((x+ct[3])^2 + ry^2) +
          (ca[4]*(x-ct[4]) + cb[4]*ry) / ((x-ct[4])^2 + ry^2) - (ca[4]*(x+ct[4]) - cb[4]*ry) / ((x+ct[4])^2 + ry^2) +
          (ca[5]*(x-ct[5]) + cb[5]*ry) / ((x-ct[5])^2 + ry^2) - (ca[5]*(x+ct[5]) - cb[5]*ry) / ((x+ct[5])^2 + ry^2) +
          (ca[6]*(x-ct[6]) + cb[6]*ry) / ((x-ct[6])^2 + ry^2) - (ca[6]*(x+ct[6]) - cb[6]*ry) / ((x+ct[6])^2 + ry^2);
    wl =  (cb[1]*(x-ct[1]) - ca[1]*ry) / ((x-ct[1])^2 + ry^2) + (cb[1]*(x+ct[1]) + ca[1]*ry) / ((x+ct[1])^2 + ry^2) +
          (cb[2]*(x-ct[2]) - ca[2]*ry) / ((x-ct[2])^2 + ry^2) + (cb[2]*(x+ct[2]) + ca[2]*ry) / ((x+ct[2])^2 + ry^2) +
          (cb[3]*(x-ct[3]) - ca[3]*ry) / ((x-ct[3])^2 + ry^2) + (cb[3]*(x+ct[3]) + ca[3]*ry) / ((x+ct[3])^2 + ry^2) +
          (cb[4]*(x-ct[4]) - ca[4]*ry) / ((x-ct[4])^2 + ry^2) + (cb[4]*(x+ct[4]) + ca[4]*ry) / ((x+ct[4])^2 + ry^2) +
          (cb[5]*(x-ct[5]) - ca[5]*ry) / ((x-ct[5])^2 + ry^2) + (cb[5]*(x+ct[5]) + ca[5]*ry) / ((x+ct[5])^2 + ry^2) +
          (cb[6]*(x-ct[6]) - ca[6]*ry) / ((x-ct[6])^2 + ry^2) + (cb[6]*(x+ct[6]) + ca[6]*ry) / ((x+ct[6])^2 + ry^2);
    return wk+1im*wl   # wk, wl
end

"
    Humlicek (1979) complex probability function:  rational approximation for y<0.85 AND |x|>18.1*y+1.65. 
"
function cpf12b(z)
    FT = eltype(real(z));
    ct  = Array{FT,1}([.3142403762544,  .9477883912402, 1.5976826351526, 2.2795070805011, 3.0206370251209, 3.88972489786978])
    ca  = Array{FT,1}([-1.393236997981977,  -0.2311524061886763,   +0.1553514656420944, -0.006218366236965554, -9.190829861057117e-5, +6.275259577e-7])
    cb  = Array{FT,1}([1.011728045548831,    -0.7519714696746353, 0.01255772699323164, 0.01002200814515897, -2.420681348155727e-4,  5.008480613664576e-7])
    cr  = FT(1.5);
    crr = FT(2.25);

    x = real(z);
    y = imag(z);
    ry   = cr+y      # yy0   = y+1.5
    y2r  = y +2cr  # y2y0  = y+3.0
    rry  = cr*ry     # y0yy0 = 1.5*(y+1.5)
    ryry = ry^2     # yy0sq = (y+1.5)^2
    wk =  ( cb[1]*((x-ct[1])^2-rry) - ca[1]*(x-ct[1])*y2r ) / (((x-ct[1])^2+crr)*((x-ct[1])^2+ryry)) +
          ( cb[1]*((x+ct[1])^2-rry) + ca[1]*(x+ct[1])*y2r ) / (((x+ct[1])^2+crr)*((x+ct[1])^2+ryry)) +
          ( cb[2]*((x-ct[2])^2-rry) - ca[2]*(x-ct[2])*y2r ) / (((x-ct[2])^2+crr)*((x-ct[2])^2+ryry)) +
          ( cb[2]*((x+ct[2])^2-rry) + ca[2]*(x+ct[2])*y2r ) / (((x+ct[2])^2+crr)*((x+ct[2])^2+ryry)) +
          ( cb[3]*((x-ct[3])^2-rry) - ca[3]*(x-ct[3])*y2r ) / (((x-ct[3])^2+crr)*((x-ct[3])^2+ryry)) +
          ( cb[3]*((x+ct[3])^2-rry) + ca[3]*(x+ct[3])*y2r ) / (((x+ct[3])^2+crr)*((x+ct[3])^2+ryry)) +
          ( cb[4]*((x-ct[4])^2-rry) - ca[4]*(x-ct[4])*y2r ) / (((x-ct[4])^2+crr)*((x-ct[4])^2+ryry)) +
          ( cb[4]*((x+ct[4])^2-rry) + ca[4]*(x+ct[4])*y2r ) / (((x+ct[4])^2+crr)*((x+ct[4])^2+ryry)) +
          ( cb[5]*((x-ct[5])^2-rry) - ca[5]*(x-ct[5])*y2r ) / (((x-ct[5])^2+crr)*((x-ct[5])^2+ryry)) +
          ( cb[5]*((x+ct[5])^2-rry) + ca[5]*(x+ct[5])*y2r ) / (((x+ct[5])^2+crr)*((x+ct[5])^2+ryry)) +
          ( cb[6]*((x-ct[6])^2-rry) - ca[6]*(x-ct[6])*y2r ) / (((x-ct[6])^2+crr)*((x-ct[6])^2+ryry)) +
          ( cb[6]*((x+ct[6])^2-rry) + ca[6]*(x+ct[6])*y2r ) / (((x+ct[6])^2+crr)*((x+ct[6])^2+ryry))
    wl =  (cb[1]*(x-ct[1]) - ca[1]*ry) / ((x-ct[1])^2 + ryry) + (cb[1]*(x+ct[1]) + ca[1]*ry) / ((x+ct[1])^2 + ryry) +
          (cb[2]*(x-ct[2]) - ca[2]*ry) / ((x-ct[2])^2 + ryry) + (cb[2]*(x+ct[2]) + ca[2]*ry) / ((x+ct[2])^2 + ryry) +
          (cb[3]*(x-ct[3]) - ca[3]*ry) / ((x-ct[3])^2 + ryry) + (cb[3]*(x+ct[3]) + ca[3]*ry) / ((x+ct[3])^2 + ryry) +
          (cb[4]*(x-ct[4]) - ca[4]*ry) / ((x-ct[4])^2 + ryry) + (cb[4]*(x+ct[4]) + ca[4]*ry) / ((x+ct[4])^2 + ryry) +
          (cb[5]*(x-ct[5]) - ca[5]*ry) / ((x-ct[5])^2 + ryry) + (cb[5]*(x+ct[5]) + ca[5]*ry) / ((x+ct[5])^2 + ryry) +
          (cb[6]*(x-ct[6]) - ca[6]*ry) / ((x-ct[6])^2 + ryry) + (cb[6]*(x+ct[6]) + ca[6]*ry) / ((x+ct[6])^2 + ryry)
    return exp(-x*x)+y*wk+1im*wl   # wk, wl
end

"
    Humlicek (1979) complex probability function: a single rational approximation. 
"
function cpf12(z)
    FT = eltype(real(z));
    # Needs to be double checked, implemented slightly differently from Schreier
    if (abs(real(z))<FT(18.1)*imag(z)+FT(1.65)) || (imag(z)>FT(0.85))
        return cpf12a(z)
    else
        return cpf12b(z)
    end
end


 #  sqrt(N)/2^(1/4) = 4.1195342878142354

a32 = np.array([ 2.6837530678616557e+00,  # a0 = L/sqrtPi
                 2.5722534081245696e+00,   2.2635372999002676e+00,   1.8256696296324824e+00,   1.3455441692345453e+00,
                 9.0192548936480144e-01,   5.4601397206393498e-01,   2.9544451071508926e-01,   1.4060716226893769e-01,
                 5.7304403529837900e-02,   1.9006155784845689e-02,   4.5195411053501429e-03,   3.9259136070122748e-04,
                -2.4532980269928922e-04,  -1.3075449254548613e-04,  -2.1409619200870880e-05,   6.8210319440412389e-06,
                 4.4015317319048931e-06,   4.2558331390536872e-07,  -4.1840763666294341e-07,  -1.4813078891201116e-07,
                 2.2930439569075392e-08,   2.3797557105844622e-08,   8.1248960947953431e-10,  -3.2080150458594088e-09,
                -5.2310170266050247e-10,   4.1537465934749353e-10,   1.1658312885903929e-10,  -5.5441820344468828e-11,
                -2.1542618451370239e-11,   8.0314997274316680e-12,   3.7424975634801558e-12,  -1.3031797863050087e-12])

"
    Complex error function --- J.A.C. Weideman: SIAM J. Num. Anal. 31 1497-1518 (1994); equation (38.I) and table I. 
"
function weideman32a(z)
    FT = eltype(real(z));
    L32 = FT(sqrt(32/sqrt(2)));
    a32 = Array{FT,1}([ 2.5722534081245696e+00,   2.2635372999002676e+00,   1.8256696296324824e+00,   1.3455441692345453e+00,
                        9.0192548936480144e-01,   5.4601397206393498e-01,   2.9544451071508926e-01,   1.4060716226893769e-01,
                        5.7304403529837900e-02,   1.9006155784845689e-02,   4.5195411053501429e-03,   3.9259136070122748e-04,
                       -2.4532980269928922e-04,  -1.3075449254548613e-04,  -2.1409619200870880e-05,   6.8210319440412389e-06,
                        4.4015317319048931e-06,   4.2558331390536872e-07,  -4.1840763666294341e-07,  -1.4813078891201116e-07,
                        2.2930439569075392e-08,   2.3797557105844622e-08,   8.1248960947953431e-10,  -3.2080150458594088e-09,
                       -5.2310170266050247e-10,   4.1537465934749353e-10,   1.1658312885903929e-10,  -5.5441820344468828e-11,
                       -2.1542618451370239e-11,   8.0314997274316680e-12,   3.7424975634801558e-12,  -1.3031797863050087e-12])

    x = real(z);
    y = imag(z);

    iz  = 1im*real(z) - imag(z)
    lpiz    = L32 + iz  # wL - y + x*complex(0.,1.)
    lmiz    = L32 - iz  # wL + y - x*complex(0.,1.)
    recLmiZ = 1 / lmiz
    Z       = lpiz * recLmiZ
    polynom = a32[1]+(a32[2]+(a32[3]+(a32[4]+(a32[5]+(a32[6]+(a32[7]+(a32[8]+(a32[9]+(a32[10]+(a32[11]+(a32[12]+(a32[13]+(a32[14]+(a32[15]+(a32[16]+(a32[17]+(a32[18]+(a32[19]+(a32[20]+(a32[21]+(a32[22]+(a32[23]+(a32[24]+(a32[25]+(a32[26]+(a32[27]+(a32[28]+(a32[29]+(a32[30]+(a32[31]+a32[32]*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z)*Z
    w       = (recSqrtPi  +  2 * polynom * recLmiZ)  *  recLmiZ
    return w
end

" Complex error function using Humlicek's and Weideman's rational approximation:
        |x|+y>15:  Humlicek (JQSRT, 1982) rational approximation for region I;
        else:      J.A.C. Weideman (SIAM-NA 1994); equation (38.I) and table I. 
"
function hum1wei32a(z)
    FT = eltype(real(z));
    if abs(real(z))+imag(z)>15
        w = 1im*FT(1/sqrt(pi)) * z / (z*z-FT(0.5))  # Humlicek (1982) approx 1 for s>15
    else
        w = weideman32a(z)
    end
    return w
end

" Complex error function using Humlicek's and Weideman's rational approximation:
        |x|+y>8:  Humlicek (JQSRT, 1982) rational approximation for region II;
        else:      J.A.C. Weideman (SIAM-NA 1994); equation (38.I) and table I. 
"
function hum2wei32a(z)
    FT = eltype(real(z));
    if abs(real(z))+imag(z)>=8
        w = humlicek2(z)  # Humlicek (1982) approx 1 for s>8
    else
        w = weideman32a(z)
    end
    return w
end

