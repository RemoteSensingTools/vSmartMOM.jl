#***************************
#...LaGrange 3- and 4-point interpolation
#...arrays A and B are the npt data points,  given aa, a value of the
#...A variable, the routine will find the corresponding bb value
#
#...input:  aa
#...output: bb
using Random

function AtoB(aa,A,B,npt)
  for I in 2:(npt+1)

    if A[I-1] >= aa
      if I < 3 || I == npt
        J = I
        if I < 3
          J = 3
        end
        if I == npt
          J = npt
        end

        # J = J - 1
        A0D1 = A[J-2] - A[J-1]

        if A0D1 == 0.0
          A0D1=0.0001
        end

        A0D2=A[J-2]-A[J]

        if A0D2 == 0.0
          A0D2=0.0000
        end

        A1D1=A[J-1]-A[J-2]

        if A1D1 == 0.0
          A1D1=0.0001
        end

        A1D2=A[J-1]-A[J]

        if A1D2 == 0.0
          A1D2=0.0001
        end

        A2D1=A[J]-A[J-2]

        if A2D1 == 0.0
          A2D1=0.0001
        end

        A2D2=A[J]-A[J-1]

        if A2D2 == 0.0
          A2D2=0.0001
        end

        A0=(aa-A[J-1])*(aa-A[J])/(A0D1*A0D2)
        A1=(aa-A[J-2])*(aa-A[J])/(A1D1*A1D2)
        A2=(aa-A[J-2])*(aa-A[J-1])/(A2D1*A2D2)

        bb = A0*B[J-2] + A1*B[J-1] + A2*B[J]

        return bb

      else

        J = I
        # J = J-1   # zero index correction
        A0D1=A[J-2]-A[J-1]

        if A0D1 == 0.0
          A0D1=0.0001
        end
        A0D2=A[J-2]-A[J]
        if A0D2 == 0.0
          A0D2=0.0001
        end
        A0D3 = (A[J-2]-A[J+1])

        if A0D3 == 0.0
          A0D3=0.0001
        end
        A1D1=A[J-1]-A[J-2]
        if A1D1 == 0.0
          A1D1=0.0001
        end
        A1D2=A[J-1]-A[J]
        if A1D2 == 0.0
          A1D2=0.0001
        end
        A1D3 = A[J-1]-A[J+1]

        if A1D3 == 0.0
          A1D3=0.0001
        end

        A2D1=A[J]-A[J-2]
        if A2D1 == 0.0
          A2D1=0.0001
        end
        A2D2=A[J]-A[J-1]
        if A2D2 == 0.0
          A2D2=0.0001
        end
        A2D3 = A[J]-A[J+1]
        if A2D3 == 0.0
          A2D3=0.0001
        end

        A3D1 = A[J+1]-A[J-2]
        if A3D1 == 0.0
          A3D1=0.0001
        end
        A3D2 = A[J+1]-A[J-1]
        if A3D2 == 0.0
          A3D2=0.0001
        end
        A3D3 = A[J+1]-A[J]
        if A3D3 == 0.0
          A3D3=0.0001
        end

        A0=(aa-A[J-1])*(aa-A[J])*(aa-A[J+1])
        A0=A0/(A0D1*A0D2*A0D3)
        A1=(aa-A[J-2])*(aa-A[J])*(aa-A[J+1])
        A1=A1/(A1D1*A1D2*A1D3)
        A2=(aa-A[J-2])*(aa-A[J-1])*(aa-A[J+1])
        A2=A2/(A2D1*A2D2*A2D3)
        A3=(aa-A[J-2])*(aa-A[J-1])*(aa-A[J])
        A3=A3/(A3D1*A3D2*A3D3)

        bb = A0*B[J-2] + A1*B[J-1] + A2*B[J] + A3*B[J+1]

        return bb

      end

      break

    end

  end

  return bb

end

function qoft(M, I, T; TT_=TIPS_2017_ISOT_HASH_CONST, TQ_ = TIPS_2017_ISOQ_HASH_CONST)

  # return 1
    # println("Recompiled 4")
    # get temperature grid
    # TT = TIPS_2017_ISOT_HASH[(M,I)]
    TT = TT_[(M,I)]
    TQ = TQ_[(M,I)]

    Tmin = minimum(TT); Tmax = maximum(TT)

    # out of temperature range
    if T<Tmin || T>Tmax
      error("TIPS2017: T ($T) must be between $Tmin K and $Tmax K.")
    end

    # interpolate partition sum for specified isotopologue

    # Lagrange interpolation from HAPI for now -- will eventually replace with
    # Interpolations.jl function

    Qt = AtoB(T,TT,TQ,length(TT))

    # You can uncomment this print line and see that the calculation is indeed
    # *Happening!*

    # println(rand() * Qt)

    return Qt

    # return rand() * Qt
    # return Qt

end

function qoft_ratio(M, I, T,T_ref; TT_=TIPS_2017_ISOT_HASH_CONST, TQ_ = TIPS_2017_ISOQ_HASH_CONST)
  # return 1
    # println("Recompiled 4")
    # get temperature grid
    # TT = TIPS_2017_ISOT_HASH[(M,I)]
    TT = TT_[(M,I)]
    TQ = TQ_[(M,I)]
    Tmin = minimum(TT); Tmax = maximum(TT)
    # out of temperature range
    if T<Tmin || T>Tmax
      error("TIPS2017: T ($T) must be between $Tmin K and $Tmax K.")
    end
    # interpolate partition sum for specified isotopologue
    # Lagrange interpolation from HAPI for now -- will eventually replace with
    # Interpolations.jl function
    Qt = AtoB(T,TT,TQ,length(TT))
    Qt2 = AtoB(T_ref,TT,TQ,length(TT))

    # println(Qt2/Qt);
    # result = Qt2/Qt;
    # You can uncomment this print line and see that the calculation is indeed
    # *Happening!*
    # println(rand() * Qt)
    return Qt2/Qt
    # return rand() * Qt
    # return Qt
end

function qoft!(M, I, T,T_ref, result; TT_=TIPS_2017_ISOT_HASH_CONST, TQ_ = TIPS_2017_ISOQ_HASH_CONST)
  # return 1
    # println("Recompiled 4")
    # get temperature grid
    # TT = TIPS_2017_ISOT_HASH[(M,I)]
    TT = TT_[(M,I)]
    TQ = TQ_[(M,I)]
    Tmin = minimum(TT); Tmax = maximum(TT)
    # out of temperature range
    if T<Tmin || T>Tmax
      error("TIPS2017: T ($T) must be between $Tmin K and $Tmax K.")
    end
    # interpolate partition sum for specified isotopologue
    # Lagrange interpolation from HAPI for now -- will eventually replace with
    # Interpolations.jl function
    Qt = AtoB(T,TT,TQ,length(TT))
    Qt2 = AtoB(T_ref,TT,TQ,length(TT))

    # println(Qt2/Qt);
    # result = Qt2/Qt;
    # You can uncomment this print line and see that the calculation is indeed
    # *Happening!*
    # println(rand() * Qt)
    result[1] = Qt2/Qt
    return nothing
    # return rand() * Qt
    # return Qt
end
