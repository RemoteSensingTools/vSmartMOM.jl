# Needs some warning if memory is getting too large !
FT = Float64
n = 20
nSpec = 10000

for i = 1:4
    # Create CPU Matrices:
    r⁻⁺_ = (randn(FT, n, n, nSpec));
    t⁺⁺_ = (randn(FT, n, n, nSpec));
    r⁺⁻_ = (randn(FT, n, n, nSpec));
    t⁻⁻_ = (randn(FT, n, n, nSpec));

    # And move to GPU as CuArray
    r⁻⁺ = CuArray(r⁻⁺_);
    t⁺⁺ = CuArray(t⁺⁺_);
    r⁺⁻ = CuArray(r⁺⁻_);
    t⁻⁻ = CuArray(t⁻⁻_);

    RadiativeTransfer.RTM.apply_D_matrix!(i, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
   
    RadiativeTransfer.RTM.apply_D_matrix!(i, r⁻⁺_, t⁺⁺_, r⁺⁻_, t⁻⁻_)


    @testset "GPU-CPU consistency $i" begin
        for (matGPU, matCPU) in ((r⁻⁺, r⁻⁺_),
                             (t⁺⁺, t⁺⁺_), 
                             (r⁺⁻, r⁺⁻_), 
                             (t⁻⁻, t⁻⁻_))
            @test Array(matGPU) ≈ matCPU    
        end
        if i == 1
            for (matGPU, matCPU) in ((r⁻⁺, r⁺⁻),
                             (t⁺⁺, t⁻⁻), 
                             (r⁻⁺_, r⁺⁻_),
                             (t⁺⁺_, t⁻⁻_))
                @test Array(matGPU) ≈ Array(matCPU)    
            end
        end
    end
end
