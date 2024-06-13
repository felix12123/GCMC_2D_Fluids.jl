using GCMC_2D_Fluids
using Test


function test_collisions()
    @testset "collisions" begin
        sys = GCMC_System(L=10, σ=1.0, μ=0.0, β=1.0, Vext=(x, y) -> 0.0)
        @test insert!(sys, (5.0, 0.0))
        @test sys.is_colliding(sys, (5.0, 0.0))
        @test sys.is_colliding(sys, (5.0, 0.4*sys.σ))
        @test sys.is_colliding(sys, (5.4, 0.0))
        @test sys.is_colliding(sys, (5.0, sys.L-sys.σ/3))
        @test !sys.is_colliding(sys, (5.0, sys.L-sys.σ*1.1))
        @test !sys.is_colliding(sys, (5.0, sys.σ*1.1))
    end
end

@testset "GCMC_2D_Fluids.jl" begin
    # Write your tests here.
    test_collisions()
end