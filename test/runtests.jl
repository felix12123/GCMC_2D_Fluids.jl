using GCMC_2D_Fluids
import GCMC_2D_Fluids: try_move!, try_insert!
using Test, Statistics


function test_collisions()
    @testset "collisions" begin
        sys = GCMC_System(L=10, σ=1.0, μ=0.0, β=1.0, Vext=(x, y) -> 0.0)
        @test insert!(sys, 5.0, 0.0)
        @test sys.is_colliding(sys, 5.0, 0.0)
        @test sys.is_colliding(sys, 5.0, 0.4*sys.σ)
        @test sys.is_colliding(sys, 5.4, 0.0)
        @test sys.is_colliding(sys, 5.0, sys.L-sys.σ/3)
        @test !sys.is_colliding(sys, 5.0, sys.L-sys.σ*1.1)
        @test !sys.is_colliding(sys, 5.0, sys.σ*1.1)
    end
end

function test_simulation()
    @testset "simulation" begin
        L = 5
        dx = 0.1
        xs = dx/2:dx:L
        xys = Iterators.product(xs, xs)
        function Vext(xy)
            if any(xy .< 1) || any(xy .> L-1)
                return Inf
            end
            return 0.0
        end
        v = Vext.(xys)

        sys = GCMC_System(L=L, σ=1.0, μ=0.0, β=1.0, Vext=Vext, dx=dx)
        a = simulate(sys, 1e5 |> Int, 4e4|>Int, repetitions=5)[1]
        @test size(a) == (L/dx|>floor, L/dx|>floor)
        
        # determines if v is inf at all corners of a cell
        function cell_is_inf(V, xy, dx)
            x, y = xy
            xy1 = [x-dx/2, y-dx/2]
            xy2 = [x+dx/2, y-dx/2]
            xy3 = [x-dx/2, y+dx/2]
            xy4 = [x+dx/2, y+dx/2]
            return min(V(xy1), V(xy2), V(xy3), V(xy4)) == Inf
        end

        v_is_inf = findall(cell_is_inf.(Ref(Vext), xys, dx))
        println(mean(a[v_is_inf] .!= 0.0)) # check if the average density is zero all cells where Vext is inf everywhere
        @test all(a[v_is_inf] .== 0.0)

        # check if particles move or not move depending on try_move! output
        sys = GCMC_System(L=L, σ=1.0, μ=0.0, β=1.0, Vext=Vext, dx=dx)
        # insert one particle randomly
        for _ in 1:1000
            GCMC_2D_Fluids.try_insert!(sys)
            if sys.N > 0
                break
            end
        end
        flag = true
        for _ in 1:1000
            old_pos = copy(sys.positions[1, :])
            moved = try_move!(sys) 
            flag = flag && (moved == !all(old_pos .== sys.positions[1, :]))
            if !flag
                @test flag
                break
            end
            flag = flag && sys.N == 1
            if !flag
                @test flag
                break
            end
        end
        if flag
            @test flag
        end
    end
end

@testset "GCMC_2D_Fluids.jl" begin
    # Write your tests here.
    # test_collisions()
    test_simulation()
end