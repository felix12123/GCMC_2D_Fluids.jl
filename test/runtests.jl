using GCMC_2D_Fluids
import GCMC_2D_Fluids: try_move!, try_insert!, step!, try_delete!
using Test, Statistics


function intact_system(sys::GCMC_System)
    for i in 1:sys.N
        if sys.Vext(sys.positions[i, 1], sys.positions[i, 2]) == Inf
            return false
        end
        for j in 1:sys.N
            if i == j
                continue
            end
            if (sys.positions[i, 1] - sys.positions[j, 1])^2 + (sys.positions[i, 2] - sys.positions[j, 2])^2 < sys.σ^2
                return false
            end
        end
    end
    return true
end

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

function test_move()
    L=10
    dx=0.1
    V_point(x,y) = (x==5 && y==5) ? 0 : Inf
    sys = GCMC_System(L=L, σ=1.0, μ=0.0, β=1.0, Vext=V_point, dx=dx)

    # insert one particle
    insert!(sys, 5.0, 5.0)
    @test sys.N == 1
    @test sys.positions[1, 1] == 5.0
    @test sys.positions[1, 2] == 5.0

    move_flag = true
    for _ in 1:1000
        move_flag &= !try_move!(sys)
    end

    @test move_flag


    # insert two particles and try to move them
    V_circ(x,y) = (x-5)^2 + (y-5)^2 < 1 ? 0 : Inf
    sys = GCMC_System(L=L, σ=1.0, μ=0.0, β=1.0, Vext=V_circ, dx=dx)
    insert!(sys, 4.1, 5.0)
    insert!(sys, 5.9, 5.0)

    move_flag = true
    for i in 1:1000
        try_move!(sys)
        move_flag &= sys.Vext(sys.positions[1, 1], sys.positions[1, 2]) == 0.0
    end

    @test move_flag
end

function test_insert()
    L=10
    dx=0.1
    V(x,y) = (x-5)^2 + (y-5)^2 < 1 ? 0 : Inf
    sys = GCMC_System(L=L, σ=1.0, μ=0.0, β=1.0, Vext=V, dx=dx)

    insert_flag = true
    for _ in 1:1000
       
        for i in 1:100
            try_insert!(sys)
        end
        for i in 1:sys.N
            insert_flag &= V(sys.positions[1, 1], sys.positions[1, 2]) == 0.0
        end
        sys = GCMC_System(L=L, σ=1.0, μ=0.0, β=1.0, Vext=V, dx=dx)
    end
    @test insert_flag
end

function test_step()
    L=10
    dx=0.1
    V(x,y) = (x-5)^2 + (y-5)^2 < 1 ? 0 : Inf
    sys = GCMC_System(L=L, σ=1.0, μ=0.0, β=1.0, Vext=V, dx=dx)
    insert!(sys, 5.0, 5.0)

    step_flag = true
    for i in 1:1000
        step!(sys)
        if sys.N == 0
            continue
        else
            step_flag &= V(sys.positions[1, 1], sys.positions[1, 2]) == 0.0
        end
    end
    @test step_flag
end

function test_simulation()
    @testset "simulation" begin
        L = 5
        dx = 0.1
        xs = dx/2:dx:L
        xys = Iterators.product(xs, xs)
        function Vext(x, y)
            if x < 1 || x > L-1 || y < 1 || y > L-1
                return Inf
            end
            return 0.0
        end

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
            return min(V(xy1...), V(xy2...), V(xy3...), V(xy4...)) == Inf
        end

        v_is_inf = findall(cell_is_inf.(Ref(Vext), xys, dx))
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
    test_move()
    test_step()
    test_collisions()
    test_simulation()
end