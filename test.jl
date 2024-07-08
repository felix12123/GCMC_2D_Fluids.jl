# include("src/GCMC_2D_Fluids.jl")
using Statistics, ImageFiltering, Plots, LinearAlgebra, PProf
import PlotlyJS: surface

include("src/potentials.jl")
include("src/structs.jl")
include("src/utils.jl")
include("src/sim_utils.jl")
include("src/simulation.jl")
include("src/training_data.jl")

# Create a system
L = 10
dx = 0.1
dx_inv = Int(1/dx)

Vext_gen = function()
	po = PotentialOptions(L=L, num_sin=(5, 8), sin_amp=(0, 1), sin_periods=(1, 3), periodic=true, wall=:box)
	generate_random_potential(po)
end

# Create and show Vext
V = Vext_gen()
x = Iterators.product(dx/2:dx:L, dx/2:dx:L)
v = V.(x)
heatmap(v) |> display

# simulate and show the result
sys = GCMC_System(L=L, dx=0.1, Vext=V, μ=0)
a = simulate(sys, 1e6 |> Int, 4e4|>Int, repetitions=5)[1]

heatmap(a) |> display