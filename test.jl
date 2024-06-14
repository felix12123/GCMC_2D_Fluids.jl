# include("src/GCMC_2D_Fluids.jl")
using ProfileView, Statistics, ImageFiltering, Plots, LinearAlgebra

include("src/structs.jl")
include("src/utils.jl")
include("src/sim_utils.jl")
include("src/simulation.jl")

# Create a system
L = 10
dx = 0.1
dx_inv = Int(1/dx)

V2(xy) = (x=xy[1]; y=xy[2]; (x < 1 || x > L-1 || y < 1 || y > L-1) ? Inf : 0)
V1(xy) = 20*cos.(xy[2] / L * 2π) + cos.(pi/4 + 2*xy[2] / L * 2π) + sin.(xy[1] / L * 2π)
V3(xy) = 0
V(xy) = V1(xy)# + V2(xy)

sys = GCMC_System(dx=dx, Vext=V, μ=0.5, L=L, move_prob=0.9)

reps = 11*1*0+1
rs = zeros(dx_inv*L, dx_inv*L, reps)
Threads.@threads for i in 1:reps
	rs[:,:,i] .= simulate(sys, 1_500_000, 10_000)[1]
end
r = mean(rs, dims=3)[:,:,1]
heatmap(r, title="r") |> display

# smooth r wigh gauss kernel
kernel = KernelFactors.gaussian((0.5, 0.5))
rsmooth = imfilter(r, kernel, "circular")
heatmap(rsmooth, title="rsmooth") |> display

using DelimitedFiles
c1 = log.(rsmooth) .- (sys.μ .- V.(Iterators.product(dx:dx:L, dx:dx:L)))
heatmap(c1, title="c1") |> display
writedlm("data1.dat", [vec(c1) vec(rsmooth)], ';')
