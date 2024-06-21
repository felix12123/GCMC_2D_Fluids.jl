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

po = PotentialOptions(L=L, num_sin=(1, 1), sin_amp=(0, 1), sin_periods=(1, 3), periodic=true, wall=:none)
V5 = generate_random_potential(po)
xys = Iterators.product(dx:dx:L, dx:dx:L)
# V5s = V5.(xys)
# heatmap(V5s, title="V5s") |> display

gauss(x, mu, sig) = 1/sqrt(2π*sig^2) * exp(-1/2 * ((x-mu)/sig)^2)
V1(xy) = 20*cos.(xy[2] / L * 2π) + cos.(pi/4 + 2*xy[2] / L * 2π) + sin.(xy[1] / L * 2π)
V2(xy) = (x=xy[1]; y=xy[2]; (x < 1 || x > L-1 || y < 1 || y > L-1) ? Inf : 0)
V3(xy) = 0
# V4(xy) = -15exp(-((xy[1]-L/2)^2 + (xy[2]-L/2)^2)/0.25) - exp(-(xy[1]-L/2)^2/0.15) - 1exp(-(xy[1]-L/2)^2/0.5)
V4(xy) = -0.15gauss(xy[1], L/2, 0.15)*gauss(xy[2], L/2, 0.15) - # small pocket to trap one particle for crystalisation
			50gauss(xy[1], L/2, 5)*gauss(xy[2], L/2, 5) - # large pocket to trap multiple particles
			0.05gauss(xy[1], L/2, 0.15) # small burrow so the crystalisation is symmetrical to x axis
#  - 0gauss(xy[1], L/2, 0.5) - 0*7^2gauss(xy[1], L/2, 100)*gauss(xy[2], L/2, 100)
# V(xy) = V4(xy)# + V2(xy)
# V = generate_random_potential(po)

# V = generate_Vext(L; num_sin=4, num_lin=0, wall=false)
# v = V.(Iterators.product(dx:dx:L, dx:dx))
# plot(v, title="V") |> display
V=V4
sys = GCMC_System(dx=dx, Vext=V, μ=130.0, L=L, move_prob=0.9)

reps = -1+Threads.nthreads()*6
# r, g = @time simulate(sys, 2*10^6, 3*10^4, repetitions=reps)
@pprof simulate(sys, 2*10^6, 0)

# smooth r wigh gauss kernel
smooth_factor = 0.075 / dx * sys.σ
kernel = KernelFactors.gaussian((smooth_factor, smooth_factor))
rsmooth = imfilter(r, kernel, "circular")

using DelimitedFiles
c1 = log.(rsmooth) .- (sys.μ .- V.(Iterators.product(dx:dx:L, dx:dx:L)))


writedlm("data1.dat", [vec(c1) vec(rsmooth)], ';')

gr()
pltV  = heatmap(V.(Iterators.product(dx:dx:L, dx:dx:L)), title="V")
pltc1 = heatmap(c1, title="c1")
pltrs = heatmap(rsmooth, title="rsmooth")
pltr  = heatmap(r, title="r")
pltg  = plot(eachindex(g) * dx, g, title="g")
plt = plot(pltV, pltc1, pltg, pltrs, layout=(2,2), size=(1000,800))
savefig(plt, "plot1.png")

writedlm("data2.dat", [dx:dx:L vec(mean(c1, dims=2)) vec(mean(r, dims=2))], '\t')

nothing;