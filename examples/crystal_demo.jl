using Statistics, ImageFiltering, Plots, LinearAlgebra, PProf
import PlotlyJS: surface

include("../src/potentials.jl")
include("../src/structs.jl")
include("../src/utils.jl")
include("../src/sim_utils.jl")
include("../src/simulation.jl")
include("../src/training_data.jl")


gauss(x, mu, sig) = 1/sqrt(2π*sig^2) * exp(-1/2 * ((x-mu)/sig)^2)
V(xy) = -0.15gauss(xy[1], L/2, 0.15)*gauss(xy[2], L/2, 0.15) - # small pocket to trap one particle for crystalisation
			100gauss(xy[1], L/2, 5)*gauss(xy[2], L/2, 5) - # large pocket to trap multiple particles
			0.15gauss(xy[1], L/2, 0.1) # small burrow so the crystalisation is symmetrical to x axis


L = 10
dx = 0.1
dx_inv = Int(1/dx)

sys = GCMC_System(dx=dx, Vext=V, μ=130.0, L=L, move_prob=0.9)

reps = -1+Threads.nthreads()*4

r, g, _ = simulate(sys, 2*10^6, 30_000, track_g=true, repetitions=reps)


using DelimitedFiles
c1 = log.(r) .- (sys.μ .- V.(Iterators.product(dx:dx:L, dx:dx:L)))


writedlm("crystal.dat", [vec(c1) vec(r)], ';')

gr()
pltV  = heatmap(V.(Iterators.product(dx:dx:L, dx:dx:L)), title="V")
pltc1 = heatmap(c1, title="c1")
pltrs = heatmap(r, title="r")
pltr  = heatmap(r, title="r")
pltg  = plot(eachindex(g) * dx, g, title="g")
plt = plot(pltV, pltc1, pltg, pltrs, layout=(2,2), size=(1000,800))
savefig(plt, "crystal.png")

# writedlm("data2.dat", [dx:dx:L vec(mean(c1, dims=2)) vec(mean(r, dims=2))], '\t')

nothing;

