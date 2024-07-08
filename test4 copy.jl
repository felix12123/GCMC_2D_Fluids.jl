# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()

using Statistics, ImageFiltering, Plots, LinearAlgebra, PProf
import PlotlyJS: surface

include("src/potentials.jl")
include("src/structs.jl")
include("src/utils.jl")
include("src/sim_utils.jl")
include("src/simulation.jl")
include("src/training_data.jl")
include("src/training_data_1D.jl")

L=10
dx=0.01

Vext_gen = function()
	po = PotentialOptions(L=L, num_sin=(5, 8), sin_amp=(0.1, 0.5), sin_periods=(0.5, 3), periodic=false, wall=:walls)
	generate_random_potential_1D(po)
end

function accept_condition(rho)
   mrho = mean(filter(!iszero, rho))
   if mrho <= 0.0
      return false
   end
   if mrho > 1
      return false
   end
   return true
end


to = GCMC_TrainingData("dx001-1e9steps_1dpot",
               L=L,
               Vext_generator=Vext_gen,
               steps=1e7,
               num_systems=5,
               repetitions=100,
               dx=dx,
               threads=(2, Threads.nthreads()-1),
               accept_condition=accept_condition,
               # rho_smooth_func=smooth_rho,
               μ_range=(-5, 10),
)



create_training_data_1D(to, true)