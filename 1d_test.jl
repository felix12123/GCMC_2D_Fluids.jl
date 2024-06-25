using Statistics

include("src/potentials.jl")
include("src/structs.jl")
include("src/utils.jl")
include("src/sim_utils.jl")
include("src/simulation.jl")
include("src/training_data.jl")
include("src/training_data_1D.jl")


L = 10
d = 1
dx = 0.01


accept_condition = function(rho)
	# check if the density is too low
	# if mean(rho) < 0
	# 	return false
	# end
	# check if the density is too high
	if mean(rho) > 0.5
		return false
	end
	srho = filter(isfinite, sort(rho |> vec))
	n = length(srho)÷20
	maxmin = mean(srho[end-n:end]) / mean(srho[1:n])
	if maxmin > 20
		return false
	end
	return true
end

Vext_gen = function()
	po = PotentialOptions(L=L, num_sin=(5, 10), sin_amp=(0, 1), sin_periods=(0.5, 4), periodic=false, wall=:walls)
	generate_random_potential_1D(po)
end

to = GCMC_TrainingData("data",
	L=L,
	Vext_generator=Vext_gen,
	steps=1e6,
	num_systems=20,
	repetitions=10,
	dx=dx,
	accept_condition=accept_condition,
	# rho_smooth_func=smooth_rho,
	μ_range=(-5, 5),
	threads=(2, Threads.nthreads()-2),
)

create_training_data_1D(to, true)

