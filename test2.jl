using Statistics, ImageFiltering, Plots, LinearAlgebra, PProf
# import PlotlyJS: surface

include("src/potentials.jl")
include("src/structs.jl")
include("src/utils.jl")
include("src/sim_utils.jl")
include("src/simulation.jl")
include("src/training_data.jl")

L = 10
dx = 0.1



Vext_gen = function()
	po = PotentialOptions(L=L, num_sin=(5, 8), sin_amp=(0, 1), sin_periods=(1, 3), periodic=true, wall=:none)
	generate_random_potential(po)
end



accept_condition = function(rho)
	# check if the density is too low
	if mean(rho) < 0
		return false
	end
	# check if the density is too high
	if mean(rho) > 0.5
		return false
	end
	srho = sort(rho |> vec)
	n = length(srho)÷100
	maxmin = mean(srho[end-n:end]) / mean(srho[1:n])
	if maxmin > 10
		return false
	end
	return true
end

function smooth_rho(rho, sys::GCMC_System)
	smooth_factor = 0.075 * sys.σ / sys.dx
	# smooth the density
	rho = imfilter(rho, Kernel.gaussian(smooth_factor))
	return rho
end

to = GCMC_TrainingData("data",
	L=L,
	Vext_generator=Vext_gen,
	steps=1e7,
	num_systems=100,
	repetitions=100,
	dx=dx,
	accept_condition=accept_condition,
	rho_smooth_func=smooth_rho,
	μ_range=(-5, 5),
)

# create_training_data(to, true)
function start()
	
	files = filter(x -> occursin(r"data.*\.dat", x), readdir(to.folder, join=true))
	accs = 0
	rejs = 0
	for file in files
		data = readdlm(file, ';')
		rho = data[:,2]
		rho = reshape(rho, (Int(sqrt(length(rho))), Int(sqrt(length(rho)))))
		
		if accept_condition(rho)
			accs += 1
			heatmap(rho)|>display
			sleep(1)
		else
			rejs += 1
		end
	end
	println("Accepted: ", accs, " Rejected: ", rejs)
end
# start()