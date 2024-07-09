using Statistics, ImageFiltering, Plots, LinearAlgebra, PProf
# import PlotlyJS: surface

include("../src/potentials.jl")
include("../src/structs.jl")
include("../src/utils.jl")
include("../src/sim_utils.jl")
include("../src/simulation.jl")
include("../src/training_data.jl")

L = 10
dx = 0.03



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
	num_systems=5,
	repetitions=100,
	dx=dx,
	threads=(2, Threads.nthreads()÷2),
	# accept_condition=accept_condition,
	# rho_smooth_func=smooth_rho,
	μ_range=(-5, 15),
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


to = GCMC_TrainingData("fwefwegweg",
               L=L,
               Vext_generator=Vext_gen,
               steps=1e7,
               num_systems=1,
               repetitions=1,
               dx=dx,
               threads=(1, 1),
               # accept_condition=accept_condition,
               # rho_smooth_func=smooth_rho,
               μ_range=(-5, 15),
)
if isdir("fwefwegweg")
	rm("fwefwegweg", recursive=true)
end
@time create_training_data(to, false)
if isdir("fwefwegweg")
	rm("fwefwegweg", recursive=true)
end

using ProfileView, Profile, BenchmarkTools
sys = GCMC_System(L=10, Vext=Vext_gen(), μ=10)
simulate_once(sys, 10^5, 1000, 1000)
Profile.clear()

@benchmark try_move!(sys) seconds=10 samples=100_000
# ProfileView.@profview sum(1:10)
ProfileView.@profview create_training_data(to, false)# war 40 sec mit unoptimized
