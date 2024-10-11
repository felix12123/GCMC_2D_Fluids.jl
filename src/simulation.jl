using Statistics

function step!(sys::GCMC_System)
	if rand() < sys.move_prob
		try_move!(sys)
	elseif rand() < sys.insert_prob
		try_insert!(sys)
	else
		try_delete!(sys)
	end
end

function update_histogram!(hist::Histogram, sys::GCMC_System)
	for i in 1:sys.N
		hist.ρ[ceil(Int, sys.positions[i,2]/hist.dx), ceil(Int, sys.positions[i,1]/hist.dx)] += 1
	end
	hist.count += 1
end
function update_g_histogram!(hist::g_Histogram, sys::GCMC_System)
	rij::Float64 = 0.0
	for i in 1:sys.N
		for j in i+1:sys.N
			rij = norm(sys.positions[i,:] .- sys.positions[j,:])
			if floor(Int, rij/hist.dx)+1 > length(hist.ρ)
				continue
			end
			hist.ρ[floor(Int, rij/hist.dx)+1] += 1
		end
	end
	hist.count += 1
end

function simulate_once(sys::GCMC_System, steps::Int64, therm_steps::Int64, sample_interval::Int64=1000; obs::Vector=[], track_g=false)
	@assert isa(sys.Vext, Function)

	obs_vals = zeros(ceil(Int, steps/sample_interval), length(obs))
	hist = Histogram(sys)
	g_hist = g_Histogram(sys)

	for _ in 1:therm_steps
		step!(sys)
	end
	for i in 1:steps
		step!(sys)
		if i % sample_interval == 0
			obs_vals[ceil(Int, i/sample_interval), :] = [o(sys) for o in obs]
			update_histogram!(hist, sys)
			if track_g; update_g_histogram!(g_hist, sys); end
		end
	end
	rho = hist.ρ/hist.count/(sys.dx^2)

	g = g_hist.ρ / g_hist.count
	return rho, g, obs_vals
end

function simulate(sys::GCMC_System, steps::Number, therm_steps::Number, sample_interval::Number=100; repetitions::Int=1, threads::Int=Threads.nthreads(), track_g=false)

	steps, therm_steps, sample_interval = ceil(Int, steps), ceil(Int, therm_steps), ceil(Int, sample_interval)
	if repetitions == 1
		if !isa(sys.Vext, Function)
			sys1 = deepcopy(sys)
			sys1.Vext = (x,y)->eval_pot(x, y, deepcopy(sys.Vext)::NamedTuple)
		end

		return simulate_once(sys1, steps, therm_steps, sample_interval; obs=[], track_g=track_g)
	end

	# initialize variables
	s = sys.L/sys.dx |> floor |> Int
	# rhos = zeros(Float64, s, s, repetitions)
	# gs = zeros(Float64, floor(Int, sys.L/2/sys.dx), repetitions)

	# create threads task packages
	tasks = 1:repetitions
	task_packages = [tasks[i:threads:end] for i in 1:threads]

	# deepcopy systems before looping and create results array
	results::Vector{Tuple{Matrix{Float64}, Vector{Float64}}} = [(zeros(Float64, s, s), zeros(Float64, floor(Int, sys.L/2/sys.dx))) for _ in 1:repetitions]
	Threads.@threads for n in 1:threads
		is = task_packages[n]
		for i in is
			system = deepcopy(sys)
			if !isa(system.Vext, Function)
				system.Vext = (x,y)->eval_pot(x, y, deepcopy(sys.Vext))
			end
			rho_i, g_i = simulate_once(system, steps, therm_steps, sample_interval; track_g=track_g)
			results[i][1] .= rho_i
			results[i][2] .= g_i
		end
	end

	rho = mean(first.(results))
	rho_uncertaincy = zeros(s, s)
	if repetitions > 1
		rho_uncertaincy = std(first.(results))
	end
	
	g = mean(last.(results))
	return rho, g, rho_uncertaincy
end