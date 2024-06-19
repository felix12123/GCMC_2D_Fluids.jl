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
	for pos in sys.positions
		hist.ρ[floor(Int, pos[1]/hist.dx)+1, floor(Int, pos[2]/hist.dx)+1] += 1
	end
	hist.count += 1
end
function update_g_histogram!(hist::g_Histogram, sys::GCMC_System)
	rij::Float64 = 0.0
	for i in eachindex(sys.positions)
		for j in i+1:sys.N
			rij = norm(sys.positions[i] .- sys.positions[j])
			if floor(Int, rij/hist.dx)+1 > length(hist.ρ)
				continue
			end
			hist.ρ[floor(Int, rij/hist.dx)+1] += 1
		end
	end
	hist.count += 1
end

function simulate_once(sys::GCMC_System, steps::Int64, therm_steps::Int64, sample_interval::Int64=1000; obs::Vector=[], track_g=false)
	obs_vals = zeros(ceil(Int, steps/sample_interval), length(obs))
	hist = Histogram(sys)
	if track_g; g_hist = g_Histogram(sys); end

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

	g = track_g ? g_hist.ρ / g_hist.count : nothing
	return rho, g, obs_vals
end

function simulate(sys::GCMC_System, steps::Int64, therm_steps::Int64, sample_interval::Int64=1000; obs::Vector=[], repetitions::Int=1, threads::Int=Threads.nthreads(), track_g=false)
	if repetitions == 1
		return simulate_once(sys, steps, therm_steps, sample_interval; obs=obs)
	end

	# initialize variables
	obss = Vector{Any}(undef, repetitions)
	s = sys.L/sys.dx |> floor |> Int
	rhos = zeros(Float64, s, s, repetitions)
	gs = zeros(Float64, floor(Int, sys.L/2/sys.dx), repetitions)

	# create threads task packages
	tasks = 1:repetitions
	task_packages = [tasks[i:threads:end] for i in 1:threads]
	Threads.@threads for is in task_packages
		for i in is
			rhos[:, :, i], gs[:, i], obss[i] = simulate_once(deepcopy(sys), steps, therm_steps, sample_interval; obs=obs, track_g=false)
		end
	end
	
	rho = mean(rhos, dims=3)[:, :, 1]
	obs = mean(obss, dims=1)
	g = mean(gs, dims=2)[:, 1]
	return rho, g, obs
end