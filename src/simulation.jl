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
	ρ = zeros(size(hist.ρ))
	for pos in sys.positions
		hist.ρ[floor(Int, pos[1]/hist.dx)+1, floor(Int, pos[2]/hist.dx)+1] += 1
	end
	hist.count += 1
end

function simulate_once(sys::GCMC_System, steps::Int64, therm_steps::Int64, sample_interval::Int64=1000; obs::Vector=[])
	obs_vals = zeros(ceil(Int, steps/sample_interval), length(obs))
	hist = Histogram(sys, dx=sys.dx)
	for i in 1:therm_steps
		step!(sys)
	end
	for i in 1:steps
		step!(sys)
		if i % sample_interval == 0
			obs_vals[ceil(Int, i/sample_interval), :] = [o(sys) for o in obs]
			update_histogram!(hist, sys)
		end
	end
	return hist.ρ/hist.count/(sys.dx^2), obs_vals
end

function simulate(sys::GCMC_System, steps::Int64, therm_steps::Int64, sample_interval::Int64=1000; obs::Vector=[], repetitions::Int=1, threads::Int=Threads.nthreads())
	if repetitions == 1
		return simulate_once(sys, steps, therm_steps, sample_interval; obs=obs)
	end

	obss = Vector{Any}(undef, repetitions)
	s = sys.L/sys.dx |> floor |> Int
	rhos = zeros(Float64, s, s, repetitions)

	# create threads task packages
	tasks = 1:repetitions
	task_packages = [tasks[i:threads:end] for i in 1:threads]
	Threads.@threads for is in task_packages
		for i in is
			rhos[:, :, i], obss[i] = simulate_once(deepcopy(sys), steps, therm_steps, sample_interval; obs=obs)
		end
	end

	rho = mean(rhos, dims=3)[:,:,1]
	obs = mean(obss, dims=1)
	return rho, obs
end