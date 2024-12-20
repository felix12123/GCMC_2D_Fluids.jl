using ProgressMeter, DelimitedFiles, Base.Threads

function get_c1(sys::GCMC_System, rho::Matrix{<:Real})::Matrix
	dx = sys.dx
	L = sys.L
	xs = Iterators.product(dx/2:dx:L, dx/2:dx:L)
	log.(rho) .- (sys.μ .- map(x->sys.Vext(x[2], x[1]), xs)) .* sys.β
end
function get_c1(sys::GCMC_System, rho::Vector{<:Real})::Vector
	dx = sys.dx
	L = sys.L
	xs = Iterators.product(dx/2:dx:L, dx/2:dx:L)
	log.(rho) .- (sys.μ .- map(x->sys.Vext(x[1], x[2]), xs)) .* sys.β
end

function create_training_data(opt::GCMC_TrainingData, verbose::Bool=false)
	data_folder = opt.folder
	N = opt.num_systems

	my_println(x...) = if verbose; println(x...); end
	my_println("Creating training data")
	if !isdir(data_folder)
		my_println("Creating folder: ", data_folder)
		mkdir(data_folder)
	end
	if !isdir("$data_folder/uncertainty")
		my_println("Creating folder: ", "$data_folder/uncertainty")
		mkdir("$data_folder/uncertainty")
	end


	acceptances = Atomic{Int}(0)
	rejections = Atomic{Int}(0)

	stop_computation = Threads.Atomic{Bool}(false)
	files = filter(x -> occursin(r"data.*\.dat", x), readdir(data_folder, join=true))
	file_number = Threads.Atomic{Int}(length(files)+1)
	if verbose; P = Progress(N, desc="Creating Training Data"); end

	tasks = 1:N*50
	task_bundels = [tasks[i:opt.threads[1]:end] for i in 1:opt.threads[1]]
	
	verbose && update!(P, 0, showvalues=[(:acc, acceptances[]), (:rej, rejections[])])
	Threads.@threads for ns in task_bundels
	# for ns in task_bundels
		for _ in ns
			if file_number[] > N # if we have enough files,
				stop_computation[] = true
				verbose && println(" we have ", count(x -> occursin(r"data.*.dat", x), readdir(data_folder)), " files, which is enough.")
				break
			end


			Vparams = generate_random_potential(opt.pot_opts)
			Vext = (x,y)->eval_pot(x, y, Vparams)
		
			μ = rand()*(opt.μ_range[2]-opt.μ_range[1]) + opt.μ_range[1]

			sys = GCMC_System(L=opt.L, σ=opt.σ, μ=μ, β=opt.β, Vext=Vext, mobility=opt.mobility, move_prob=opt.move_prob, insert_prob=opt.insert_prob, dx=opt.dx)

			rho, _, rho_uncertainty = simulate(sys, opt.steps, opt.therm_steps, opt.sample_interval, track_g=false, repetitions=opt.repetitions, threads=opt.threads[2])
			rho = opt.rho_smooth_func(rho, sys)

			if opt.accept_condition(rho)
				atomic_add!(acceptances, 1)
			else
				atomic_add!(rejections, 1)
				verbose ? update!(P, showvalues=[(:acc, acceptances[]), (:rej, rejections[])]) : nothing
				continue
			end

			# save the data
			c1 = get_c1(sys, rho)
			data = [vec(c1) vec(rho)]
			file_number_tmp = Threads.atomic_add!(file_number, 1)
			writedlm("$data_folder/data$(file_number_tmp).dat", data, ';')
			writedlm("$data_folder/uncertainty/uncert$(file_number_tmp).dat", vec(rho_uncertainty), ';')
			verbose ? next!(P, showvalues=[(:acc, acceptances[]), (:rej, rejections[])]) : nothing
		end
	end
	my_println("Training data created.")
	my_println("Acceptance rate: ", acceptances[]/(acceptances[]+rejections[]))

	data_files = filter(x -> occursin(r"data.*\.dat", x), readdir(data_folder, join=true))
	if length(data_files) > N
		println("Warning: More files than systems. Deleting $(length(data_files) - N) files.")
		for num_file in N+1:length(data_files)
			file = data_folder * "/data" * string(num_file) * ".dat"
			rm(file)
		end
	end
end