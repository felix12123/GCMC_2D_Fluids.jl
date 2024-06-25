# we can create a 1D projection of training data with potentials that are constant in one dimension.
function generate_random_potential_1D(po::PotentialOptions)
	# generate sinus functions
	num_sin = rand(po.num_sin[1]:po.num_sin[2])
	sin_amps = rand(num_sin) .* (po.sin_amp[2] - po.sin_amp[1]) .+ po.sin_amp[1]
	if po.periodic
		sin_periods = rand(po.sin_periods[1]:po.sin_periods[2], num_sin)
	else
		sin_periods = rand(num_sin) .* (po.sin_periods[2] - po.sin_periods[1]) .+ po.sin_periods[1]
	end

	# roll directions of sinus wave, either x, y, or x-y direction. depends on periodicity
	directions = hcat(ones(Int, num_sin), zeros(Int, num_sin))

	# generate wall function
	if !(po.wall == :box || po.wall == :none || po.wall == :walls)
		@warn "Wall type \"$(po.wall)\" not implemented, no wall is being used."
	end

	wall_thickness = rand() * (po.wall_thickness[2] - po.wall_thickness[1]) + po.wall_thickness[1]
	phase_shift = rand(num_sin) .* 2pi
	function V(xy)
		s = 0.0
		for i in 1:num_sin
			s += sin(xy ⋅ directions[i, :] * 2pi / po.L * sin_periods[i] + phase_shift[i]) * sin_amps[i]
		end
		
		if po.wall == :box
			if xy[1] < wall_thickness || xy[1] > po.L-wall_thickness || xy[2] < wall_thickness || xy[2] > po.L-wall_thickness
				s += Inf
			end
		elseif po.wall == :walls
			if xy[1] < wall_thickness || xy[1] > po.L-wall_thickness
				s += Inf
			end
		end
		return s
	end
end

function create_training_data_1D(opt::GCMC_TrainingData, verbose::Bool=false)
	my_println(x...) = if verbose; println(x...); end

	# testing if the potential is 1D
	Vext = opt.Vext_generator()
	for _ in 1:20
		x = rand() * opt.L
		Vs = [Vext([x, rand() * opt.L]) for _ in 1:100]
		if !all(diff(filter(isfinite, Vs)) .≈ 0)
			error("The potential is not 1D")
		end
	end

	my_println("Creating 2D training data")
	create_training_data(opt, verbose)

	my_println("Projecting to 1D...")

	folder = opt.folder
	files = filter(x -> occursin(".dat", x), readdir(folder, join=true))

	for file in files
		data = readdlm(file, ';')
		c1 = data[:, 1]
		rho = data[:, 2]

		if length(c1) < opt.L^2 / opt.dx^2
			@warn "Data file $file is too short, skipping."
			continue
		end

		L = Int(sqrt(length(c1)))
		c1 = reshape(c1, L, L)
		rho = reshape(rho, L, L)

		# project to 1D
		c1_1D = mean(c1, dims=2)[:, 1]
		rho_1D = mean(rho, dims=2)[:, 1]

		# save the data
		data_1D = [vec(c1_1D) vec(rho_1D)]
		writedlm(file, data_1D, ';')
	end

	
end