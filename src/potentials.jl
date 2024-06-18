mutable struct PotentialOptions
	L::Float64
	num_sin::Tuple{Int, Int}
	sin_amp::Tuple{Float64, Float64}
	sin_periods::Tuple{Float64, Float64}
	periodic::Bool
	wall::Symbol
	wall_thickness::Tuple{Float64, Float64}

	function PotentialOptions(;L::Int, num_sin::Union{Int, Tuple{Int, Int}}=(1:5), sin_amp::Union{<:Real, Tuple{<:Real, <:Real}}=(0,1), sin_periods::Union{<:Real, Tuple{<:Real, <:Real}}=(1,3), periodic::Bool=true, wall::Symbol=:none, wall_thickness::Union{<:Real, Tuple{<:Real, <:Real}}=(0,1))
		if num_sin isa Int
			num_sin = (1, num_sin)
		end
		if sin_amp isa Real
			sin_amp = (0, sin_amp)
		end
		if sin_periods isa Real
			sin_periods = (1, sin_periods)
		end
		if periodic
			sin_periods = round.(Int, sin_periods)
			if sin_periods[1] < 1
				sin_periods[1] = 1
			end
			if sin_periods[2] < 1
				sin_periods[2] = 1
			end
		end
		if wall_thickness isa Real
			wall_thickness = (0, wall_thickness)
		end
		new(L, num_sin, sin_amp, sin_periods, periodic, wall, wall_thickness)
	end
end

function generate_random_potential(po::PotentialOptions)
	# generate sinus functions
	num_sin = rand(po.num_sin[1]:po.num_sin[2])
	sin_amps = rand(num_sin) .* (po.sin_amp[2] - po.sin_amp[1]) .+ po.sin_amp[1]
	if po.periodic
		sin_periods = rand(po.sin_periods[1]:po.sin_periods[2], num_sin)
	else
		sin_periods = rand(num_sin) .* (po.sin_periods[2] - po.sin_periods[1]) .+ po.sin_periods[1]
	end

	# roll directions of sinus wave, either x, y, or x-y direction. depends on periodicity
	directions = zeros(Int, num_sin, 2)
	if po.periodic
		for i in 1:num_sin
			(directions[i,:] .= rand([[1, 0], [0, 1]]))
		end
	else
		for i in 1:num_sin
			directions[i,:] .= rand([[1, 0], [0, 1], [1, 1]])
		end
	end

	# generate wall function
	if po.wall == :box || po.wall == :none
	else
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
		end
		return s
	end
end