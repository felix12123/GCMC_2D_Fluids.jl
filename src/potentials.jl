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


Vext_sin(x; n::Int, A::Number, φ::Number, L::Number) = A * sin(2π * x * n / L + φ)

Vext_lin(x; x1::Number, x2::Number, E1::Number, E2::Number) = x > x1 && x < x2 ? E1 + (x - x1) * (E2 - E1) / (x2 - x1) : 0

Vext_wall(x; xw::Number, L::Number) = x < xw || x > L - xw ? Inf : 0

function generate_Vext(L::Number; num_sin=4, num_lin=rand(1:5), wall=true)
    Avar = 1.0
    sin_parameters = []
    for n in 1:num_sin  # Generate random parameters for periodic sine functions with increasing frequency
        push!(sin_parameters, (n = n, A = randn() * Avar, φ = rand() * 2π, L = L))
    end
    Evar = 1.0
    lin_parameters = []
    for _ in 1:num_lin  # Generate random parameters for discontinuous linear segments
        push!(lin_parameters, (x1 = round(rand() * L, digits=2), x2 = round(rand() * L, digits=2), E1 = randn() * Evar, E2 = randn() * Evar))
    end
    xwmax = 1.0
    wall_params = (xw = round(rand() * xwmax, digits=2), L = L)  # Set a random wall width
    function (xy)  # Return a method which evaluates a combination of all functions with the chosen parameters above
        result = 0.0
		  x = xy[1]
        for sin_params in sin_parameters
            result += Vext_sin(x; sin_params...)
        end
        for lin_params in lin_parameters
            result += Vext_lin(x; lin_params...)
        end
        if wall
            result += Vext_wall(x; wall_params...)
        end
        result
    end
end
