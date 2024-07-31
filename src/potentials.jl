mutable struct PotentialOptions
	L::Float64
	dx::Float64
	num_sin::Tuple{Int, Int}
	sin_amp::Tuple{Float64, Float64}
	sin_periods::Tuple{Float64, Float64}
	periodic::Bool
	wall::Symbol
	wall_thickness::Tuple{Float64, Float64}
	num_plateaus::Tuple{Int, Int}
	plateau_size::Tuple{Float64, Float64}
	plateau_potential::Tuple{Float64, Float64}

	function PotentialOptions(;L::Int, num_sin::Union{Int, Tuple{Int, Int}}=(1, 5), sin_amp::Union{<:Real, Tuple{<:Real, <:Real}}=(0,1), sin_periods::Union{<:Real, Tuple{<:Real, <:Real}}=(1,3), periodic::Bool=true, wall::Symbol=:none, wall_thickness::Union{<:Real, Tuple{<:Real, <:Real}}=(0,1), dx::Number=0.1, num_plateaus::Union{Int, Tuple{Int, Int}}=(0, 0), plateau_size::Union{<:Real, Tuple{<:Real, <:Real}}=(0,L/2), plateau_potential::Union{<:Real, Tuple{<:Real, <:Real}}=(-1,1))
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
			wall_thickness = (wall_thickness/2, wall_thickness)
		end
		if num_plateaus isa Int
			num_plateaus = (0, num_plateaus)
		end
		if plateau_size isa Real
			plateau_size = (0, plateau_size)
		end
		if plateau_potential isa Real
			plateau_potential = (-plateau_potential, plateau_potential)
		end
		new(L, dx, num_sin, sin_amp, sin_periods, periodic, wall, wall_thickness, 	num_plateaus, plateau_size, plateau_potential)
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
			directions[i,:] .= rand([[1, 0], [0, 1]])
		end
	else
		for i in 1:num_sin
			directions[i,:] .= rand([[1, 0], [0, 1], [1, 1]])
		end
	end
	phase_shift = rand(num_sin) .* 2pi

	# generate wall function
	if po.wall == :box || po.wall == :none
	else
		@warn "Wall type \"$(po.wall)\" not implemented, no wall is being used."
	end
	wall_thickness = rand() * (po.wall_thickness[2] - po.wall_thickness[1]) + po.wall_thickness[1]
	wall_thickness = max(wall_thickness, wall_thickness % po.dx)

	num_plateaus = rand(po.num_plateaus[1]:po.num_plateaus[2])
	plateau_height = rand(num_plateaus) .* (po.plateau_size[2] - po.plateau_size[1]) .+ po.plateau_size[1]
	plateau_width = rand(num_plateaus) .* (po.plateau_size[2] - po.plateau_size[1]) .+ po.plateau_size[1]
	plateau_corner = [rand(2) .* [po.L .- plateau_width[i], po.L .- plateau_height[i]] for i in 1:num_plateaus]
	plateau_potential = rand(num_plateaus) .* (po.plateau_potential[2] - po.plateau_potential[1]) .+ po.plateau_potential[1]

	function V(x, y)
		s = 0.0
		for i in 1:num_sin
			s += sin((x * directions[i, 1] + y * directions[i, 2]) * 2pi / po.L * sin_periods[i] + phase_shift[i]) * sin_amps[i]
		end
		
		if po.wall == :box
			if x < wall_thickness || x > po.L-wall_thickness || y < wall_thickness || y > po.L-wall_thickness
				s += Inf
			end
		end

		for i in 1:num_plateaus
			if x > plateau_corner[i][1] && x < plateau_corner[i][1] + plateau_width[i] && y > plateau_corner[i][2] && y < plateau_corner[i][2] + plateau_height[i]
				s += plateau_potential[i]
			end
		end
		return s
	end
end


# Vext_sin(x; n::Int, A::Number, φ::Number, L::Number) = A * sin(2π * x * n / L + φ)

# Vext_lin(x; x1::Number, x2::Number, E1::Number, E2::Number) = x > x1 && x < x2 ? E1 + (x - x1) * (E2 - E1) / (x2 - x1) : 0

# Vext_wall(x; xw::Number, L::Number) = x < xw || x > L - xw ? Inf : 0

# function generate_Vext(L::Number; num_sin=4, num_lin=rand(1:5), wall=true)
#     Avar = 1.0
#     sin_parameters = []
#     for n in 1:num_sin  # Generate random parameters for periodic sine functions with increasing frequency
#         push!(sin_parameters, (n = n, A = randn() * Avar, φ = rand() * 2π, L = L))
#     end
#     Evar = 1.0
#     lin_parameters = []
#     for _ in 1:num_lin  # Generate random parameters for discontinuous linear segments
#         push!(lin_parameters, (x1 = round(rand() * L, digits=2), x2 = round(rand() * L, digits=2), E1 = randn() * Evar, E2 = randn() * Evar))
#     end
#     xwmax = 1.0
#     wall_params = (xw = round(rand() * xwmax, digits=2), L = L)  # Set a random wall width
#     function (xy)  # Return a method which evaluates a combination of all functions with the chosen parameters above
#         result = 0.0
# 		  x = x
#         for sin_params in sin_parameters
#             result += Vext_sin(x; sin_params...)
#         end
#         for lin_params in lin_parameters
#             result += Vext_lin(x; lin_params...)
#         end
#         if wall
#             result += Vext_wall(x; wall_params...)
#         end
#         result
#     end
# end
