using LinearAlgebra

mutable struct GCMC_Simulation
	L::Float64 # Box length
	σ::Float64 # Particle diameter
	μ_range::Tuple{Float64, Float64} # Chemical potential range
	β::Float64 # Inverse temperature
	Vext::Function # External potential
	mobility::Float64 # Particle mobility (standard deviation of random displacement)
end

mutable struct GCMC_System
	L::Float64 # Box length
	N::Int64 # Number of particles
	σ::Float64 # Particle diameter
	μ::Float64 # Chemical potential
	β::Float64 # Inverse temperature
	Vext::Function # External potential
	positions::Vector{Tuple{Float64, Float64}}
	is_colliding::Function # Function to check if a particle is colliding. might be switched depending on amount of particles
	mobility::Float64 # Particle mobility (standard deviation of random displacement)


	function GCMC_System(;L::Real=10, σ::Real=1.0, μ::Real=0.0, β::Real=1.0, Vext::Function=(x, y) -> 0.0)
		new(L, 0, σ, μ, β, Vext, Float64[], is_colliding_small)
	end
end


# Energy of a certain particle i
function energy(sys::GCMC_System, i::Int)
	Vext(sys.positions[i])
end
# Energy of a 
function energy(sys::GCMC_System, is::AbstractRange{Int})
	sum(Vext.(sys.positions[is]))
end


# check if a particle can be inserted
function is_colliding_small(sys::GCMC_System, pos::Tuple{Float64, Float64})
	if isempty(sys.positions)
		return false
	end

	if any([norm(pos .- pos_i) < sys.σ for pos_i in sys.positions])
		return true
	end

	# if particle is partly outside the box, check periodic boundary conditions
	if pos[1] < sys.σ/2
		new_pos = pos .+ [sys.L, 0]
		if any([norm(new_pos .- pos_i) < sys.σ for pos_i in sys.positions])
			return true
		end
	end
	if pos[1] > sys.L - sys.σ/2
		new_pos = pos .- [sys.L, 0]
		if any([norm(new_pos .- pos_i) < sys.σ for pos_i in sys.positions])
			return true
		end
	end
	if pos[2] < sys.σ/2
		new_pos = pos .+ [0, sys.L]
		if any([norm(new_pos .- pos_i) < sys.σ for pos_i in sys.positions])
			return true
		end
	end
	if pos[2] > sys.L - sys.σ/2
		new_pos = pos .- [0, sys.L]
		if any([norm(new_pos .- pos_i) < sys.σ for pos_i in sys.positions])
			return true
		end
	end
	return false
end

import Base.insert!
function insert!(sys::GCMC_System, pos::Tuple{Float64, Float64}, check_col::Bool=true)
	if check_col && sys.is_colliding(sys, pos)
		return false
	end
	push!(sys.positions, pos)
	sys.N += 1
	return true
end

import Base.delete!
function delete!(sys::GCMC_System, i::Int)
	deleteat!(sys.positions, i)
	sys.N -= 1
end

function try_insert!(sys::GCMC_System)
	x, y = rand(2) * sys.L # Random position
	if sys.is_colliding(sys, (x, y)) # Check if the position is already occupied
		return false
	end

	insert!(sys, (x, y), false) # Insert the particle

	dE = energy(sys, sys.N) - sys.μ # Energy difference
	# FIXME this might be wrong:
	α_insert = min(1, 1 - 2 * sys.L^2 / (sys.L * sys.σ^2 + 1) * exp(-sys.β * dE)) # Acceptance probability
	
	if rand() < α_insert # Accept the insertion
		return true
	end

	# If the insertion is rejected, delete the particle
	delete!(sys, sys.N)
	return false
end


function try_delete!(sys::GCMC_System)
	i = rand(1:sys.N) # Random particle
	dE = -energy(sys, i) + sys.μ # Energy difference
	# FIXME this might be wrong:
	α_delete = min(1, sys.σ^2 * sys.N / (2*sys.L^2) * exp(-sys.β * dE)) # Acceptance probability

	if rand() < α_delete # Accept the deletion
		delete!(sys, i)
		return true
	end
	return false
end


function try_move!(sys::GCMC_System)
	i = rand(1:sys.N) # Random particle

	x, y = mod.(sys.positions[i] + randn(2) * sys.mobility, sys.L) # Random move

	if sys.is_colliding(sys, (x, y)) # Check if the position is already occupied
		return false
	end

	#...
	
end