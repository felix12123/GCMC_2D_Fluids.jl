using LinearAlgebra


# Energy of a certain particle i
function energy(sys::GCMC_System, i::Int)
	if i > sys.N
		@warn "Particle $i does not exist in the system."
		return 0
	end
	sys.Vext(sys.positions[i, :]...)
end
# Energy of a 
function energy(sys::GCMC_System, is::AbstractRange{Int}=1:sys.N)
	sum(sys.Vext.(sys.positions[is,1], sys.positions[is,2]))
end


# check if a particle can be inserted
function is_colliding_small(sys::GCMC_System, x::Real, y::Real; exclude_particle::Int=0)
	@assert exclude_particle <= sys.N
	for i in 1:sys.N
		if i == exclude_particle
			continue
		end
		if (mod(x - sys.positions[i,1] + sys.σ, sys.L) .- sys.σ)^2 + (mod(y - sys.positions[i,2] + sys.σ, sys.L) .- sys.σ)^2 < sys.σ^2
			return true
		end
	end
	return false
end




import Base.insert!
function insert!(sys::GCMC_System, x::Real, y::Real, check_col::Bool=true)
	if check_col && sys.is_colliding(sys, x, y)
		return false
	end
	sys.positions[sys.N+1, :] .= Float64(x), Float64(y)
	sys.N += 1
	return true
end

import Base.delete!
function delete!(sys::GCMC_System, i::Int)
	if i > sys.N
		return false
	end
	sys.positions[i,:] .= sys.positions[sys.N,:]
	sys.positions[sys.N,:] .= [NaN64, NaN64]
	sys.N -= 1
end
function move!(sys::GCMC_System, i::Int, x::Real, y::Real)
	if i > sys.N
		return false
	end
	sys.positions[i,1] = Float64(x)
	sys.positions[i,2] = Float64(y)
end

function try_insert!(sys::GCMC_System)
	x, y = rand(2) * sys.L # Random position
	if sys.is_colliding(sys, x, y) # Check if the position is already occupied
		return false
	end

	insert!(sys, x, y, false) # Insert the particle

	dE = energy(sys, sys.N) - sys.μ # Energy difference


	α_insert = clamp(sys.L^2 / (sys.N+1) * exp(-sys.β * dE), 0, 1) # Acceptance probability

	if rand() < α_insert # Accept the insertion
		return true
	end

	# If the insertion is rejected, delete the particle
	delete!(sys, sys.N)
	return false
end


function try_delete!(sys::GCMC_System)
	if sys.N == 0
		return false
	end

	i = rand(1:sys.N) # Random particle
	dE = -energy(sys, i) + sys.μ # Energy difference

	
	α_delete = clamp(sys.N / (sys.L^2) * exp(-sys.β * dE), 0, 1) # Acceptance probability

	if rand() < α_delete # Accept the deletion
		delete!(sys, i)
		return true
	end
	return false
end


function try_move!(sys::GCMC_System)
	if sys.N == 0
		return false
	end

	i = rand(1:sys.N) # Random particle

	x, y = sys.positions[i,:]
	dx = randn() * sys.mobility
	dy = randn() * sys.mobility
	new_x = mod(x+dx, sys.L) # Random move
	new_y = mod(y+dy, sys.L)

	# check if new position is colliding, if yes we can save time and not calculate the energy
	if sys.is_colliding(sys, new_x, new_y, exclude_particle=i)
		return false
	end

	E1 = energy(sys, i) # Energy before the move

	# Move the particle
	move!(sys, i, new_x, new_y)

	E2 = energy(sys, i) # Energy after the move
	dE = E2-E1 # Energy difference
	α_move = min(1, exp(-sys.β * dE)) # Acceptance probability

	if rand() < α_move # Accept the move
		return true
	end

	# If the move is rejected, move back to the old position
	move!(sys, i, x, y)

	return false # Move rejected
end


