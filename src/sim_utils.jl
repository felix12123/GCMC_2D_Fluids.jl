using LinearAlgebra


# Energy of a certain particle i
function energy(sys::GCMC_System, i::Int)
	sys.Vext(sys.positions[i])
end
# Energy of a 
function energy(sys::GCMC_System, is::AbstractRange{Int})
	sum(sys.Vext.(sys.positions[is]))
end


# check if a particle can be inserted
function is_colliding_small(sys::GCMC_System, pos::Tuple{<:Real, <:Real})
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
function insert!(sys::GCMC_System, pos::Tuple{<:Real, <:Real}, check_col::Bool=true)
	if check_col && sys.is_colliding(sys, pos)
		return false
	end
	push!(sys.positions, Float64.(pos))
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
	α_insert = min(1, 1 - (sys.N * sys.σ^2 + 1) / sys.L^2 * exp(-sys.β * dE)) # Acceptance probability
	
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
	# FIXME this might be wrong:
	α_delete = min(1, sys.σ^2 * sys.N / (2*sys.L^2) * exp(-sys.β * dE)) # Acceptance probability

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

	pos = sys.positions[i]
	new_pos = mod.(pos .+ randn(2) * sys.mobility, sys.L) |> Tuple # Random move

	E1 = energy(sys, i) # Energy before the move

	delete!(sys, i)

	if insert!(sys, new_pos)
		E2 = energy(sys, sys.N) # Energy after the move
		dE = E2-E1 # Energy difference

		# FIXME this might be wrong:
		α_move = min(1, exp(-sys.β * dE)) # Acceptance probability

		if rand() < α_move # Accept the move
			return true
		end
		delete!(sys, sys.N)
	end

	insert!(sys, pos, false) # Reinsert the particle
	return false # Move rejected
end


