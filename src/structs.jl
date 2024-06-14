
mutable struct GCMC_Simulation
	L::Float64 # Box length
	σ::Float64 # Particle diameter
	μ_range::Tuple{Float64, Float64} # Chemical potential range
	β::Float64 # Inverse temperature
	Vext::Function # External potential
	mobility::Float64 # Particle mobility (standard deviation of random displacement)

	# we first roll the move probability. if it fails, we roll the insert probability.
	# if both fail, we try to delete.
	move_prob::Float64 # Probability of moving a particle
	insert_prob::Float64 # Probability of inserting a particle
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

	# we first roll the move probability. if it fails, we roll the insert probability.
	# if both fail, we try to delete.
	move_prob::Float64 # Probability of moving a particle
	insert_prob::Float64 # Probability of inserting a particle

	function GCMC_System(;L::Real=10, σ::Real=1.0, μ::Real=0.0, β::Real=1.0, Vext::Function=(x) -> 0.0, mobility=0.25σ, move_prob::Float64=0.8, insert_prob::Float64=0.5)
		new(L, 0, σ, μ, β, Vext, Float64[], is_colliding_small, mobility, move_prob, insert_prob)
	end
end

# we need to keep track of the average density
mutable struct Histogram
	dx::Float64
	ρ::Matrix{Float64}
	count::Int
	function Histograms(sys::GCMC_Simulation; dx=0.05)
		new(floor(sys.L/dx)^2, dx, zeros(floor(sys.L/dx), floor(sys.L/dx)), 0)
	end
end