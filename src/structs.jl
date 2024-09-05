# predefine function is_colliding_small
is_colliding_small() = error("predefined function used")

mutable struct GCMC_TrainingData
	L::Float64 # Box length
	σ::Float64 # Particle diameter
	μ_range::Tuple{Float64, Float64} # Chemical potential range
	β::Float64 # Inverse temperature
	mobility::Float64 # Particle mobility (standard deviation of random displacement)

	dx::Float64 # Discretization step for histograms

	steps::Int64 # Number of steps
	therm_steps::Int64 # Number of thermalization steps
	sample_interval::Int64 # Interval to sample observables

	threads::Tuple{Int64, Int64} # Number of threads for outer and inner loop (repetitions)
	repetitions::Int64 # Number of repetitions per potential
	num_systems::Int # Number of systems to generate

	# we first roll the move probability. if it fails, we roll the insert probability.
	# if both fail, we try to delete.
	move_prob::Float64 # Probability of moving a particle
	insert_prob::Float64 # Probability of inserting a particle

	folder::String # Folder to save data

	accept_condition::Function # Function to check if a density distribution is accepted
	Vext_generator::Function # Function to generate a new external potential
	rho_smooth_func::Function # Function to smooth the density distribution

	function GCMC_TrainingData(folder::String;
			L::Real=10.0,
			σ::Real=1,
			μ_range::Tuple{Real, Real}=(-1,10),
			β::Real=1,
			mobility::Real=0.15σ,
			dx::Real=0.05,
			steps::Real=10^6,
			therm_steps::Real=3*10^4,
			sample_interval::Real=500,
			threads::Union{<:Real, Tuple{<:Real, <:Real}}=Threads.nthreads(),
			repetitions::Real=1,
			num_systems::Real=20,
			move_prob::Real=0.9,
			insert_prob::Real=0.5,
			accept_condition::Function=(x...)->true,
			Vext_generator::Function=(x...)->0.0,
			rho_smooth_func::Function=(x, y...)->x
		)
		if isa(threads, Real)
			half_threads = floor(Int, threads/2)
			threads = (half_threads, threads - half_threads)
		end
		new(L, σ, μ_range, β, mobility, dx, ceil(Int, steps), ceil(Int, therm_steps), ceil(Int, sample_interval), threads, ceil(Int, repetitions), ceil(Int, num_systems), move_prob, insert_prob, folder, accept_condition, Vext_generator, rho_smooth_func)
	end
end


mutable struct GCMC_System
	L::Float64 # Box length
	N::Int64 # Number of particles
	σ::Float64 # Particle diameter
	μ::Float64 # Chemical potential
	β::Float64 # Inverse temperature
	Vext::Function # External potential
	positions::Matrix{Float64}
	is_colliding::Function # Function to check if a particle is colliding. might be switched depending on amount of particles
	mobility::Float64 # Particle mobility (standard deviation of random displacement)

	# we first roll the move probability. if it fails, we roll the insert probability.
	# if both fail, we try to delete.
	move_prob::Float64 # Probability of moving a particle
	insert_prob::Float64 # Probability of inserting a particle

	dx::Float64 # Discretization step for histogram
	particle_shape::Symbol # Shape of the particle
	function GCMC_System(;L::Real=10, σ::Real=1.0, μ::Real=0.0, β::Real=1.0, Vext::Function=(x) -> 0.0, mobility=0.1σ, move_prob::Float64=0.8, insert_prob::Float64=0.5, dx::Float64=0.05, particle_shape::Symbol=:circle)
		max_N = ceil(Int, L^2/(σ^2 * (particle_shape == :circle ? π/4 : 1)))  * 1000 #FIXME remove the 1000
		positions = zeros(Float64, max_N, 2)
		positions .= NaN64
		new(L, 0, σ, μ, β, Vext, positions, is_colliding_small, mobility, move_prob, insert_prob, dx, particle_shape)
	end
end

# we need to keep track of the average density
mutable struct Histogram
	dx::Float64
	ρ::Matrix{Float64}
	count::Int
	function Histogram(sys::GCMC_System)
		hist_len = floor(Int, sys.L/sys.dx)
		rho = zeros(hist_len, hist_len)
		new(sys.dx, rho, 0)
	end
end


mutable struct g_Histogram
	dx::Float64
	ρ::Vector{Int}
	count::Int
	function g_Histogram(sys::GCMC_System)
		new(sys.dx, zeros(Int, floor(Int, sys.L/2/sys.dx)), 0)
	end
end