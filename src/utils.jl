using Plots
import Base: show
function show(io::IO, sys::GCMC_System) # FIXME does not work properly yet with periodic boundary conditions
	p = plot(legend=false, xlim=(-0, sys.L), ylim=(-0, sys.L), size=(500, 500))
	for pos in sys.positions
		scatter!(p, [pos[1]], [pos[2]], markercolor=:blue, markersize=sys.σ/sys.L*230)
		if pos[1] < sys.σ/2
			scatter!(p, [pos[1]+sys.L], [pos[2]], markercolor=:blue, markersize=sys.σ/sys.L*230)
		end
		if pos[1] > sys.L - sys.σ/2
			scatter!(p, [pos[1]-sys.L], [pos[2]], markercolor=:blue, markersize=sys.σ/sys.L*230)
		end
		if pos[2] < sys.σ/2
			scatter!(p, [pos[1]], [pos[2]+sys.L], markercolor=:blue, markersize=sys.σ/sys.L*230)
		end
		if pos[2] > sys.L - sys.σ/2
			scatter!(p, [pos[1]], [pos[2]-sys.L], markercolor=:blue, markersize=sys.σ/sys.L*230)
		end
	end
	display(p)
end
