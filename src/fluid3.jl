"""
	D3Q9 <: AbstractLBConfig{3, 9}

A lattice Boltzmann configuration for 3D, 19-velocity model.
"""

# println("Hello from fluid3.jl")
# pts = Point[]
# rev_pts = Point[]
# for (i, j, k) in Base.product(repeat([[-1, 0, 1]], 3)...)
# 	if !(Point(i, j, k) in rev_pts) && (abs(i) + abs(j) + abs(k)) != 3
# 		push!(pts, Point(i, j, k))
# 		push!(rev_pts, Point(-i, -j, -k))
# 	end
# end

# length(pts)
# length(rev_pts)

# for pt in pts
# 	println("Point(", pt[1], ",", pt[2], ",", pt[3], "),")
# end

# for pt in reverse(rev_pts)
# 	println("Point(", pt[1], ",", pt[2], ",", pt[3], "),")
# end

struct D3Q19 <: AbstractLBConfig{3, 19} end
function directions(::D3Q19)
	return (
		Point(0, -1, -1),
		Point(-1, 0, -1),
		Point(0, 0, -1),
		Point(1, 0, -1),
		Point(0, 1, -1),
		Point(-1, -1, 0),
		Point(0, -1, 0),
		Point(1, -1, 0),
		Point(-1, 0, 0),
		Point(0, 0, 0),
		Point(1, 0, 0),
		Point(-1, 1, 0),
		Point(0, 1, 0),
		Point(1, 1, 0),
		Point(0, -1, 1),
		Point(-1, 0, 1),
		Point(0, 0, 1),
		Point(1, 0, 1),
		Point(0, 1, 1),
	)
end

# directions[k] is the opposite of directions[flip_direction_index(k)
function flip_direction_index(::D3Q19, i::Int)
	return 20 - i
end

# streaming step
function stream!(
	lb::AbstractLBConfig{3, N},  # lattice configuration
	newgrid::AbstractArray{D,3}, # the updated grid
	grid::AbstractArray{D,3}, # the original grid
	barrier::AbstractArray{Bool,3}, # the barrier configuration
) where {N, T, D <: Cell{N, T}}
	ds = directions(lb)
	@inbounds for ci in CartesianIndices(newgrid)
		i, j, l = ci.I
		newgrid[ci] = Cell(
			ntuple(N) do k # collect the densities
				ei = ds[k]
				m, n, q = size(grid)
				i2, j2, l2 = mod1(i - ei[1], m), mod1(j - ei[2], n), mod1(k - ei[3], q)
				if barrier[i2, j2, l2]
					# if the cell is a barrier, the fluid flows back
					density(grid[i, j, l], flip_direction_index(lb, k))
				else
					# otherwise, the fluid flows to the neighboring cell
					density(grid[i2, j2, l2], k)
				end
			end,
		)
	end
end

# https://arxiv.org/pdf/1511.04633.pdf
# the distribution of the 19 velocities at the equilibrium state
weights(::D3Q19) = (1 / 36,
	1 / 36,
	1 / 18, 1 / 36,
	1 / 36,
	1 / 36,
	1 / 18,
	1 / 36, 1 / 18,
	1 / 3,
	1 / 18,
	1 / 36,
	1 / 18,
	1 / 36,
	1 / 36,
	1 / 36,
	1 / 18,
	1 / 36,
	1 / 36)

# https://www.egr.msu.edu/~kutay/LBsite/
function _equilibrium_density(u, ei)
	# the equilibrium density of the fluid with a specific mean momentum
	return (1 + 3 * dot(ei, u) + 9 / 2 * dot(ei, u)^2 - 3 / 2 * dot(u, u))
end

"""
	curl(u::AbstractMatrix{Point3D{T}})

Compute the curl of the momentum field in 3D, which is defined as:
```math
∂u_y/∂x−∂u_x/∂y,
∂u_y/∂x−∂u_x/∂y,
∂u_y/∂x−∂u_x/∂y,

```
"""
function curl(u::Array{Point3D{T},3}) where T
	return map(CartesianIndices(u)) do ci
		i, j, k = ci.I
		m, n, p = size(u)
		uy = u[mod1(i + 1, m), j, k][2] - u[mod1(i - 1, m), j, k][2]
		ux = u[i, mod1(j + 1, n), k][1] - u[i, mod1(j - 1, n), k][1]
		uz = u[i, j, mod1(k + 1, p)][3] - u[i, j, mod1(k - 1, p)][3]

		return (uy - ux, ux - uz, uz - uy)
	end
end

function example_d3q19(;
	height = 10, width = 30, depth = 10,
	u0 = Point(0.0, 0.1, 0.2)) # initial and in-flow speed
	# Initialize all the arrays to steady rightward flow:
	rho = equilibrium_density(D3Q19(), 1.0, u0)
	rgrid = fill(rho, height, width,depth)

	# Initialize barriers:
	barrier = falses(height, width, depth)  # True wherever there's a barrier
	mid = div(height, 2)
	barrier[(mid - 1):(mid + 1), div(height, 2), div(depth, 2)] .= true  # simple linear barrier

	return LatticeBoltzmann(D3Q19(), rgrid, barrier)
end

