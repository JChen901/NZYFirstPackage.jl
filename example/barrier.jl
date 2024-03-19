using Makie: RGBA # for visualization
using Makie, GLMakie
using MyFirstPackage # our package

# Set up the visualization with Makie:
lb = example_d3q19()
vorticity = Observable(Vec3f.(reshape(curl(momentum.(Ref(lb.config), lb.grid)),:)))
xlim, ylim, zlim = size(lb.grid)

pts = reshape([Point3f(i*4/xlim, j*4/ylim, k*4/zlim) for i in 1:xlim, j in 1:ylim, k in 1:zlim], :)
fig, ax, plot = arrows(pts, vorticity, fxaa=true, # turn on anti-aliasing
    linecolor = :grey, arrowcolor = :blue,
    linewidth = 0.1, arrowsize = Vec3f(0.1, 0.1, 0.1),
    align = :center, axis=(type=Axis3,)
)

xlims!(ax, 0, 4)
ylims!(ax, 0, 4)
zlims!(ax, 0, 4)
# Recording the simulation
record(fig, joinpath(@__DIR__, "barrier.mp4"), 1:100; framerate = 10) do i
    for i=1:20
        step!(lb)
    end
    vorticity[] = Vec3f.(reshape(curl(momentum.(Ref(lb.config), lb.grid)),:))
end

using BenchmarkTools
@benchmark step!($(deepcopy(lb)))


