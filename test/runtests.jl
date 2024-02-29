using MyFirstPackage
using Test

@testset "flip_direction_index" begin
    d=D2Q9()
    @test flip_direction_index(d, 1) == 9
end
