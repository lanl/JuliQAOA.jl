
@testset "Cost Functions" begin
    g = Graphs.SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 1, 4)
    add_edge!(g, 2, 3)
    add_edge!(g, 2, 5)
    add_edge!(g, 3, 5)
    add_edge!(g, 4, 5)
    
    weights = Dict((1,2)=>4,(1,4)=>2,(2,3)=>1,(2,5)=>3,(3,5)=>3,(4,5)=>1)
    
    @testset "maxcut" begin
        @test maxcut(g, [0, 0, 1, 1, 1]) == 3
        @test maxcut(g, [0, 1, 0, 0, 1]; weights=weights) == 9
    end
    
    @testset "bisection" begin
        @test bisection(g, [0, 1, 0, 1, 0]) == 5
    end
    
    @testset "densest k-subgraph" begin
        @test densest_ksubgraph(g, [0, 1, 1, 0, 1]) == 3
        @test densest_ksubgraph(g, [1, 1, 0, 1, 0]) == 2
    end
    
    @testset "k-vertex cover" begin
        @test kvertex_cover(g, [1, 0, 1, 0, 1]) == 6
        @test kvertex_cover(g, [0, 1, 1, 0, 0]) == 4
    end
    
    @testset "kSAT" begin
        instance = [[2, -1, -4], [-4, 5, -2], [6, 5, 4], [-3, -4, 2], [3, -6, -5]]
        @test kSAT(instance, [0, 1, 1, 0, 1, 1]) == 5
        @test kSAT(instance, [1, 0, 1, 0, 0, 0]) == 4
    end
    
    @testset "spin energy" begin
        model = Dict((1,)=>-1.2, (2,)=>2.7, (3,)=>1.4, (1,2)=>2.2, (1,3)=>-5.9, (2,3)=>-2.1, 
        (1,2,3)=>-4.3)
        @test spin_energy(model, [0, 1, 1]) ≈ -8
        @test spin_energy(model, [1, 0, 1]) ≈ -7.8
    end
end
