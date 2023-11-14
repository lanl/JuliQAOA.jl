
@testset "Utilities" begin

    @testset "states" begin
        @test collect(states(1)) == [[0], [1]]
        @test collect(states(2)) == [[0, 0], [1, 0], [0, 1], [1, 1]]
        @test collect(states(3)) == [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    end

    @testset "gospers_hack" begin
        @test JuliQAOA.gospers_hack(6) == 9
        @test JuliQAOA.gospers_hack(3) == 5
        @test JuliQAOA.gospers_hack(3,5) == [3,5,6,9,10]
    end

    @testset "iterate" begin
        iter = JuliQAOA.DickeIterator(4, 2, 3, 16)
        @test iterate(iter) == ([1, 1, 0, 0], 5)
        @test iterate(iter, 5) == ([1, 0, 1, 0], 6)
    end

    @testset "length" begin 
        @test length(JuliQAOA.DickeIterator(5, 2, 3, 16)) == 10
        @test length(JuliQAOA.DickeIterator(5, 3, 3, 16)) == 10
    end

    @testset "dicke_states" begin
        @test collect(dicke_states(4,0)) == [[0, 0, 0, 0]]
        @test sort(collect(dicke_states(4,1))) == sort([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        @test sort(collect(dicke_states(4,2))) == sort([[1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1], [0, 1, 1, 0], [0, 1, 0, 1], [0, 0, 1, 1]])
        @test sort(collect(dicke_states(4,3))) == sort([[1, 1, 1, 0], [1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 1, 1]])
        @test collect(dicke_states(4,4)) == [[1, 1, 1, 1]]
    end

end