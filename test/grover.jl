
@testset "Grover" begin

    @testset "grover mixer -- unconstrained" begin
        n = 6

        g = erdos_renyi(n, 0.5)
        obj_vals = [maxcut(g, x) for x in states(n)]
    
        m1 = mixer_x(n,0:n)/2^n
        m2 = mixer_grover(n)
        
        p = 3
        angles = rand(2*p)

        # test statevector
        @test statevector(angles, m1, obj_vals) ≈ statevector(angles, m2, obj_vals)

        # test probabilities
        @test probabilities(angles, m1, obj_vals) ≈ probabilities(angles, m2, obj_vals)

        # test expectation value
        @test exp_value(angles, m1, obj_vals) ≈ exp_value(angles, m2, obj_vals)

        # test expectation value with different measurement
        measure = rand(2^n)
        @test exp_value(angles, m1, obj_vals, measure) ≈ exp_value(angles, m2, obj_vals, measure)

        # Note: don't need to test with different initial state, as Grover mixers work from 
        # the uniform superposition by default. Can do the x-mixer version if you want to 
        # try a non-standard initial state
        init_state = random_sv(2^n)
        @test !(statevector(init_state, angles, m1, obj_vals) ≈ statevector(init_state, angles, m2, obj_vals))
        @test statevector(init_state, angles, m2, obj_vals) ≈ statevector(angles, m2, obj_vals)
        @test exp_value(init_state, angles, m2, obj_vals) ≈ exp_value(angles, m2, obj_vals)
    end


    @testset "grover mixer -- consrained" begin
        n = 6
        k = 3
        N = binomial(n,k)

        g = erdos_renyi(n, 0.5)
        obj_vals = [densest_ksubgraph(g, x) for x in dicke_states(n,k)]
    
        mixer = mixer_grover(n,k)
        mixing_mat = ones(N,N)/N # = |F><F| where |F> = 1/sqrt(N) * |1,1,...,1>
        
        p = 3
        angles = rand(2*p)

       # test statevector
       svtest = statevector(angles, mixer, obj_vals)

       init_state = ones(Complex, N)/sqrt(N)
       svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

       @test svtest ≈ svtarget

       # test probabilities
       probtest = probabilities(angles, mixer, obj_vals)
       probtarget = abs2.(svtarget)

       @test probtest ≈ probtarget

       # test expectation value
       expvaltest = exp_value(angles, mixer, obj_vals)
       expvaltarget = sum(abs2.(svtarget) .* obj_vals)

       @test expvaltest ≈ expvaltarget

       # test expectation value with different measurement
       measure = rand(N)
       measuretest = exp_value(angles, mixer, obj_vals, measure)
       measuretarget = sum(abs2.(svtarget) .* measure)

       @test measuretest ≈ measuretarget
    end
end