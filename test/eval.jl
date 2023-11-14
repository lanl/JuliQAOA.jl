
@testset "Statevector Simulation" begin

    @testset "3-round transverse field maxcut" begin
        n = 6

        g = erdos_renyi(n, 0.5)
        obj_vals = [maxcut(g, x) for x in states(n)]
    
        mixer = mixer_x(n)
        
        p = 3
        angles = rand(2*p)

        # test statevector
        svtest = statevector(angles, mixer, obj_vals)

        init_state = ones(Complex, 2^n)/sqrt(2^n)
        mixing_mat = sum_symmetric_krons([PX,I,I,I,I,I])
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

        @test svtest ≈ svtarget

        # change initial state
        init_state = random_sv(2^n)

        svtest = statevector(init_state, angles, mixer, obj_vals)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

        @test svtest ≈ svtarget

        # check probabilities

        probtest = probabilities(init_state, angles, mixer, obj_vals)
        probtarget = abs2.(svtarget)

        @test probtest ≈ probtarget

        # check expectation value

        expvaltest = exp_value(init_state, angles, mixer, obj_vals)
        expvaltarget = sum(abs2.(svtarget) .* obj_vals)

        @test expvaltest ≈ expvaltarget

        # check expectation value with different measurement

        measure = rand(2^n)
        measuretest = exp_value(init_state, angles, mixer, obj_vals, measure)
        measuretarget = sum(abs2.(svtarget) .* measure)

        @test measuretest ≈ measuretarget

        # check expectation value with different measurement but traditional init state

        measure = rand(2^n)
        measuretest = exp_value(angles, mixer, obj_vals, measure)

        init_state = ones(Complex, 2^n)/sqrt(2^n)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)
        measuretarget = sum(abs2.(svtarget) .* measure)

        @test measuretest ≈ measuretarget

    end

    @testset "3-round multi-X mixer with Sherrington-Kirkpatrick" begin
        n = 6

        model = sk_model(n)
        obj_vals = [spin_energy(model, x) for x in states(n)]
    
        mixer = mixer_x(n, [2]) + .3*mixer_x(n,[3])
        mixing_mat = sum_symmetric_krons([PX,PX,I,I,I,I]) + .3*sum_symmetric_krons([PX,PX,PX,I,I,I])

        p = 3
        angles = rand(2*p)

        # test statevector
        svtest = statevector(angles, mixer, obj_vals)

        init_state = ones(Complex, 2^n)/sqrt(2^n)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

        @test svtest ≈ svtarget

        # change initial state
        init_state = random_sv(2^n)

        svtest = statevector(init_state, angles, mixer, obj_vals)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

        @test svtest ≈ svtarget

        # check probabilities

        probtest = probabilities(init_state, angles, mixer, obj_vals)
        probtarget = abs2.(svtarget)

        @test probtest ≈ probtarget

        # check expectation value

        expvaltest = exp_value(init_state, angles, mixer, obj_vals)
        expvaltarget = sum(abs2.(svtarget) .* obj_vals)

        @test expvaltest ≈ expvaltarget

        # check expectation value with different measurement

        measure = rand(2^n)
        measuretest = exp_value(init_state, angles, mixer, obj_vals, measure)
        measuretarget = sum(abs2.(svtarget) .* measure)

        @test measuretest ≈ measuretarget

        # check expectation value with different measurement but traditional init state

        measure = rand(2^n)
        measuretest = exp_value(angles, mixer, obj_vals, measure)

        init_state = ones(Complex, 2^n)/sqrt(2^n)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)
        measuretarget = sum(abs2.(svtarget) .* measure)

        @test measuretest ≈ measuretarget

    end

    @testset "3-round clique densest k-subgraph" begin
        n = 6
        k = 3

        g = erdos_renyi(n, 0.5)
        obj_vals = [densest_ksubgraph(g, x) for x in dicke_states(n,k)]
    
        mixer = mixer_clique(n,k)
        
        p = 3
        angles = rand(2*p)

        mixing_mat_full = 1/2*(sum_symmetric_krons([PX,PX,I,I,I,I])+sum_symmetric_krons([PY,PY,I,I,I,I]))
        target_locs = [8, 12, 14, 15, 20, 22, 23, 26, 27, 29, 36, 38, 39, 42, 43, 45, 50, 51, 53, 57]
        mixing_mat = mixing_mat_full[target_locs, target_locs]

        N = binomial(n,k)

        # test statevector
        svtest = statevector(angles, mixer, obj_vals)

        init_state = ones(Complex, N)/sqrt(N)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

        @test svtest ≈ svtarget

        # change initial state
        init_state = random_sv(N)

        svtest = statevector(init_state, angles, mixer, obj_vals)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

        @test svtest ≈ svtarget

        # check probabilities

        probtest = probabilities(init_state, angles, mixer, obj_vals)
        probtarget = abs2.(svtarget)

        @test probtest ≈ probtarget

        # check expectation value

        expvaltest = exp_value(init_state, angles, mixer, obj_vals)
        expvaltarget = sum(abs2.(svtarget) .* obj_vals)

        @test expvaltest ≈ expvaltarget

        # check expectation value with different measurement

        measure = rand(N)
        measuretest = exp_value(init_state, angles, mixer, obj_vals, measure)
        measuretarget = sum(abs2.(svtarget) .* measure)

        @test measuretest ≈ measuretarget

        # check expectation value with different measurement but traditional init state

        measure = rand(N)
        measuretest = exp_value(angles, mixer, obj_vals, measure)

        init_state = ones(Complex, N)/sqrt(N)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)
        measuretarget = sum(abs2.(svtarget) .* measure)

        @test measuretest ≈ measuretarget

    end

    @testset "3-round ring k-vertex cover" begin
        n = 6
        k = 3

        g = erdos_renyi(n, 0.5)
        obj_vals = [kvertex_cover(g, x) for x in dicke_states(n,k)]
    
        mixer = mixer_ring(n,k)
        
        p = 3
        angles = rand(2*p)

        mixing_matx = 1/2*(kron(PX,PX,I,I,I,I)+kron(I,PX,PX,I,I,I)+kron(I,I,PX,PX,I,I)+kron(I,I,I,PX,PX,I)+kron(I,I,I,I,PX,PX))
        mixing_maty = 1/2*(kron(PY,PY,I,I,I,I)+kron(I,PY,PY,I,I,I)+kron(I,I,PY,PY,I,I)+kron(I,I,I,PY,PY,I)+kron(I,I,I,I,PY,PY))
        mixing_mat_full = mixing_matx + mixing_maty
        target_locs = [8, 12, 14, 15, 20, 22, 23, 26, 27, 29, 36, 38, 39, 42, 43, 45, 50, 51, 53, 57]
        mixing_mat = mixing_mat_full[target_locs, target_locs]

        N = binomial(n,k)

        # test statevector
        svtest = statevector(angles, mixer, obj_vals)

        init_state = ones(Complex, N)/sqrt(N)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

        @test svtest ≈ svtarget

        # change initial state
        init_state = random_sv(N)

        svtest = statevector(init_state, angles, mixer, obj_vals)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)

        @test svtest ≈ svtarget

        # check probabilities

        probtest = probabilities(init_state, angles, mixer, obj_vals)
        probtarget = abs2.(svtarget)

        @test probtest ≈ probtarget

        # check expectation value

        expvaltest = exp_value(init_state, angles, mixer, obj_vals)
        expvaltarget = sum(abs2.(svtarget) .* obj_vals)

        @test expvaltest ≈ expvaltarget

        # check expectation value with different measurement

        measure = rand(N)
        measuretest = exp_value(init_state, angles, mixer, obj_vals, measure)
        measuretarget = sum(abs2.(svtarget) .* measure)

        @test measuretest ≈ measuretarget

        # check expectation value with different measurement but traditional init state

        measure = rand(N)
        measuretest = exp_value(angles, mixer, obj_vals, measure)

        init_state = ones(Complex, N)/sqrt(N)
        svtarget = qaoa_simulate_bruteforce(init_state, angles, mixing_mat, obj_vals)
        measuretarget = sum(abs2.(svtarget) .* measure)

        @test measuretest ≈ measuretarget

    end
end