function finite_diff(f, x0; eps=1e-10)
    new_g = similar(x0)
    for i in eachindex(x0)
        eps_vec = zeros(length(x0))
        eps_vec[i] += eps
        new_g[i] = (f(x0 .+ eps_vec) - f(x0))/eps
    end
    return new_g
end

@testset "Angle Finding" begin

    @testset "Ancillary functions" begin

        @testset "get_operator_period" begin
            @test JuliQAOA.get_operator_period([2, 4]) ≈ π
            @test JuliQAOA.get_operator_period([2.000001, 4.0]) == Inf
            @test JuliQAOA.get_operator_period([2.000001, 4.0]; tol=1e-5) ≈ π
            @test JuliQAOA.get_operator_period([2.5, 4.0]) == Inf
        end

        @testset "clean_angles" begin
            mixer = mixer_x(6)
            @test JuliQAOA.clean_angles([1.5π, 2.5π, 3.5π, 4.5π], mixer, [1.0, 4.0]) == [0.5π, 0.5π, 1.5π, 0.5π]
            mixer = mixer_ring(6,3)
            @test JuliQAOA.clean_angles([1.5π, 2.5π, 3.5π, 4.5π], mixer, [2.5, 4.0]) == [1.5π, 2.5π, 3.5π, 4.5π]
        end

        @testset "guess_angles" begin
            @test JuliQAOA.guess_angles([[1, 2], [1, 2, 3, 4]]) == [1, 2, 2, 3, 4, 4]
        end

    end

    @testset "Gradient" begin
        n = 6
        p = 5
        x0 = rand(2*p)
        
        g = erdos_renyi(n, 0.5)
        
        @testset "unconstrained" begin
            obj_vals = [maxcut(g, x) for x in states(n)]
            measure = rand(2^n)
            init_state = random_sv(2^n)
            ref_state = deepcopy(init_state)
            mixer = mixer_x(n)
            
            # default init state, measure
            fd_res = finite_diff(x->exp_value(x, mixer, obj_vals), x0)
            ad_res = grad(x0, mixer, obj_vals)
            
            @test all(abs.(fd_res - ad_res) .< 1e-4)
            
            # different init state, default measure
            fd_res = finite_diff(x->exp_value(init_state, x, mixer, obj_vals), x0)
            ad_res = grad(init_state, x0, mixer, obj_vals)
            
            @test all(abs.(fd_res - ad_res) .< 1e-4)
            @test ref_state == init_state

            # default init state, different measure
            fd_res = finite_diff(x->exp_value(x, mixer, obj_vals, measure), x0)
            ad_res = grad(x0, mixer, obj_vals, measure)
            
            @test all(abs.(fd_res - ad_res) .< 1e-4)
            
            # different init state + different measure
            fd_res = finite_diff(x->exp_value(init_state, x, mixer, obj_vals, measure), x0)
            ad_res = grad(init_state, x0, mixer, obj_vals, measure)
            
            @test all(abs.(fd_res - ad_res) .< 1e-4)
            @test init_state == ref_state
        end

        @testset "constrained" begin
            k = 3
            obj_vals = [densest_ksubgraph(g, x) for x in dicke_states(n,k)] 
            mixer = mixer_clique(n, k)
            
            fd_res = finite_diff(x->exp_value(x, mixer, obj_vals), x0)
            ad_res = grad(x0, mixer, obj_vals)
            
            @test all(abs.(fd_res - ad_res) .< 1e-4)
        end

        # NOTE: grad() doesn't work with the fast Grover mixer implementation
        # instead, you have to create a general mixer of all ones(N,N)/N 
        #@testset "grover" begin
            #k = 3
            #obj_vals = [densest_ksubgraph(g, x) for x in dicke_states(n,k)] 
            #mixer = mixer_grover(n, k)
            
            #fd_res = finite_diff(x->exp_value(x, mixer, obj_vals), x0)
            #ad_res = grad(x0, mixer, obj_vals)
            
            #@test all(abs.(fd_res - ad_res) .< 1e-4)
        #end
    end
end