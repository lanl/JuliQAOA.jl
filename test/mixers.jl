
@testset "Mixers" begin

    @testset "x mixer" begin
        mixer = mixer_x(4)
        full_mixer = kron(H,H,H,H)*Diagonal(mixer.d)*kron(H,H,H,H)
        target = sum_symmetric_krons([PX,I,I,I])
        @test full_mixer ≈ target

        mixer = mixer_x(4, [2,3])
        full_mixer = kron(H,H,H,H)*Diagonal(mixer.d)*kron(H,H,H,H)
        target = sum_symmetric_krons([PX,PX,I,I]) + sum_symmetric_krons([PX,PX,PX,I])
        @test full_mixer ≈ target

        mixer = mixer_x(5, 0:5)
        full_mixer = kron(H,H,H,H,H)*Diagonal(mixer.d)*kron(H,H,H,H,H)
        target = ones(2^5,2^5)
        @test full_mixer ≈ target
    end

    @testset "clique mixer" begin
        mixer = mixer_clique(4,2)
        full_mixer = mixer.v*Diagonal(mixer.d)*mixer.vinv
        target = 1/2*(sum_symmetric_krons([PX,PX,I,I])+sum_symmetric_krons([PY,PY,I,I]))
        locs = [4, 6, 7, 10, 11, 13]
        @test full_mixer ≈ target[locs, locs]
    end

    @testset "ring mixer" begin
        mixer = mixer_ring(5,3)
        full_mixer = mixer.v*Diagonal(mixer.d)*mixer.vinv
        targetx = 1/2*(kron(PX,PX,I,I,I)+kron(I,PX,PX,I,I)+kron(I,I,PX,PX,I)+kron(I,I,I,PX,PX))
        targety = 1/2*(kron(PY,PY,I,I,I)+kron(I,PY,PY,I,I)+kron(I,I,PY,PY,I)+kron(I,I,I,PY,PY))
        target = targetx + targety
        locs = [8, 12, 14, 15, 20, 22, 23, 26, 27, 29]
        @test full_mixer ≈ target[locs, locs]
    end

    @testset "warm start" begin
        r = rand(2,2)
        r = (r + r')/2
        target = kron(r,I,I,I)+kron(I,r,I,I)+kron(I,I,r,I)+kron(I,I,I,r)
        rs = fill(r, 4)
        mixer = mixer_warmstart(rs)
        full_mixer = kron(mixer.v...)*Diagonal(mixer.d)*kron(mixer.vinv...)
        @test full_mixer ≈ target
    end

end