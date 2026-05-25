using Test
using NMsolver
using LinearAlgebra

@testset "NMsolver.jl" begin

    @testset "スカラーニュートン法 (解析的導関数)" begin
        # x^2 - 2 = 0 → x = sqrt(2)
        f = x -> x^2 - 2
        df = x -> 2x
        x, iter, conv = newton(f, 1.0; df=df)
        @test conv
        @test abs(x - sqrt(2)) < 1e-8
        @test iter < 20

        # x^3 - x - 2 = 0 → x ≈ 1.5214
        g = x -> x^3 - x - 2
        dg = x -> 3x^2 - 1
        x2, _, conv2 = newton(g, 1.5; df=dg)
        @test conv2
        @test abs(g(x2)) < 1e-9
    end

    @testset "スカラーニュートン法 (数値微分)" begin
        # x^2 - 2 = 0 → x = sqrt(2)（導関数なし）
        f = x -> x^2 - 2
        x, iter, conv = newton(f, 1.0)
        @test conv
        @test abs(x - sqrt(2)) < 1e-8

        # sin(x) = 0 → x ≈ π (初期値 3.0)
        x2, _, conv2 = newton(sin, 3.0)
        @test conv2
        @test abs(x2 - π) < 1e-8

        # exp(x) - 3 = 0 → x = log(3)
        h = x -> exp(x) - 3
        x3, _, conv3 = newton(h, 1.0)
        @test conv3
        @test abs(x3 - log(3)) < 1e-8
    end

    @testset "収束しない場合" begin
        # 最大反復回数を超えた場合
        f = x -> x^2 - 2
        x, iter, conv = newton(f, 1.0; maxiter=2)
        @test !conv
        @test iter == 2
    end

    @testset "連立方程式ニュートン法 (数値微分)" begin
        # x^2 + y^2 = 1, x = y → (1/√2, 1/√2)
        F = v -> [v[1]^2 + v[2]^2 - 1, v[1] - v[2]]
        x, iter, conv = newton_system(F, [1.0, 0.5])
        @test conv
        @test abs(x[1] - 1/sqrt(2)) < 1e-8
        @test abs(x[2] - 1/sqrt(2)) < 1e-8

        # x + y = 3, x*y = 2 → (1, 2) または (2, 1)
        G = v -> [v[1] + v[2] - 3, v[1]*v[2] - 2]
        x2, _, conv2 = newton_system(G, [0.5, 2.5])
        @test conv2
        @test abs(G(x2)[1]) < 1e-9
        @test abs(G(x2)[2]) < 1e-9
    end

    @testset "連立方程式ニュートン法 (解析的ヤコビアン)" begin
        F = v -> [v[1]^2 + v[2]^2 - 1, v[1] - v[2]]
        Jfunc = v -> [2v[1] 2v[2]; 1.0 -1.0]

        x, iter, conv = newton_system(F, [1.0, 0.5]; J=Jfunc)
        @test conv
        @test abs(x[1] - 1/sqrt(2)) < 1e-8
        @test abs(x[2] - 1/sqrt(2)) < 1e-8
        @test iter < 20
    end

    @testset "3変数連立方程式" begin
        # x + y + z = 6, x^2 + y^2 + z^2 = 14, x*y*z = 6 → (1, 2, 3)
        F = v -> [
            v[1] + v[2] + v[3] - 6,
            v[1]^2 + v[2]^2 + v[3]^2 - 14,
            v[1]*v[2]*v[3] - 6
        ]
        x, _, conv = newton_system(F, [0.5, 1.5, 3.5])
        @test conv
        @test norm(F(x)) < 1e-9
    end

    # 4〜10変数: F(x)[i] = x[i]^2 - i, 解: x[i] = sqrt(i)
    @testset "$(n)変数連立方程式" for n in 4:10
        F = v -> [v[i]^2 - Float64(i) for i in 1:n]
        x0 = 2.0 * ones(n)
        x, _, conv = newton_system(F, x0)
        @test conv
        @test norm(F(x)) < 1e-9
        @test norm(x - [sqrt(Float64(i)) for i in 1:n]) < 1e-8
    end

    # 10変数: 三対角型連立方程式, 解: x[i] = 1
    @testset "10変数三対角型連立方程式" begin
        n = 10
        c = vcat([2.0], fill(3.0, n-2), [2.0])
        function F_tri(v)
            r = zeros(n)
            r[1] = v[1]^2 + v[2] - c[1]
            for i in 2:n-1
                r[i] = v[i-1] + v[i]^2 + v[i+1] - c[i]
            end
            r[n] = v[n-1] + v[n]^2 - c[n]
            return r
        end
        x, _, conv = newton_system(F_tri, 0.5 * ones(n))
        @test conv
        @test norm(F_tri(x)) < 1e-9
        @test norm(x - ones(n)) < 1e-8
    end

end
