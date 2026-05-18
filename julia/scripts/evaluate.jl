using NMsolver
using LinearAlgebra
using Statistics
using Printf
using Plots

mkpath(joinpath(@__DIR__, "..", "results"))
const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

# ============================================================
# テスト関数の定義
# ============================================================

# 分離型二乗系: F(x)[i] = x[i]^2 - i, 解: x[i] = sqrt(i)
function make_decoupled(n::Int)
    F = v -> [v[i]^2 - Float64(i) for i in 1:n]
    x_true = [sqrt(Float64(i)) for i in 1:n]
    x0 = 2.0 * ones(n)
    return F, x0, x_true
end

# 三対角型連立方程式, 解: x[i] = 1
function make_tridiagonal(n::Int)
    c = vcat([2.0], fill(3.0, n - 2), [2.0])
    function F(v)
        r = zeros(n)
        r[1] = v[1]^2 + v[2] - c[1]
        for i in 2:n-1
            r[i] = v[i-1] + v[i]^2 + v[i+1] - c[i]
        end
        r[n] = v[n-1] + v[n]^2 - c[n]
        return r
    end
    x_true = ones(n)
    x0 = 0.5 * ones(n)
    return F, x0, x_true
end

# ============================================================
# ベンチマーク
# ============================================================

const N_TRIALS = 500
const NS = 1:10

println("=" ^ 68)
println("NMsolver 性能評価")
println("=" ^ 68)

errors_dec  = Float64[]
iters_dec   = Int[]
times_dec   = Float64[]

errors_tri  = Float64[]
iters_tri   = Int[]
times_tri   = Float64[]

for n in NS
    # --- 分離型 ---
    if n == 1
        f1 = x -> x^2 - 1.0
        local x1, it1, cv1
        t1 = minimum(@elapsed((x1, it1, cv1) = newton(f1, 2.0)) for _ in 1:N_TRIALS)
        push!(errors_dec, abs(x1 - 1.0))
        push!(iters_dec,  it1)
        push!(times_dec,  t1)
    else
        F, x0, x_true = make_decoupled(n)
        local x2, it2, cv2
        t2 = minimum(@elapsed((x2, it2, cv2) = newton_system(F, copy(x0))) for _ in 1:N_TRIALS)
        push!(errors_dec, norm(x2 - x_true))
        push!(iters_dec,  it2)
        push!(times_dec,  t2)
    end

    # --- 三対角型 (n >= 2) ---
    if n >= 2
        F, x0, x_true = make_tridiagonal(n)
        local x3, it3, cv3
        t3 = minimum(@elapsed((x3, it3, cv3) = newton_system(F, copy(x0))) for _ in 1:N_TRIALS)
        push!(errors_tri, norm(x3 - x_true))
        push!(iters_tri,  it3)
        push!(times_tri,  t3)
    end
end

# ============================================================
# 結果テーブルの出力
# ============================================================

println("\n[分離型二乗系: F(x)[i] = x[i]^2 - i]")
println("-" ^ 68)
@printf("%-6s | %-8s | %-16s | %-12s\n", "変数数", "反復回数", "解の誤差(L2ノルム)", "計算時間(μs)")
println("-" ^ 68)
for (i, n) in enumerate(NS)
    @printf("  %-4d | %-8d | %-16.3e | %-12.3f\n",
        n, iters_dec[i], errors_dec[i], times_dec[i] * 1e6)
end

println("\n[三対角型連立方程式: x[i] = 1 が解]")
println("-" ^ 68)
@printf("%-6s | %-8s | %-16s | %-12s\n", "変数数", "反復回数", "解の誤差(L2ノルム)", "計算時間(μs)")
println("-" ^ 68)
for (i, n) in enumerate(2:10)
    @printf("  %-4d | %-8d | %-16.3e | %-12.3f\n",
        n, iters_tri[i], errors_tri[i], times_tri[i] * 1e6)
end

# ============================================================
# グラフの作成
# ============================================================

ns_vec     = collect(NS)
ns_tri_vec = collect(2:10)

# 1. 誤差 vs 変数数
p1 = plot(ns_vec, errors_dec,
    xlabel  = "変数数 n",
    ylabel  = "解の誤差 (L2ノルム)",
    title   = "ニュートン法: 解の誤差 vs 変数数",
    label   = "分離型二乗系",
    marker  = :circle,
    yscale  = :log10,
    lw = 2, ms = 6)
plot!(p1, ns_tri_vec, errors_tri,
    label  = "三対角型",
    marker = :square,
    lw = 2, ms = 6)
savefig(p1, joinpath(RESULTS_DIR, "error_vs_n.png"))
println("\n誤差グラフを保存: results/error_vs_n.png")

# 2. 反復回数 vs 変数数
p2 = plot(ns_vec, iters_dec,
    xlabel  = "変数数 n",
    ylabel  = "反復回数",
    title   = "ニュートン法: 反復回数 vs 変数数",
    label   = "分離型二乗系",
    marker  = :circle,
    lw = 2, ms = 6)
plot!(p2, ns_tri_vec, iters_tri,
    label  = "三対角型",
    marker = :square,
    lw = 2, ms = 6)
savefig(p2, joinpath(RESULTS_DIR, "iters_vs_n.png"))
println("反復回数グラフを保存: results/iters_vs_n.png")

# 3. 計算時間 vs 変数数
p3 = plot(ns_vec, times_dec .* 1e6,
    xlabel  = "変数数 n",
    ylabel  = "計算時間 (μs)",
    title   = "ニュートン法: 計算時間 vs 変数数",
    label   = "分離型二乗系",
    marker  = :circle,
    lw = 2, ms = 6)
plot!(p3, ns_tri_vec, times_tri .* 1e6,
    label  = "三対角型",
    marker = :square,
    lw = 2, ms = 6)
savefig(p3, joinpath(RESULTS_DIR, "time_vs_n.png"))
println("計算時間グラフを保存: results/time_vs_n.png")

# 4. 全グラフ合成
p_all = plot(p1, p2, p3, layout = (3, 1), size = (700, 900))
savefig(p_all, joinpath(RESULTS_DIR, "summary.png"))
println("サマリーグラフを保存: results/summary.png")

println("\n評価完了!")
