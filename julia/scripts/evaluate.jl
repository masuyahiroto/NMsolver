using NMsolver
using LinearAlgebra
using Statistics
using Printf

ENV["GKSwstype"] = "100"
using Plots

mkpath(joinpath(@__DIR__, "..", "results"))
const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

# ============================================================
# テスト関数の定義
# ============================================================

# 分離型二乗系: F(x)[i] = x[i]^2 - i, 解: x[i] = sqrt(i)
function make_decoupled(n::Int)
    F = v -> [v[i]^2 - Float64(i) for i in 1:n]
    J = v -> Diagonal(2 .* v)
    x_true = [sqrt(Float64(i)) for i in 1:n]
    x0 = 2.0 * ones(n)
    return F, J, x0, x_true
end

# 三対角型連立方程式, 解: x[i] = 1
function make_tridiagonal(n::Int)
    c = vcat([2.0], fill(3.0, n - 2), [2.0])
    function F(v)
        r = zeros(eltype(v), n)
        r[1] = v[1]^2 + v[2] - c[1]
        for i in 2:n-1
            r[i] = v[i-1] + v[i]^2 + v[i+1] - c[i]
        end
        r[n] = v[n-1] + v[n]^2 - c[n]
        return r
    end
    function J(v)
        M = zeros(eltype(v), n, n)
        M[1, 1] = 2v[1]
        M[1, 2] = 1
        for i in 2:n-1
            M[i, i-1] = 1
            M[i, i] = 2v[i]
            M[i, i+1] = 1
        end
        M[n, n-1] = 1
        M[n, n] = 2v[n]
        return M
    end
    x_true = ones(n)
    x0 = 0.8 * ones(n)
    return F, J, x0, x_true
end

# ============================================================
# ベンチマーク
# ============================================================

const N_TRIALS = 500
const NS = 1:10
const METHODS = (:finite, :analytic, :forwarddiff, :damped)

println("=" ^ 68)
println("NMsolver 性能評価")
println("=" ^ 68)

function solve_with(method, F, J, x0)
    if method == :finite
        return newton_system(F, x0)
    elseif method == :analytic
        return newton_system(F, x0; J=J)
    elseif method == :forwarddiff
        return newton_system(F, x0; jacobian=:forwarddiff)
    elseif method == :damped
        return newton_system(F, x0; J=J, damping=true)
    else
        error("未対応の評価方法です: $method")
    end
end

function evaluate_case(problem, n, method, F, J, x0, x_true)
    local x, iter, converged
    elapsed = minimum(@elapsed((x, iter, converged) = solve_with(method, F, J, copy(x0))) for _ in 1:N_TRIALS)
    return (
        problem=problem,
        n=n,
        method=method,
        error=norm(x - x_true),
        residual=norm(F(x)),
        iter=iter,
        converged=converged,
        time=elapsed,
    )
end

results = NamedTuple[]

for n in NS
    # --- 分離型 ---
    F, J, x0, x_true = make_decoupled(n)
    for method in METHODS
        push!(results, evaluate_case("分離型二乗系", n, method, F, J, x0, x_true))
    end

    # --- 三対角型 (n >= 2) ---
    if n >= 2
        F, J, x0, x_true = make_tridiagonal(n)
        for method in METHODS
            push!(results, evaluate_case("三対角型", n, method, F, J, x0, x_true))
        end
    end
end

# ============================================================
# 結果テーブルの出力
# ============================================================

function print_results(problem, ns)
    println("\n[$problem]")
    println("-" ^ 92)
    @printf("%-6s | %-12s | %-8s | %-16s | %-16s | %-12s\n",
        "変数数", "方法", "反復回数", "解の誤差(L2)", "残差(L2)", "時間(μs)")
    println("-" ^ 92)
    for n in ns
        for method in METHODS
            row = only(filter(r -> r.problem == problem && r.n == n && r.method == method, results))
            @printf("  %-4d | %-12s | %-8d | %-16.3e | %-16.3e | %-12.3f\n",
                row.n, String(row.method), row.iter, row.error, row.residual, row.time * 1e6)
        end
    end
end

print_results("分離型二乗系", NS)
print_results("三対角型", 2:10)

csv_path = joinpath(RESULTS_DIR, "comparison.csv")
open(csv_path, "w") do io
    println(io, "problem,n,method,converged,iter,error_l2,residual_l2,time_us")
    for row in results
        @printf(io, "%s,%d,%s,%s,%d,%.16e,%.16e,%.6f\n",
            row.problem, row.n, String(row.method), string(row.converged),
            row.iter, row.error, row.residual, row.time * 1e6)
    end
end
println("\n比較CSVを保存: results/comparison.csv")

# ============================================================
# グラフの作成
# ============================================================

ns_vec     = collect(NS)
ns_tri_vec = collect(2:10)

function series(problem, method, field, ns)
    return [getfield(only(filter(r -> r.problem == problem && r.n == n && r.method == method, results)), field) for n in ns]
end

plot_error(values) = max.(values, eps(Float64))

# 1. 誤差 vs 変数数
p1 = plot(
    xlabel  = "number of variables n",
    ylabel  = "solution error (L2 norm)",
    title   = "Decoupled quadratic system: error comparison",
    yscale  = :log10,
    lw = 2, ms = 6)
for method in METHODS
    plot!(p1, ns_vec, plot_error(series("分離型二乗系", method, :error, NS)),
        label=String(method), marker=:circle, lw=2, ms=6)
end
savefig(p1, joinpath(RESULTS_DIR, "error_vs_n.png"))
println("\n誤差グラフを保存: results/error_vs_n.png")

# 2. 反復回数 vs 変数数
p2 = plot(
    xlabel  = "number of variables n",
    ylabel  = "iterations",
    title   = "Decoupled quadratic system: iteration comparison",
    lw = 2, ms = 6)
for method in METHODS
    plot!(p2, ns_vec, series("分離型二乗系", method, :iter, NS),
        label=String(method), marker=:circle, lw=2, ms=6)
end
savefig(p2, joinpath(RESULTS_DIR, "iters_vs_n.png"))
println("反復回数グラフを保存: results/iters_vs_n.png")

# 3. 計算時間 vs 変数数
p3 = plot(
    xlabel  = "number of variables n",
    ylabel  = "time (microseconds)",
    title   = "Decoupled quadratic system: time comparison",
    lw = 2, ms = 6)
for method in METHODS
    plot!(p3, ns_vec, series("分離型二乗系", method, :time, NS) .* 1e6,
        label=String(method), marker=:circle, lw=2, ms=6)
end
savefig(p3, joinpath(RESULTS_DIR, "time_vs_n.png"))
println("計算時間グラフを保存: results/time_vs_n.png")

# 4. 三対角型の誤差比較
p4 = plot(
    xlabel  = "number of variables n",
    ylabel  = "solution error (L2 norm)",
    title   = "Tridiagonal system: error comparison",
    yscale  = :log10,
    lw = 2, ms = 6)
for method in METHODS
    plot!(p4, ns_tri_vec, plot_error(series("三対角型", method, :error, 2:10)),
        label=String(method), marker=:square, lw=2, ms=6)
end
savefig(p4, joinpath(RESULTS_DIR, "tridiagonal_error_vs_n.png"))
println("三対角型誤差グラフを保存: results/tridiagonal_error_vs_n.png")

# 5. 全グラフ合成
p_all = plot(p1, p2, p3, p4, layout = (4, 1), size = (800, 1200))
savefig(p_all, joinpath(RESULTS_DIR, "summary.png"))
println("サマリーグラフを保存: results/summary.png")

println("\n評価完了!")
