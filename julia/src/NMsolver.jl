module NMsolver

using LinearAlgebra

export newton, newton_system

"""
    newton(f, x0; df=nothing, tol=1e-10, maxiter=100) -> (x, iter, converged)

スカラー方程式 f(x) = 0 をニュートン法で解く。

# 引数
- `f`: 解くべき関数 f(x) = 0
- `x0`: 初期値
- `df`: f の導関数（省略時は数値微分を使用）
- `tol`: 収束判定の許容誤差
- `maxiter`: 最大反復回数

# 返り値
- `x`: 収束解
- `iter`: 反復回数
- `converged`: 収束したかどうか
"""
function newton(f, x0::Real; df=nothing, tol=1e-10, maxiter=100)
    x = float(x0)

    for i in 1:maxiter
        fx = f(x)

        if abs(fx) < tol
            return x, i, true
        end

        if df !== nothing
            dfx = df(x)
        else
            h = sqrt(eps(x)) * max(1.0, abs(x))
            dfx = (f(x + h) - f(x - h)) / (2h)
        end

        if abs(dfx) < eps(Float64)
            error("導関数がゼロに近い: x = $x")
        end

        x = x - fx / dfx
    end

    return x, maxiter, false
end

"""
    newton_system(F, x0; J=nothing, tol=1e-10, maxiter=100) -> (x, iter, converged)

連立方程式 F(x) = 0 をニュートン法で解く。

# 引数
- `F`: 解くべき関数 F(x) = 0（ベクトル値）
- `x0`: 初期値（ベクトル）
- `J`: ヤコビアン行列（省略時は数値微分を使用）
- `tol`: 収束判定の許容誤差（ノルムで評価）
- `maxiter`: 最大反復回数

# 返り値
- `x`: 収束解
- `iter`: 反復回数
- `converged`: 収束したかどうか
"""
function newton_system(F, x0::AbstractVector; J=nothing, tol=1e-10, maxiter=100)
    x = float.(x0)

    for i in 1:maxiter
        Fx = F(x)

        if norm(Fx) < tol
            return x, i, true
        end

        Jx = J !== nothing ? J(x) : _numerical_jacobian(F, x)

        # ニュートンステップ: J * dx = -F(x)
        dx = Jx \ (-Fx)
        x = x + dx
    end

    return x, maxiter, false
end

function _numerical_jacobian(F, x)
    n = length(x)
    J = zeros(n, n)

    for j in 1:n
        h = sqrt(eps(x[j])) * max(1.0, abs(x[j]))
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        J[:, j] = (F(xp) - F(xm)) / (2h)
    end

    return J
end

end # module NMsolver
