module NMsolver

using LinearAlgebra
using ForwardDiff

export newton, newton_system

"""
    newton(f, x0; df=nothing, derivative=:finite, tol=1e-10, maxiter=100,
           damping=false, return_history=false, verbose=false)

スカラー方程式 f(x) = 0 をニュートン法で解く。

# 引数
- `f`: 解くべき関数 f(x) = 0
- `x0`: 初期値
- `df`: f の導関数（省略時は数値微分を使用）
- `derivative`: `df` を省略した場合の微分方法（`:finite` または `:forwarddiff`）
- `tol`: 収束判定の許容誤差
- `maxiter`: 最大反復回数
- `damping`: `true` の場合はバックトラックにより更新幅を調整
- `return_history`: `true` の場合は反復履歴も返す
- `verbose`: `true` の場合は各反復の情報を表示

# 返り値
- `x`: 収束解
- `iter`: 反復回数
- `converged`: 収束したかどうか
- `history`: `return_history=true` の場合のみ返す反復履歴
"""
function newton(f, x0::Real; df=nothing, derivative=:finite, tol=1e-10,
                maxiter=100, damping=false, min_alpha=1e-8,
                return_history=false, verbose=false)
    x = float(x0)
    history = NamedTuple[]

    for i in 1:maxiter
        fx = f(x)

        if abs(fx) < tol
            return _newton_return(x, length(history), true, history, return_history)
        end

        if df !== nothing
            dfx = df(x)
        elseif derivative == :forwarddiff
            dfx = ForwardDiff.derivative(f, x)
        elseif derivative == :finite
            dfx = _numerical_derivative(f, x)
        else
            error("未対応の微分方法です: derivative = $derivative")
        end

        if abs(dfx) < eps(Float64)
            error("導関数がゼロに近い: x = $x")
        end

        dx = -fx / dfx
        alpha = damping ? _line_search_scalar(f, x, fx, dx, min_alpha) : 1.0
        x = x + alpha * dx
        residual = abs(f(x))

        push!(history, (iter=i, x=x, residual=residual, step=abs(alpha * dx), alpha=alpha))
        verbose && println("iter=$i, x=$x, residual=$residual, step=$(abs(alpha * dx)), alpha=$alpha")

        if residual < tol && abs(alpha * dx) < tol
            return _newton_return(x, i, true, history, return_history)
        end
    end

    return _newton_return(x, maxiter, false, history, return_history)
end

"""
    newton_system(F, x0; J=nothing, jacobian=:finite, tol=1e-10, maxiter=100,
                  damping=false, return_history=false, verbose=false)

連立方程式 F(x) = 0 をニュートン法で解く。

# 引数
- `F`: 解くべき関数 F(x) = 0（ベクトル値）
- `x0`: 初期値（ベクトル）
- `J`: ヤコビアン行列（省略時は数値微分を使用）
- `jacobian`: `J` を省略した場合のヤコビアン計算方法（`:finite` または `:forwarddiff`）
- `tol`: 収束判定の許容誤差（ノルムで評価）
- `maxiter`: 最大反復回数
- `damping`: `true` の場合はバックトラックにより更新幅を調整
- `return_history`: `true` の場合は反復履歴も返す
- `verbose`: `true` の場合は各反復の情報を表示

# 返り値
- `x`: 収束解
- `iter`: 反復回数
- `converged`: 収束したかどうか
- `history`: `return_history=true` の場合のみ返す反復履歴
"""
function newton_system(F, x0::AbstractVector; J=nothing, jacobian=:finite,
                       tol=1e-10, maxiter=100, damping=false, min_alpha=1e-8,
                       return_history=false, verbose=false)
    x = float.(x0)
    history = NamedTuple[]

    for i in 1:maxiter
        Fx = F(x)
        residual = norm(Fx)

        if residual < tol
            return _newton_return(x, length(history), true, history, return_history)
        end

        Jx = if J !== nothing
            J(x)
        elseif jacobian == :forwarddiff
            ForwardDiff.jacobian(F, x)
        elseif jacobian == :finite
            _numerical_jacobian(F, x)
        else
            error("未対応のヤコビアン計算方法です: jacobian = $jacobian")
        end

        # ニュートンステップ: J * dx = -F(x)
        dx = Jx \ (-Fx)
        alpha = damping ? _line_search_system(F, x, residual, dx, min_alpha) : 1.0
        x = x + alpha * dx
        next_residual = norm(F(x))
        step_norm = norm(alpha * dx)

        push!(history, (iter=i, x=copy(x), residual=next_residual, step=step_norm, alpha=alpha))
        verbose && println("iter=$i, residual=$next_residual, step=$step_norm, alpha=$alpha")

        if next_residual < tol && step_norm < tol
            return _newton_return(x, i, true, history, return_history)
        end
    end

    return _newton_return(x, maxiter, false, history, return_history)
end

function _newton_return(x, iter, converged, history, return_history)
    return_history ? (x, iter, converged, history) : (x, iter, converged)
end

function _numerical_derivative(f, x)
    h = sqrt(eps(x)) * max(1.0, abs(x))
    return (f(x + h) - f(x - h)) / (2h)
end

function _numerical_jacobian(F, x)
    n = length(x)
    Fx = F(x)
    J = zeros(eltype(Fx), length(Fx), n)

    for j in 1:n
        h = sqrt(eps(x[j])) * max(1.0, abs(x[j]))
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        J[:, j] = (F(xp) - F(xm)) / (2h)
    end

    return J
end

function _line_search_scalar(f, x, fx, dx, min_alpha)
    alpha = 1.0
    current = abs(fx)

    while alpha >= min_alpha
        if abs(f(x + alpha * dx)) < current
            return alpha
        end
        alpha /= 2
    end

    return alpha
end

function _line_search_system(F, x, residual, dx, min_alpha)
    alpha = 1.0

    while alpha >= min_alpha
        if norm(F(x + alpha * dx)) < residual
            return alpha
        end
        alpha /= 2
    end

    return alpha
end

end # module NMsolver
