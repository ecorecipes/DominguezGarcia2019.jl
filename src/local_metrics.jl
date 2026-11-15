# Stability metrics that need only the Jacobian (community matrix) J at an equilibrium, and for
# the community-level ("median") versions the equilibrium biomasses. Definitions follow SI §10
# of Domínguez-García et al. (2019), Neubert & Caswell (1997) and Arnoldi et al. (2016, 2018).
# All functions return the raw quantity, in the sign convention of the dataset's raw columns;
# `prepare_metrics` turns them into the paper's "larger is more stable" versions.

"""
    reactivity(J)

Maximal instantaneous amplification rate of a pulse, `λ_max((J + Jᵀ)/2)` (Neubert & Caswell
1997). Positive values mean a reactive equilibrium. The dataset's raw `R0` column holds this
value; SI eq. 6 writes its negative, which the analysis applies later.
"""
reactivity(J::AbstractMatrix) = eigmax(Symmetric((J + J') / 2))

"""
    resilience(J)

Asymptotic return rate `−Re λ_dom(J)` (SI eq. 10); positive for a stable equilibrium.
"""
resilience(J::AbstractMatrix) = -maximum(real, eigvals(J))

"""
    amplification_envelope(J, t; p = 2)

`A(t) = ‖exp(J t)‖ₚ`, the operator `p`-norm of the propagator (SI eq. 7). With `p = 2` it is the
largest amplification at time `t` of any pulse; `p = Inf` (the maximum absolute row sum) is the
norm behind the published `Amax` (see [`community_record`](@ref)).
"""
amplification_envelope(J::AbstractMatrix, t::Real; p::Real = 2) = opnorm(exp(J * t), p)

"""
    max_amplification(J; dt = 0.1, tcap = 1e4, p = 2)

The maximum of the amplification envelope ([`amplification_envelope`](@ref), operator `p`-norm)
and the time at which it occurs, `(Amax, tmax)` (SI eqs. 8–9), on a time grid of step `dt` (the
dataset's `tmax` values are multiples of 0.1). A non-reactive equilibrium (`R0 ≤ 0`) gives
`(1.0, 0.0)` for every `p`, as in the published data, where `Amax = 1` exactly when `R0 ≤ 0`.
Otherwise the grid is followed until the envelope falls back below 1, or to `tcap`.
"""
function max_amplification(J::AbstractMatrix; dt::Real = 0.1, tcap::Real = 1e4, p::Real = 2)
    reactivity(J) ≤ 0 && return (1.0, 0.0)
    E = exp(J * dt)
    Φ = Matrix(E)
    best, tbest, k = 1.0, 0.0, 1
    while k * dt ≤ tcap
        A = opnorm(Φ, p)
        A > best && ((best, tbest) = (A, k * dt))
        A < 1 && break
        Φ = Φ * E
        k += 1
    end
    return (best, round(tbest; digits = 10))
end

"""
    pulse_covariance(B)

The covariance of pulse perturbations proportional to equilibrium biomasses, `Diagonal(B.^2)`,
the choice of Domínguez-García et al. (2019) following Arnoldi et al. (2018).
"""
pulse_covariance(B::AbstractVector) = Diagonal(B .^ 2)

"""
    median_reactivity(J, C)

The median initial rate of displacement, `Tr(C J)/Tr(C)` for perturbation covariance `C`
(SI eq. 11; Arnoldi et al. 2018).
"""
median_reactivity(J::AbstractMatrix, C::AbstractMatrix) = tr(C * J) / tr(C)

"""
    median_displacement(J, C, t)

The median displacement at time `t` relative to its initial value,
`sqrt(Tr(C exp(Jᵀt) exp(Jt)) / Tr(C))` (SI eq. 12, as a norm ratio).
"""
function median_displacement(J::AbstractMatrix, C::AbstractMatrix, t::Real)
    E = exp(J * t)
    return sqrt(max(tr(C * (E' * E)), 0.0) / tr(C))
end

"""
    median_max_amplification(J, C; dt = 0.1, threshold = 0.005, tcap = 1e4)

`(MAmax, Mtmax)`: the maximum over a grid of step `dt` of the
[`median_displacement`](@ref) and the time it is reached, or `(1.0, 0.0)` if the displacement
never exceeds its initial value. The scan runs until the displacement falls below `threshold`
(the SI's asymptotic criterion, as in [`median_resilience`](@ref)) or `t > tcap`. It does not
stop when the displacement first falls: a median displacement can shrink at first and grow
later, even when the median reactivity is negative.
"""
function median_max_amplification(J::AbstractMatrix, C::AbstractMatrix; dt::Real = 0.1, threshold::Real = 0.005,
                                  tcap::Real = 1e4)
    E = exp(J * dt)
    Φ = Matrix(E)
    best, tbest, k, trC = 1.0, 0.0, 1, tr(C)
    while k * dt ≤ tcap
        M = sqrt(max(tr(C * (Φ' * Φ)), 0.0) / trC)
        M > best && ((best, tbest) = (M, k * dt))
        M < threshold && break
        Φ = Φ * E
        k += 1
    end
    return (best, round(tbest; digits = 10))
end

"""
    median_resilience(J, C; threshold = 0.005, dt = 0.1)

The median return rate `−ln(Tr(C exp(Jᵀt) exp(Jt))/Tr(C)) / (2t)` (SI eq. 13), evaluated at the
first time `t` on a grid of step `dt` at which the [`median_displacement`](@ref) is below
`threshold` (SI: "the asymptotic limit is reached if the displacement is below 0.005"). The time
is found by doubling and then bisection on the grid. Returns `NaN` for an unstable `J`.
"""
function median_resilience(J::AbstractMatrix, C::AbstractMatrix; threshold::Real = 0.005, dt::Real = 0.1)
    resilience(J) > 0 || return NaN
    below(k) = median_displacement(J, C, k * dt) < threshold
    hi = 1
    while !below(hi)
        hi *= 2
        hi > 10^12 && return NaN
    end
    lo = hi ÷ 2
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        below(mid) ? (hi = mid) : (lo = mid)
    end
    t = hi * dt
    return -log(median_displacement(J, C, t)^2) / (2t)
end

"""
    stochastic_invariability(J; dense_max = 30, tol = 1e-10, maxiter = 2000)

Intrinsic stochastic invariability `Is = 1 / (2 ‖Ĵ⁻¹‖₂)` with `Ĵ = J ⊗ I + I ⊗ J` (SI eq. 14;
Arnoldi et al. 2016). For `size(J, 1) ≤ dense_max` the Kronecker sum is formed and its smallest
singular value computed directly. Otherwise `‖Ĵ⁻¹‖₂` is found by power iteration on
`Ĵ⁻ᵀĴ⁻¹`, each product a pair of Sylvester solves (`Ĵ vec X = vec(JX + XJᵀ)`), so memory stays
O(N²).
"""
function stochastic_invariability(J::AbstractMatrix; dense_max::Integer = 30, tol::Real = 1e-10,
                                  maxiter::Integer = 2000)
    n = size(J, 1)
    if n ≤ dense_max
        In = Matrix{Float64}(I, n, n)
        Jh = kron(J, In) + kron(In, J)
        return minimum(svdvals(Jh)) / 2
    end
    Jt = Matrix(J')
    X = ones(n, n) ./ n
    λ = 0.0
    for _ in 1:maxiter
        Y = sylvester(J, Jt, -X)        # J Y + Y Jᵀ = X
        Z = sylvester(Jt, J, -Y)        # Jᵀ Z + Z J = Y
        λnew = norm(Z)
        X = Z ./ λnew
        abs(λnew - λ) ≤ tol * λnew && (λ = λnew; break)
        λ = λnew
    end
    return 1 / (2sqrt(λ))
end

"""
    sensitivity_matrix(J)

The sensitivity matrix `−J⁻¹`: the equilibrium response of every species to a small sustained
press on each (Bender et al. 1984).
"""
sensitivity_matrix(J::AbstractMatrix) = -inv(J)

"""
    press_sensitivity(J)

Community sensitivity `S = Σᵢ |Σⱼ sᵢⱼ|` of the sensitivity matrix (SI eq. 15); the raw column
`L1DBRMt_s`.
"""
press_sensitivity(J::AbstractMatrix) = sum(abs, sum(sensitivity_matrix(J); dims = 2))

"""
    total_sensitivity(J)

`Σᵢⱼ |sᵢⱼ|`, the raw column `absumS`. SI eq. 16's mean element strength `<sᵢⱼ>` divides it by
N², which `prepare_metrics` does.
"""
total_sensitivity(J::AbstractMatrix) = sum(abs, sensitivity_matrix(J))
