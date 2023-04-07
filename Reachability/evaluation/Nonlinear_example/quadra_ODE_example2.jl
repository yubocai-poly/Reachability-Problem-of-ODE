# =================
# Dependencies
# =================

using ReachabilityAnalysis, CarlemanLinearization
using Plots, LaTeXStrings, LinearAlgebra, SparseArrays
using Plots.PlotMeasures
using LazySets: center
using CarlemanLinearization: _error_bound_specabs_R

include("../utils.jl")

# =================
# Model definition (The second example)
# =================

function system_carlin(a)

  F1 = zeros(2, 2)
  F1[1, 1] = -a
  F1[2, 2] = -2 * a

  F2 = zeros(2, 4) # [x, x⊗x]
  F2[1, 2] = -1
  F2[2, 1] = -2 * a
  F2[2, 4] = -2

  return F1, F2
end

# =================
# Changing Condition function
# =================
# See Definition (2.2) in [2]. These bounds use the spectral norm (p = 2)
function _error_bound_specabs_R_(x₀, F₁, F₂; check=true)
  nx₀ = norm(x₀, 2)
  nF₂ = opnorm(F₂, 2)

  # compute eigenvalues and sort them by increasing real part
  λ = eigvals(F₁, sortby=real)
  λ₁ = last(λ)
  Re_λ₁ = real(λ₁)
  if check
      @assert Re_λ₁ <= 0 "expected Re(λ₁) ≤ 0, got $Re_λ₁"
  end
  R = nx₀ * nF₂ / abs(Re_λ₁)
  return (R, Re_λ₁, nx₀, nF₂)
end


# =================
# Solution method
# =================

## Solution with CARLIN
function _solve_system_carlin(; N=4, T=30.0, δ=0.1, radius0=0, bloat=false, resets=nothing, a, x0)
  x0c = [x0, x0^2 - a]

  F1, F2 = system_carlin(a)
  R, Re_lambda1 = _error_bound_specabs_R_(x0c, F1, F2; check=true)

  n = 2
  dirs = _template(n=n, N=N)
  alg = LGG09(δ=δ, template=dirs)

  if radius0 == 0
    X0 = convert(Hyperrectangle, Singleton(x0c))
  else
    X0 = Hyperrectangle(x0c, radius0)
  end

  if isnothing(resets)
    @time sol = _solve_CARLIN(X0, F1, F2; alg=alg, N=N, T=T, bloat=bloat)
  else
    @time sol = _solve_CARLIN_resets(X0, F1, F2; resets=resets, alg=alg, N=N, T=T, bloat=bloat)
  end

  return sol
end

@taylorize function system_equation(dx, x, a)
  x1, x2 = x # y = x2

  dx[1] = - x1 * x2 - a * x1
  dx[2] = -2 * x2^2 - 2 * a * x1^2 - 2 * a * x2

end

# Ploting the results
Tmax = 10.0
rr0 = 0.0

# figure with error bounds with comparison
solN4_a1_bloat = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=true, a=1.0, x0=0.5)

function figure_System_withError()

  fig = plot(legend=:topright, xlab = L"\textrm{Time t}", ylab = L"\textrm{x(t) or y(t)} ", title="Numerial solution of the system of equations of x(t) (with error bounds, \n with comparison of different values of a)",
  legendfontsize=12,
  tickfont=font(10, "Times"),
  guidefontsize=10,
  xguidefont=font(10, "Times"),
  yguidefont=font(10, "Times"),
  xlims= (1.0, 6.0),
  ylims = (-0.0, 0.2),
  bottom_margin=5mm,
  left_margin=5mm,
  right_margin=5mm,
  top_margin=5mm,
  size=(800, 600))

  # carleman linearization solution with error bounds
  plot!(fig, solN4_a3_bloat,  vars=(0, 1), color=:lightblue, lc=:lightblue, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=3.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a2_bloat,  vars=(0, 1), color=:lightgrey, lc=:lightgrey, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=2.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a1_bloat,  vars=(0, 1), color=:green, lc=:green, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=4, } \alpha=1.0")

  # carleman linearization solution 
  plot!(fig, solN4_a1,  vars=(0, 1), color=:orange, lc=:orange, linewidth=2, label=L"x'=-x+axy \textrm{, a=1.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a2,  vars=(0, 1), color=:red, lc=:red, linewidth=2, label=L"y'=-2y+2ay^{2} \textrm{, a=2.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a1,  vars=(0, 1), color=:blue, lc=:blue, linewidth=2, label=L"x'=-x+axy \textrm{, a=3.0, N=4, } \alpha=1.0")

return fig
end