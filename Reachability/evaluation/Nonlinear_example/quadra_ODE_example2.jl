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
  F1[2, 2] = 2 * a

  F2 = zeros(2, 4) # [x, x⊗x]
  F2[1, 2] = -1
  F2[2, 1] = 2 * a
  F2[2, 4] = -2

  return F1, F2
end

# =================
# Solution method
# =================

## Solution with CARLIN
function _solve_system_carlin(; N=4, T=30.0, δ=0.1, radius0=0, bloat=false, resets=nothing, a)
  x0c = [0.1, 0.01]

  F1, F2 = system_carlin(a)
  R, Re_lambda1 = _error_bound_specabs_R(x0c, F1, F2; check=true)

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