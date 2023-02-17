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
# Model definition
# =================

function system_carlin(alpha)

  a = 1.0
  F1 = zeros(2, 2)
  F1[1, 1] = -1
  F1[2, 2] = -2 / alpha

  F2 = zeros(2, 4) # [x, x⊗x]
  F2[1, 4] = 2 * a / 4
  F2[2, 2] = a

  print("F1 = ", F1, '\n')
  return F1, F2
end

# =================
# Solution method
# =================

## Solution with CARLIN


function _solve_system_carlin(; N=4, T=30.0, δ=0.1, radius0=0, bloat=false, resets=nothing, alpha)
  x0c = [0.1, 0.01]

  F1, F2 = system_carlin(alpha)
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

## Solution with TMJETS

@taylorize function system_equation(dx, x, alpha)
  a = 1.0
  x1, x2 = x # y = x2

  dx[1] = -2 / alpha * x2 + 2 * a * x2^2 / (alpha) ^ 2
  dx[2] = -x1 + a * x1 * x2

end

# function _solve_system_carlin_TM(; T=30.0, radius0=0, trajectories=-1)
#   x0c = [0.0, 0.0]

#   if radius0 == 0
#     X0 = convert(Hyperrectangle, Singleton(x0c))
#   else
#     X0 = Hyperrectangle(x0c, radius0)
#   end

#   prob = @ivp(x' = system_equation(x), x(0) ∈ X0, dim=2)

#   if trajectories == -1
#     sol = solve(prob, T=T, alg=TMJets())
#   else
#     sol = solve(prob, T=T, alg=TMJets(), trajectories=trajectories)
#   end
  

#   return sol
# end

# ===============
# Results
# ===============

# parameters
Tmax = 20.0
rr0 = 0.0

# # taylor models solution
# _solve_system_carlin_TM(T=30.0, radius0=rr0, trajectories=-1)

# no error bounds, N = 2
_solve_system_carlin(N=2, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=2.0)
time_NoError_N2_alpha1 = @elapsed _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=2.0)

_solve_system_carlin(N=2, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=5.0)
time_NoError_N2_alpha2 = @elapsed _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=5.0)

_solve_system_carlin(N=2, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=10.0)
time_NoError_N2_alpha3 = @elapsed _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=10.0)

# including error bounds, N = 5
# _solve_system_carlin(N=5, T=Tmax, δ=0.1, radius0=rr0, bloat=true, alpha=2.0)
# time_Error_N5 = @elapsed _solve_system_carlin(N=5, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=2.0)

print("result of N=2, Alpha=2.0, No Error: ", (time_NoError_N2_alpha1), '\n')
print("result of N=2, Alpha=5.0, No Error: ", (time_NoError_N2_alpha2), '\n')
print("result of N=2, Alpha=10.0, No Error: ", (time_NoError_N2_alpha3), '\n')
# print("result of N=5, Error", (time_Error_N5), '\n')
# print("result of TM", (time_TM), '\n')

# print(io, "The first Quadratization Model, Carleman, no error bound, N=2, $(time_NoError_N2)\n")
# print(io, "The first Quadratization Model, Carleman, error bound, N=5, $(time_Error_N5)\n")
# print(io, "The first Quadratization Model, Taylor Models, $(time_TM)\n")

# figure with NO error bounds
function figure_System_NoError()

  Tmax = 10.0
  rr0 = 0.0
  solN2_alpha1 = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=2.0)
  solN2_alpha2 = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=5.0)
  solN2_alpha3 = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=10.0)
  # solN4 = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=false)
  # solN6 = _solve_system_carlin(N=6, T=Tmax, δ=0.1, radius0=rr0, bloat=false)

  fig = plot(legend=:topright, xlab = L"x", ylab = L"y",
              legendfontsize=25,
              tickfont=font(25, "Times"),
              guidefontsize=25,
              xguidefont=font(15, "Times"),
              yguidefont=font(15, "Times"),
              bottom_margin=5mm,
              left_margin=5mm,
              right_margin=5mm,
              top_margin=5mm,
              size=(800, 600))
  

  plot!(fig, solN2_alpha1,  vars=(0, 2), color=:red, linewidth=2, label="Second variable Alpha = 2.0")
  plot!(fig, solN2_alpha2,  vars=(0, 2), color=:blue, linewidth=2, label="Second variable Alpha = 5.0")
  plot!(fig, solN2_alpha3,  vars=(0, 2), color=:green, linewidth=2, label="Second variable Alpha = 10.0")
  

  return fig
              
end

fig = figure_System_NoError()
# savefig(fig, joinpath(TARGET_FOLDER, "figure_1a_non error.pdf"))

# figure with error bounds
function figure_System_withError()

  Tmax = 10.0
  rr0 = 0.0
  solN4_alpha1 = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=2.0)
  solN4_alpha2 = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=5.0)
  solN4_alpha3 = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=10.0)

  solN4_alpha1_bloat = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=true, resets=[4.0], alpha=2.0)
  solN4_alpha2_bloat = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=true, resets=[4.0], alpha=5.0)
  solN4_alpha3_bloat = _solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=true, resets=[4.0], alpha=10.0)

  fig = plot(legend=:topright, xlab = L"x", ylab = L"y",
              legendfontsize=10,
              tickfont=font(25, "Times"),
              guidefontsize=25,
              xguidefont=font(15, "Times"),
              yguidefont=font(15, "Times"),
              bottom_margin=5mm,
              left_margin=5mm,
              right_margin=5mm,
              top_margin=5mm,
              size=(800, 600))

  # carleman linearization solution

  # plot!(fig, solN4_alpha1,  vars=(0, 1), color=:red, linewidth=2, label="Second variable Alpha = 2.0")

  # plot!(fig, solN4_alpha2,  vars=(0, 1), color=:blue, linewidth=2, label="Second variable Alpha = 5.0")

  plot!(fig, solN4_alpha3,  vars=(0, 1), color=:green, linewidth=2, label="Second variable Alpha = 10.0")

  # carleman linearization solution with error bounds

  # plot!(fig, solN4_alpha1_bloat,  vars=(0, 1), color=:aquamarine, linewidth=2, linestyle=:dash, label="Second variable Alpha = 2.0")

  # plot!(fig, solN4_alpha2_bloat,  vars=(0, 1), color=:darksalmon, linewidth=2, linestyle=:dash, label="Second variable Alpha = 5.0")

  plot!(fig, solN4_alpha3_bloat,  vars=(0, 1), color=:darkseagreen, linewidth=2, linestyle=:dash, label="Second variable Alpha = 10.0")

  return fig
end

# fig = figure_System_withError()
# savefig(fig, joinpath(TARGET_FOLDER, "figure_1b_error.pdf"))


