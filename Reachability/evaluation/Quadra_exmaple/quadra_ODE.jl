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

function system_carlin(a, alpha)

  F1 = zeros(2, 2)
  F1[1, 1] = -1
  F1[2, 2] = -2 / alpha

  F2 = zeros(2, 4) # [x, x⊗x]
  F2[1, 2] = 2 * a / alpha
  F2[2, 4] = a

  print("F1 = ", F1, '\n')
  print("F2 = ", F2, '\n')
  return F1, F2
end

# =================
# Solution method
# =================

## Solution with CARLIN


function _solve_system_carlin(; N=4, T=30.0, δ=0.1, radius0=0, bloat=false, resets=nothing, alpha, a)
  x0c = [0.1, 0.1]

  F1, F2 = system_carlin(a, alpha)
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

@taylorize function system_equation(dx, x, alpha, a)
  x1, x2 = x # y = x2

  dx[1] = -x1 + a * x1 * x2
  dx[2] = -2 / alpha * x2 + 2 * a * x2^2 / (alpha) ^ 2

end

function _solve_system_carlin_TM(; T=30.0, radius0=0, trajectories=-1)
  x0c = [0.0, 0.0]

  if radius0 == 0
    X0 = convert(Hyperrectangle, Singleton(x0c))
  else
    X0 = Hyperrectangle(x0c, radius0)
  end

  prob = @ivp(x' = system_equation(x), x(0) ∈ X0, dim=2)

  if trajectories == -1
    sol = solve(prob, T=T, alg=TMJets())
  else
    sol = solve(prob, T=T, alg=TMJets(), trajectories=trajectories)
  end
  

  return sol
end

# ===============
# Results
# ===============

# parameters
Tmax = 20.0
rr0 = 0.0

# # taylor models solution
# _solve_system_carlin_TM(T=30.0, radius0=rr0, trajectories=-1)

# no error bounds, N = 2
_solve_system_carlin(N=2, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=1.0, a=1.0)
time_NoError_N2_a1 = @elapsed _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=2.0, a=1.0)

_solve_system_carlin(N=2, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=1.0, a=5.0)
time_NoError_N2_a2 = @elapsed _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=2.0, a=5.0)

# no error bounds, N = 4
_solve_system_carlin(N=4, T=Tmax, δ=0.1, radius0=rr0, bloat=false, alpha=1.0, a=5.0)
time_NoError_N4_a2 = @elapsed _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=2.0, a=5.0)

# including error bounds, N = 5
_solve_system_carlin(N=5, T=Tmax, δ=0.1, radius0=rr0, bloat=true, alpha=1.0, a=1.0)
time_Error_N5 = @elapsed _solve_system_carlin(N=5, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=1.0, a=1.0)

print(io, "result of N=2, Alpha=1.0, a=1.0, No Error, Computation time: ", (time_NoError_N2_a1), '\n')
print(io, "result of N=2, Alpha=1.0, a=5.0, No Error, Computation time: ", (time_NoError_N2_a2), '\n')
print(io, "result of N=4, Alpha=1.0, a=5.0, No Error, Computation time: ", (time_NoError_N4_a2), '\n')


# figure with NO error bounds
function figure_System_NoError()

  Tmax = 10.0
  rr0 = 0.0
  solN4_a1 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=1.0)
  solN4_a2 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=5.0)
  solN4_a3 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=10.0)


  fig = plot(legend=:topright, xlab = L"\textrm{Time t}", ylab = L"\textrm{x(t) or y(t)} ", title="Numerial solution of the system of equations (with no error bounds)",
              legendfontsize=12,
              tickfont=font(10, "Times"),
              guidefontsize=10,
              xguidefont=font(10, "Times"),
              yguidefont=font(10, "Times"),
              bottom_margin=5mm,
              left_margin=5mm,
              right_margin=5mm,
              top_margin=5mm,
              size=(800, 600))
  
  plot!(fig, solN4_a1,  vars=(0, 1), color=:red, lc=:red, linewidth=2, label=L"x'=-x+axy \textrm{, a=1.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a2,  vars=(0, 1), color=:blue, lc=:blue, linewidth=2, label=L"x'=-x+axy \textrm{, a=5.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a3,  vars=(0, 1), color=:aquamarine, lc=:aquamarine, linewidth=2, linestyle=:dash, label=L"x'=-x+axy \textrm{, a=10.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a1,  vars=(0, 2), color=:green, lc=:green, linewidth=2, label=L"y'=-2y+2ay^{2} \textrm{, a=1.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a2,  vars=(0, 2), color=:orange, lc=:orange, linewidth=2, label=L"y'=-2y+2ay^{2} \textrm{, a=5.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a3,  vars=(0, 2), color=:darksalmon, lc=:darksalmon, linewidth=2, linestyle=:dash, label=L"y'=-2y+2ay^{2} \textrm{, a=10.0, N=4, } \alpha=1.0")

  return fig
              
end

fig = figure_System_NoError()
display(fig)
savefig(fig, joinpath(TARGET_FOLDER, "figure_1a_non_error.pdf"))


# figure with error bounds, no comparison
function figure_System_withError()

  Tmax = 10.0
  rr0 = 0.0
  solN4_a1 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=1.0)
  solN4_a1_bloat = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=true, resets=[4.0], alpha=1.0, a=1.0)

  fig = plot(legend=:topright, xlab = L"\textrm{Time t}", ylab = L"\textrm{x(t) or y(t)} ", title="Numerial solution of the system of equations (with error bounds, \n no comparison)",
              legendfontsize=12,
              tickfont=font(10, "Times"),
              guidefontsize=10,
              xguidefont=font(10, "Times"),
              yguidefont=font(10, "Times"),
              xlims= (1.0, 5.0),
              bottom_margin=5mm,
              left_margin=5mm,
              right_margin=5mm,
              top_margin=5mm,
              size=(800, 600))

  # carleman linearization solution with error bounds
  plot!(fig, solN4_a1_bloat,  vars=(0, 1), color=:green, lc=:green, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a1_bloat,  vars=(0, 2), color=:green, lc=:green, linewidth=2, linestyle=:dash, label=L"\textrm{error(y), a=1.0, N=4, } \alpha=1.0")

  # carleman linearization solution 
  plot!(fig, solN4_a1,  vars=(0, 1), color=:orange, lc=:orange, linewidth=2, label=L"x'=-x+axy \textrm{, a=1.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a1,  vars=(0, 2), color=:orange, lc=:orange, linewidth=2, label=L"y'=-2y+2ay^{2} \textrm{, a=1.0, N=4, } \alpha=1.0")

  return fig
end

fig = figure_System_withError()
display(fig)
savefig(fig, joinpath(TARGET_FOLDER, "figure_1b_error.pdf"))


# figure with error bounds with comparison
function figure_System_withError()

  Tmax = 10.0
  rr0 = 0.0
  solN4_a1 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=1.0)
  solN4_a1_bloat = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=true, resets=[4.0], alpha=1.0, a=1.0)

  solN4_a2 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=2.0)
  solN4_a2_bloat = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=true, resets=[4.0], alpha=1.0, a=2.0)

  fig = plot(legend=:topright, xlab = L"\textrm{Time t}", ylab = L"\textrm{x(t) or y(t)} ", title="Numerial solution of the system of equations (with error bounds, \n with comparison of different values of a)",
              legendfontsize=12,
              tickfont=font(10, "Times"),
              guidefontsize=10,
              xguidefont=font(10, "Times"),
              yguidefont=font(10, "Times"),
              xlims= (1.0, 6.0),
              ylims = (-0.02, 0.05),
              bottom_margin=5mm,
              left_margin=5mm,
              right_margin=5mm,
              top_margin=5mm,
              size=(800, 600))

  # carleman linearization solution with error bounds
  plot!(fig, solN4_a2_bloat,  vars=(0, 1), color=:lightgrey, lc=:lightgrey, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=2.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a2_bloat,  vars=(0, 2), color=:lightgrey, lc=:lightgrey, linewidth=2, linestyle=:dash, label=L"\textrm{error(y), a=2.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a1_bloat,  vars=(0, 1), color=:green, lc=:green, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a1_bloat,  vars=(0, 2), color=:green, lc=:green, linewidth=2, linestyle=:dash, label=L"\textrm{error(y), a=1.0, N=4, } \alpha=1.0")

  # carleman linearization solution 
  plot!(fig, solN4_a1,  vars=(0, 1), color=:orange, lc=:orange, linewidth=2, label=L"x'=-x+axy \textrm{, a=1.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a1,  vars=(0, 2), color=:orange, lc=:orange, linewidth=2, label=L"y'=-2y+2ay^{2} \textrm{, a=1.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a2,  vars=(0, 1), color=:red, lc=:red, linewidth=2, label=L"y'=-2y+2ay^{2} \textrm{, a=2.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a2,  vars=(0, 2), color=:red, lc=:red, linewidth=2, label=L"y'=-2y+2ay^{2} \textrm{, a=2.0, N=4, } \alpha=1.0")

  return fig
end

fig = figure_System_withError()
display(fig)
savefig(fig, joinpath(TARGET_FOLDER, "figure_1c_error.pdf"))


