# =================
# Dependencies
# =================

using ReachabilityAnalysis, CarlemanLinearization
using Plots, LaTeXStrings, LinearAlgebra, SparseArrays
using Plots.PlotMeasures
using LazySets: center
using CarlemanLinearization: _error_bound_specabs_R

include("../utils.jl")
include("quadra_ODE.jl")

# Ploting the results
Tmax = 10.0
rr0 = 0.0
solN4_a1 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=1.0)
solN4_a5 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=10.0)
solN4_a10 = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=false, alpha=1.0, a=20.0)

# figure with NO error bounds, plot for x
function figure_System_NoError()
  
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
  plot!(fig, solN4_a5,  vars=(0, 1), color=:blue, lc=:blue, linewidth=2, label=L"x'=-x+axy \textrm{, a=5.0, N=4, } \alpha=1.0")
  plot!(fig, solN4_a10,  vars=(0, 1), color=:aquamarine, lc=:aquamarine, linewidth=2, linestyle=:dash, label=L"x'=-x+axy \textrm{, a=10.0, N=4, } \alpha=1.0")
  
  return fig
              
end