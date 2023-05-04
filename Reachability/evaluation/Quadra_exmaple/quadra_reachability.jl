# =================
# Dependencies
# =================

using ReachabilityAnalysis, CarlemanLinearization
using Plots, LaTeXStrings, LinearAlgebra, SparseArrays
using Plots.PlotMeasures
using LazySets: center
using CarlemanLinearization: _error_bound_specabs_R

include("../utils.jl")
include("quadra_ODE_example1.jl")

# Ploting the results
Tmax = 10.0
rr0 = 0.0

solN2_minimum_bloat = _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=2*sqrt(2), a=0.1, x0=0.5)
solN2_alpha2_bloat = _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=2.0, a=0.1, x0=0.5)
solN2_alpha1_bloat = _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=1.0, a=0.1, x0=0.5)

function figure_System_withError_minimum()

  fig = plot(legend=:bottomright, xlab = L"\textrm{Time t}", ylab = L"\textrm{x(t)} ", title="Numerial solution of the system of equations of x(t) \n (with error bounds, N=2, x0=0.5)",
  legendfontsize=12,
  tickfont=font(10, "Times"),
  guidefontsize=10,
  xguidefont=font(10, "Times"),
  yguidefont=font(10, "Times"),
  xlims = (2.0, 6.0),
  ylims = (-0.02, 0.02),
  bottom_margin=5mm,
  left_margin=5mm,
  right_margin=5mm,
  top_margin=5mm,
  size=(800, 600))

  # carleman linearization solution with error bounds
  plot!(fig, solN2_alpha1_bloat,  vars=(0, 1), color=:lightblue, lc=:lightblue, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=2, } x_{0}=0.5, \alpha=1.0")
  plot!(fig, solN2_alpha2_bloat,  vars=(0, 1), color=:grey, lc=:grey, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=2, } x_{0}=0.5, \alpha=2.0")
  plot!(fig, solN2_minimum_bloat,  vars=(0, 1), color=:green, lc=:green, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=2, } x_{0}=0.5, \alpha=\frac{\sqrt{N}}{x_0}=2\sqrt{2}(minimum)")


  # carleman linearization solution 
  # plot!(fig, solN4_a1,  vars=(0, 1), color=:orange, lc=:orange, linewidth=2, label=L"x'=-x+axy \textrm{, a=1.0, N=4, } \alpha=1.0")
  # plot!(fig, solN4_a2,  vars=(0, 1), color=:red, lc=:red, linewidth=2, label=L"y'=-2y+2ay^{2} \textrm{, a=2.0, N=4, } \alpha=1.0")
  # plot!(fig, solN4_a1,  vars=(0, 1), color=:blue, lc=:blue, linewidth=2, label=L"x'=-x+axy \textrm{, a=3.0, N=4, } \alpha=1.0")

return fig
end


fig = figure_System_withError_minimum()
display(fig)
savefig(fig, joinpath(TARGET_FOLDER, "figure_minimum_error1.pdf"))

# Try with different value in N
solN3_minimum_bloat = _solve_system_carlin(N=3, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=2*sqrt(3), a=0.1, x0=0.5)
solN3_alpha2_bloat = _solve_system_carlin(N=3, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=10.0, a=0.1, x0=0.5)
solN3_alpha1_bloat = _solve_system_carlin(N=3, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=1.0, a=0.1, x0=0.5)

function figure_System_withError_minimum2()

  fig = plot(legend=:topright, xlab = L"\textrm{Time t}", ylab = L"\textrm{x(t)} ", title="Numerial solution of the system of equations of x(t) \n (with error bounds, N=3, x0=0.5)",
  legendfontsize=12,
  tickfont=font(10, "Times"),
  guidefontsize=10,
  xguidefont=font(10, "Times"),
  yguidefont=font(10, "Times"),
  xlims = (3.0, 6.0),
  ylims = (0, 0.01),
  bottom_margin=5mm,
  left_margin=5mm,
  right_margin=5mm,
  top_margin=5mm,
  size=(800, 600))

  # carleman linearization solution with error bounds
  plot!(fig, solN3_alpha1_bloat,  vars=(0, 1), color=:lightblue, lc=:lightblue, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=3, } x_{0}=0.5, \alpha=1.0")
  plot!(fig, solN3_alpha2_bloat,  vars=(0, 1), color=:grey, lc=:grey, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=3, } x_{0}=0.5, \alpha=10.0")
  plot!(fig, solN3_minimum_bloat,  vars=(0, 1), color=:green, lc=:green, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=3, } x_{0}=0.5, \alpha=\frac{\sqrt{N}}{x_0}=2\sqrt{3}(minimum)")


return fig
end

fig = figure_System_withError_minimum2()
display(fig)
savefig(fig, joinpath(TARGET_FOLDER, "figure_minimum_error2.pdf"))

# Try with different value in x0
solN2_minimum_bloat_x0 = _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=4*sqrt(2), a=0.1, x0=0.25)
solN2_alpha2_bloat_x0 = _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=10.0, a=0.1, x0=0.25)
solN2_alpha1_bloat_x0 = _solve_system_carlin(N=2, T=Tmax, δ=0.01, radius0=rr0, bloat=true, alpha=1.0, a=0.1, x0=0.25)

function figure_System_withError_minimum3()

  fig = plot(legend=:topright, xlab = L"\textrm{Time t}", ylab = L"\textrm{x(t)} ", title="Numerial solution of the system of equations of x(t)\n  (with error bounds, N=2, x0=0.25)",
  legendfontsize=12,
  tickfont=font(10, "Times"),
  guidefontsize=10,
  xguidefont=font(10, "Times"),
  yguidefont=font(10, "Times"),
  xlims = (3.0, 6.0),
  ylims = (0, 0.01),
  bottom_margin=5mm,
  left_margin=5mm,
  right_margin=5mm,
  top_margin=5mm,
  size=(800, 600))

  # carleman linearization solution with error bounds
  plot!(fig, solN2_alpha1_bloat_x0,  vars=(0, 1), color=:lightblue, lc=:lightblue, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=2, } x_{0}=0.25, \alpha=1.0")
  plot!(fig, solN2_alpha2_bloat_x0,  vars=(0, 1), color=:grey, lc=:grey, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=10.0, N=2, } x_{0}=0.25, \alpha=10.0")
  plot!(fig, solN2_minimum_bloat_x0,  vars=(0, 1), color=:green, lc=:green, linewidth=2, linestyle=:dash, label=L"\textrm{error(x), a=1.0, N=2, } x_{0}=0.25, \alpha=\frac{\sqrt{N}}{x_0}=4\sqrt{2}(minimum)")


return fig
end

fig = figure_System_withError_minimum3()
display(fig)
savefig(fig, joinpath(TARGET_FOLDER, "figure_minimum_error3.pdf"))