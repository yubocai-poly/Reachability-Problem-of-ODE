# Research Internship

## Quadratization preprocessing for the reachability Problem of ODE

### 1. Overview

--- 
One of the fundamental problems for dynamical systems is the reachability problem. The problem is: given a dynamical system (defined by a system of differential equations) and some range of the initial conditions, find a (rigorous and guaranteed) bound for the state of the system at some other time. This is an important ingredient, for example, of many approaches to verification of dynamical systems. The problem is very well studied for linear systems but is much less understood for the nonlinear case. A [[recent paper]](https://arxiv.org/pdf/2108.10390.pdf) proposes to reduce the nonlinear case to the linear one by using Carleman linearization. The approach presented in the paper requires the system to have at most quadratic nonlinearities. On the other hand, Gleb Pogudin, the supervisor of the project, and his colleagues have recently designed an algorithm for transforming any system into such an at-most-[[quadratic form]](https://arxiv.org/abs/2103.08013).

The goal of the project would be to consider specific systems with high-degree nonlinearities and apply the composition of the algorithms described above to solve the reachability problem. The main challenge is to adjust the transformed at most quadratic system in a way that would make the subsequent Carleman linearization as numerically stable as possible. Preliminary experiments show that such adjustment may have a tremendous impact on the final result. 

### 2. Useful Links

--- 
Papers:
- [Reachability of weakly nonlinear systems using Carleman linearization](https://arxiv.org/pdf/2108.10390.pdf)
- [Optimal monomial quadratization for ODE systems](https://arxiv.org/abs/2103.08013)

Code and Library:
- [Repeatibility evaluation for "Reachability of weakly nonlinear systems using Carleman linearization" (RP'21) by Marcelo Forets and Christian Schilling.](https://github.com/JuliaReach/RP21_RE)
- [Qbee](https://github.com/AndreyBychkov/QBee/)
