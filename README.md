# Personal Research Internship

## Quadratization in the reachability problem for ODEs

### 1. Overview

--- 
I would like to thank the kind professor [Gleb Pogudin](http://www.lix.polytechnique.fr/Labo/Gleb.POGUDIN/) for hosting me during this project in the [MAX](http://www.lix.polytechnique.fr/max/max-web/max/max-home.en.html) team at the [Laboratoire d’informatique de l’École Polytechnique (LIX)](https://www.lix.polytechnique.fr/) and [Centre national de la recherche scientifique (CNRS)](https://www.cnrs.fr/) the Spring of 2023. I am grateful for the guidance from my supervisor Prof. Gleb Pogudin on my research direction and the ideas I received when I encountered a bottleneck in my research.

One of the fundamental problems for dynamical systems is the reachability problem. The problem is: given a dynamical system (defined by a system of differential equations) and some range of the initial conditions, find a (rigorous and guaranteed) bound for the state of the system at some other time. This is an important ingredient, for example, of many approaches to the verification of dynamical systems. The problem is very well studied for linear systems but is much less understood for the nonlinear case. Forets and Schilling [4] propose to reduce the nonlinear case to the linear one by using Carleman linearization. The approach presented in the paper requires the system to have at most quadratic nonlinearities. On the other hand, Gleb Pogudin, the supervisor of the project, and his colleagues have recently designed an algorithm for transforming any system into such an at most quadratic form [1].

This project considers specific systems with high-degree nonlinearities and applies the composition of the algorithms described above to solve the reachability problem. The main challenge is to adjust the transformation at most quadratic systems in a way that would make the subsequent Carleman linearization as numerically stable as possible. In the project, we try to first explore 3 typical ODE systems with different characteristics, and then we try to prove the existence of quadratic method generation in a $2 \times n$ (where $2$ is the number of variables and $n$ is the highest dimension) ODE system with linear part of negative eigenvalue to a quadratic form with a linear part of negative eigenvalue and expand this conclusion into general cases.

**Note: In this project, we mainly use Python and Julia to perform numerical operations.**


### 2. Useful Links

--- 
Papers:
- [Reachability of weakly nonlinear systems using Carleman linearization](https://arxiv.org/pdf/2108.10390.pdf)
- [Optimal monomial quadratization for ODE systems](https://arxiv.org/abs/2103.08013)

Code and Library:
- [Repeatibility evaluation for "Reachability of weakly nonlinear systems using Carleman linearization" (RP'21) by Marcelo Forets and Christian Schilling.](https://github.com/JuliaReach/RP21_RE)
- [Qbee](https://github.com/AndreyBychkov/QBee/)

### 3. Bibliography

---

[1] Andrey Bychkov and Gleb Pogudin. Optimal monomial quadratization for ODE systems. 2021. arXiv: 2103.08013 [cs.SC].

[2] D. C. Carothers et al. “SOME PROPERTIES OF SOLUTIONS TO POLYNOMIAL SYSTEMS OF DIFFERENTIAL EQUATIONS”. In: 2005.

[3] John E Dennis Jr and Robert B Schnabel. Numerical methods for unconstrained optimization and nonlinear equations. SIAM, 1996.

[4] Marcelo Forets and Christian Schilling. “Reachability of Weakly Nonlinear Systems Using Carleman Linearization”. In: Lecture Notes in Computer Science. Springer International Publishing, 2021, pp. 85–99. DOI: 10.1007/978-3-030-89716-1_6. URL: https://doi.org/10.1007%2F978-3-030-89716-1_6.

[5] Roger A Horn and Charles R Johnson. Matrix analysis. Cambridge university press, 2012.

[6] Adolf Hurwitz et al. “On the conditions under which an equation has only roots with negative real parts”. In: Selected papers on mathematical trends in control theory 65 (1964), pp. 273–284.

[7] Ivana Kovacic and Michael J Brennan. The Duﬀing equation: nonlinear oscillators and their behaviour. John Wiley & Sons, 2011.
