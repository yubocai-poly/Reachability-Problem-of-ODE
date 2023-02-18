## Week 1 Report - Apply Reachability Analysis to our first Example System

### 1. Useful informations

- **Reachability package**: [[link]](https://github.com/JuliaReach)
- **Carleman linearization**: [[link]](https://github.com/JuliaReach/RP21_RE)
- **Reachability of weakly nonlinear systems using Carleman linearization**: [[link]](https://arxiv.org/abs/2108.10390)

---

### 2. Introduction to Reachability Analysis of Carleman Linearization for quadratized system
Here we use the [[Reachability package]](https://github.com/JuliaReach) which is in julia and the system code [here](https://github.com/JuliaReach/RP21_RE). Here we bring our first example system to the reachability analysis. The system is:

$$
x'=-x+ax^{3}
$$

From quadratization, we know that we can introduce a new variable $w_{0}=y=x^2$ to quadratize the system. Then for the original equation, we have:

$$
\begin{cases}
x' = -x + axy\\
y' = -2y+2ay^{2}
\end{cases}
$$

Then we use the Carleman linearization to linearize the system. Follow the algorithm in the paper [[Reachability of weakly nonlinear systems using Carleman linearization]](https://arxiv.org/abs/2108.10390), we have the $F_1 \in \mathbb{R}^{n \times n}$ matrix for linear part:

$$
F_1=
\begin{matrix}
x \\ y
\end{matrix}
\begin{bmatrix}
-1 & 0 \\
0 & -2  
\end{bmatrix}.
$$

And for the nonlinear part we have $F_2 \in \mathbb{R}^{n \times n^{2}}$ matrix:

$$
F_2=
\begin{matrix}
x \\ y
\end{matrix}
\begin{bmatrix}
0 & 2a & 0 & 0\\
0 & 0 & 0 & a
\end{bmatrix}.
$$

where we have the basis for the nonlinear part is $[x^{2}, xy, yx, y^{2}]$. Since we use quadratization algorithsm to quadratize the system, we should at most have this 4 basis.

Therefore, we can use the algorithsm from the paper **Epidemic model(SEIR)** to analyse this emaple system. The code is [here]().