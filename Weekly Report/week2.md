## Week 2 Report - Optimal Quadratization ($\alpha$ value) for Reachability Analysis of the first Example System

### 1. Introduction to the Problem

From previous week, we study the ODE $x'=-x+ax^{3}$ with the following quadratization:

$$
\begin{cases}
x' = -x + axy\\
y' = -2y+2ay^{2}
\end{cases}
$$

where we introduced a new variable $w_{0}=y=x^{2}$ to quadratize the system. Then, we studied the reachability analysis of the system using the Carleman linearization. 

Now we want to take a step further and study the optimal quadratization for the system which also produce the least error bound in the reachability analysis. For optimal quadratization, we have several options: $w_{0}=y=x^{2}$ or $w_{0}=y=\alpha x^{2}$ which are equivalent. If we take $w_{0}=y=\alpha x^{2}$, then we have the following computation:

$$
\begin{aligned}
y' &= 2\alpha x x' \\
&= 2\alpha x (-x + ax^3) \\
&= -2\alpha x^{2} + 2\alpha ax^4 \\
\end{aligned}
$$

Since we have $y=\alpha x^2 \Rightarrow x^2 = \frac{1}{\alpha}y$, therefore we have:

$$
-2\alpha x^{2} + 2\alpha ax^4 = -2y + \frac{2a}{\alpha} y^2
$$

Then we have the following system:

$$
\begin{cases}
x' = -x + \frac{a}{\alpha}xy\\
y' = -2y + \frac{2a}{\alpha} y^2
\end{cases}
$$

If we denote the new parameter $b = \frac{a}{\alpha}$, then the system goes back to the original format. Therefore, we can only study the reachability analysis with $b$.

$$
\begin{cases}
x' = -x + =bxy\\
y' = -2y+2by^{2}
\end{cases}
$$

---

### 2. Computing the Error Bound and Reachability

From previous analysis,  we have the $F_1 \in \mathbb{R}^{n \times n}$ matrix for linear part:

$$
F_1=
\begin{bmatrix}
-1 & 0 \\
0 & -2  
\end{bmatrix}.
$$

And for the nonlinear part we have $F_2 \in \mathbb{R}^{n \times n^{2}}$ matrix:

$$
F_2=
\begin{bmatrix}
0 & 2b & 0 & 0\\
0 & 0 & 0 & b
\end{bmatrix}.
$$

For matrix $F_1$ we have the eigenvalues $\lambda_1 = -1$ and $\lambda_2 = -2$. Then we got $\Re\left(\lambda_1\right)=-1$ (the real part of $\lambda_1$). We cite the [paper](https://arxiv.org/pdf/2108.10390.pdf) for reachabiliy analysis, we have following part:

**Definition 1.** System is said to be weakly nonlinear if the ratio

$$
R:=\frac{\left\|x_0\right\|\left\|F_2\right\|}{\left|\Re\left(\lambda_1\right)\right|}
$$

satisfies $R<1$.

**Definition 2.** System (1) is said to be dissipative if $\Re\left(\lambda_1\right)<0$ (i.e., the real part of all eigenvalues is negative).

$$
\text{The conditions } \Re\left(\lambda_1\right)<0 \text{ and } R<1 \text{ ensure arbitrary-time convergence.}
$$

**Theorem 1 ([30, Corollary 1])**. Assuming that (1) is weakly nonlinear and dissipative, the error bound associated with the linearized problem (2) truncated at order $N$ satisfies

$$
\left\|\eta_1(t)\right\| \leq \varepsilon(t):=\left\|x_0\right\| R^N\left(1-e^{\Re\left(\lambda_1\right) t}\right)^N,
$$

with $R$ as defined in (5). This error bound holds for all $t \geq 0$.

Then we compute the error bound for our example system:

$$
\left\|F_2\right\| = sup(F_2) = 2b
$$

Then we have for $R$:

$$
R:=\frac{\left\|x_0\right\|\left\|F_2\right\|}{\left|\Re\left(\lambda_1\right)\right|} = 2b \left\|x_0\right\| 
$$

Therefore, we have the formula for error bound:

$$
\begin{aligned}
\left\|\eta_1(t)\right\| \leq \varepsilon(t)& :=\left\|x_0\right\| R^N\left(1-e^{\Re\left(\lambda_1\right) t}\right)^N \\
& = \left\|x_0\right\|^{N+1} 2b^{N} \left(1-e^{-t}\right)^{N}
\end{aligned}
$$

$$
\lVert x \rVert 
$$