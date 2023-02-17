## Week 0 Report - Familiarization with the QBee package and quadratization of the system


### 1. Useful informations

- **QBee package**: [[link]](https://github.com/AndreyBychkov/QBee/)
- **QBee documentation**: [[link]](https://qbee.readthedocs.io/)
- **Examples**: [[link]](https://github.com/AndreyBychkov/QBee/tree/master/examples)
- **Optimal monomial quadratization for ODE systems**: [[link]](https://arxiv.org/abs/2103.08013)

---

### 2. Introduction to Quadratization Problem
For a nonlinear ODE system, on the right hand, we have varaible with degree $\geq 2$. For example, the following system is nonlinear:

$$
x' = -x + x^3 + y^2
$$

However, we find that for this system, there exists $x^3$. Quadratization is to introduce new variables in order to decrease the degree of the system into at most 2. Notice that $xy$ is in degree 2 but $xy^2$ is in degree 3.

### 3. Quadratization with QBee
Here we use the QBee package to quadratize the system, for example, the following system:
$$
x'=-x+ax^{3}
$$


```python
a = parameters("a")
x = functions("x")

eq = - x + a * x ** 3
system = [(x, eq)]

quadr_system = polynomialize_and_quadratize(system, 
                                            new_vars_name = "w")

==================================================
Quadratization result
==================================================
Number of introduced variables: 1
Nodes traversed: 3
Introduced variables:
w{0} = x**2                                         
```
From the result, we can introduce a new variable $w_{0}=y=x^2$ to quadratize the system. Then for the original equation, we have:
$$
\begin{aligned}
x' &= x'=-x+ax^{3} \\
&= -x +axy
\end{aligned}
$$
However, since we also introduce a new variable $w_{0}=y=x^2$, we can also need to compute the derivative of $w_{0}$ in order to find the ODE of it.
$$
\begin{aligned}
y' &= 2x x' = 2x(-x+ax^{3}) \\
&= -2x^{2}+2ax^{4}\\
&= -2y +2ay^2
\end{aligned}
$$
Then we can combine the two equations to get the quadratized system:
$$
\begin{cases}
y' = -2y+2ay^{2}\\
x' = -x + axy
\end{cases}
$$
Here, on the right hand side, we can see that all the variables are in degree at most 2.

---

### 4. Other Examples
For more examples and detailed code, please refer to the [examples](https://github.com/yubocai-poly/Reachability-Problem-of-ODE/blob/main/quadra_code/quadratization.ipynb). Here I only show the quadratization of the following system:


We have the third system of following form:

$$
x'' = kx+ax^{3}
$$
(we set $x_{0}=x$ and $x_{1}=x'$). Then we have the system of the ODEs:
$$
\begin{cases}
x_{0}'=x_1\\
x_{1}'=kx_{0}+ax_{0}^{3}
\end{cases}
$$


```python
k, a = parameters("k, a")
x0, x1 = functions("x0, x1")

eq1 = -k * x0 - a * x0**3
system = [
    (x1, eq1),
    (x0, x1),
]

quadr_system = polynomialize_and_quadratize(system, input_der_orders={x0: 2})

==================================================
Quadratization result
==================================================
Number of introduced variables: 4
Nodes traversed: 83
Introduced variables:
w_{0} = x0**3
w_{1} = x0**3 # Actually, I don't know why it is the same as w_{0}
w_{2} = x0*x0'
w_{3} = x0**2
```

From the program we introduce 3 variables:
$$
\begin{cases}
w_{0} = x_{0}^3\\
w_{1} = x_{0} x_{0}'\\
w_{2} = x_{0}^2
\end{cases}
$$
Then we got the following:
$$
\begin{cases}
x_{0}'=x_1\\
x_{1}'=kx_{0}+aw_{0}
\end{cases}
$$
Now we need to dealing with $w_{0}$, $w_{1}$ and $w_{2}$, we have the following computations:
$$
\begin{aligned}
w_{0}' &= 3x_{0}^2 x_{0}' \\
&= 3 w_{2} x_{1} \\
\end{aligned}
$$
$$
\begin{aligned}
w_{1}' &= x_{0} x_{0}'' + x_{0}' x_{0}' \\
&= x_{1}^2 + x_{0} x_{1}' = x_{1}^2 + x_{0} (kx_{0}+ax_{0}^{3})\\
& = x_{1}^2 + kx_{0}^{2}+ak x_{0}^{4}\\
& = x_{1}^2 + k w_{2} + ak w_{2}^2\\
\end{aligned}
$$
$$
w_{2}' = 2x_{0} x_{0}' = 2 w_{1}
$$
Therefore, we have the following system:
$$
\begin{cases}
x_{0}'=x_1\\
x_{1}'=kx_{0}+aw_{0}\\
w_{0}'=3 w_{2} x_{1}\\
w_{1}'=x_{1}^2 + k w_{2} + ak w_{2}^2\\
w_{2}'=2 w_{1}
\end{cases}
$$