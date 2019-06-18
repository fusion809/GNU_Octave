# GNU Octave scripts
This repository contains my GNU Octave scripts. 

| Script           | Description                                                                  |
|------------------|------------------------------------------------------------------------------|
| airy-newton.m    | Determines the first root, on the negative x axis, of the Airy Ai(x) function.|
| am-calc-aix-prob-interp.m | Determines the coefficients in a Chebyshev series for a function f(x) using inner products. Root mean square errors are calculated too. |
| bash.m           | Starts a bash shell with the command 'bash' |
| definite-integration.m | Definite integration using lsode compared to using a Chebyshev spectral method. Usually lsode outperforms the spectral method. |
| first.m          | Solves Schrodinger equation with V = k x on [0,200] (which is meant to approximate infinity) with N = 1000. |
| jerk.m           | Jerk equation (from [Wikipedia](https://en.wikipedia.org/wiki/Chaos_theory#Jerk_systems)) function script. |
| jerk-ex.m        | Solves the aforementioned jerk equation and plots it.                        |
| lorenz-ex.m      | Solves the Lorenz equations.                                                 |
| lorenz.m         | Lorenz equation function script.                                             |
| lorenz-plot.m    | Plots lorenz-ex.m solution.                                                  |
| ode-xy-cos-y2.m  | Solves ODE specified by f.m                                                  |
| push.m           | A failed attempt at making a function that commits git changes with commit message input, x. |
| prayer.m         | Redundant script I used to determine how close I was to 99 prayer on RuneScape  |
| randex.m         | Contains randex function, which is the diff function for a random ODE.       |
| RandomExample2.m | Solves d2y/dx2 = (C^2)/(y^2) - lambda^2 * (y(1))^(n-1) over [0,10] with y(0)=1; dy/dx(0) = 1 |
| simpen.m         | Differential function for the simple pendulum.                               |
| simpen-ex.m      | Solves the problem of the simple pendulum.                                   |
| SLEq.m           | Solves 1D Schrodinger equation with linear potential.                        |
| vim.m            | Starts vim editor with bizarre rendering issues. |
| zsh.m            | Starts a zsh shell, but it has some weird output created on every new line. |
