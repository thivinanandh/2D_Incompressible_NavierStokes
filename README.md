[Github Link](https://github.com/thivinanandh/2D_Incompressible_NavierStokes_Solver)

# Incompressible Navier-Stokes Solver(2D) using Projection Scheme
---

## Compilation 
---

run the make file using 

```
make
```
command and then run the executable 

## Features
---
* Staggered Grid Solver for Stabilty
* Eigen Value based matrix generation for faster computation
* Pressure Projection scheme implementaton
* Iterative Solvers (BiCGSTAB)
* Output file format for tecplot visualisation


This is a 2D Incompressible Navier Stokes Solver written using C++. Here the code is configured for the example "Lid Driven Cavity" and the results are compared with the Ghia. et. al results. 


## Governing Equations
================


$$\textbf{u}_{\text{t}}  + (\textbf{u}\cdot\nabla) \textbf{u} + \nabla{P}  - \frac{1}{\text{Re}}( \nabla^{\text{2}} \textbf{u})    = \textbf{f}$$

$$\nabla \cdot \textbf{u}  = 0$$


## Methodology - Projection Scheme
===============================

In this Project , we have used a method called \"Projection Method\"
, which involves a *predictor and corrector* approach to finding the
velocity.At first , a predictor velocity $u^*$ is calculated by ignoring
the pressure gradient term in the Navier stokes equation . The corrector
velocity is then $u^{n+1}$ is calculated by adding the pressure
correction term in such a way that $u^{n+1}$ would have solved the
continuity equation .

$$\textbf{u}^{\text{*}}   = \textbf{u}^{\text{n}} - \Delta{t}(\textbf{u}\cdot\nabla \textbf{u})^{\text{n}} + \frac{\Delta{t}}{\text{Re}}( \nabla^{\text{2}} \cdot \textbf{u}^{\text{*}})$$

$$\textbf{u}^{\text{n+1}}   = \textbf{u}^{\text{*}} - \Delta{t} \;\nabla{P^{n+1}}$$


This is achieved by taking the divergence of the equation(4) which would
result in a continuity term and a Laplacian of pressure term. This
equation is commonly referred as the *Pressure Poisson* equation

$\nabla\cdot {\textbf{u}^{\text{n+1}}}   = \nabla\cdot {\textbf{u}^{\text{*}}} - \Delta{t} \;\nabla^{2}P$

Eliminate RHS since , 
$\nabla\cdot {\textbf{u}^{\text{n+1}}}  = 0$
asper continuity, then the equation becomes $\;\nabla^{2}P   = \frac{1}{\Delta{t}}\nabla\cdot {\textbf{u}^{\text{*}}}$. 
The Corrector equation to solve this system can be obtained from

$\textbf{u}^{\text{n+1}}   = \textbf{u}^{\text{*}} - \nabla{t} \;\nabla{P^{n+1}}$

### Boundary Conditions
-------------------


For lid driven cavity we will be using the following boundary
Conditions\
**Dirichlet type Boundary** for velocities at all the edges of the
system where the u on the top part will be the velocity of the lid and
all the other velocities are considered as zero\

$V_{lid}= 1.0 \;\;m/s$

**Zero Neumann Boundary** is provided for pressure
at all the edges of the system. This is given since the velocities at
the edges are zero, the pressure gradient which would be suitable at
edges to provide zero velocity at those edges should be zero.

## Implementation
==============

### Grid Structure
--------------

To solve this problem , we will be using a fully staggered grid , in
which the velocity and the pressure points will be staggered away from
the actual grid points which will make the computations more stable.The figure (1) \[He, Ping. (2016)\] of
staggered grid is given below

![ Control Volume for u](Images/Stag.png)

![ Control Volume for u](Images/Ugrid.jpg)

\

### Boundary Condition Implementation in Staggered Grid
---------------------------------------------------

The velocities are now calculated at staggered grid. This means at the
end of each computation the velocity and pressures at each grid point
will be computed by the interpolation of those quantities around the
actual grid points. So in similar way the boundary conditions are also
implemented in such a way that after the interpolation , the actual grid
points will have the desired Boundary conditions.

###  Velocity  

The u velocity has to be 1 at the top and 0 on all the other edges of
the system. Since only the interior grid points will be calculated in
the process,The cells at the corner are made negative of the cells
adjacent to them in certain directions (direction based on which the
interpolation is made for that particular velocity) so that after
interpolation they will produce zero at the actual grid points.

$$\begin{aligned}
u_{s [i,0]}  &= - u_{s [i,1]}  \forall i = 0,..,N \\
u_{s [i,N+1]},u_{s [i,N]} &=  V_{lid}  \forall  i = 0,..,N\end{aligned}$$

Similarly for v grid 

$$\begin{aligned}
v_{s [1,j]} \;\;\;\; &= \;\;\;\;- v_{s [0,j]}  \;\;\;\; \forall\;\;\;\; j = 0,..,N \\
v_{s [N+1,j]}\;\;\;\; &=\;\;\;\;  -v_{s [N,j]}  \;\;\;\; \forall\;\;\;\; j = 0,..,N\end{aligned}$$


### Pressure

The Pressure at staggered grids have to be chosen in such a way that
pressure difference at the edges have to be zero. This compliments the
zero velocity boundary condition that we have imposed at the edges.

$$\begin{aligned}
p_{s [i,0]}  &= \;\;\;\;p_{s [i,1]}  \;\;\;\; \forall\;\;\;\; i = 0,..,N+1 \\
p_{s [i,N+1]} &= \;\;\;\; p_{s [i,N]}  \;\;\; \forall\;\;\;\; i = 0,..,N+1\end{aligned}$$
$$\begin{aligned}
p_{s [0,j]} \;\;\;\; &= \;\;\;\; p_{s [1,j]}  \;\;\;\; \forall\;\;\;\; j = 0,..,N \\
p_{s [N+1,j]}\;\;\;\; &=\;\;\;\;  p_{s [N,j]}  \;\;\; \forall\;\;\;\; j = 0,..,N\end{aligned}$$

Discretizing Schemes 
---------------------

### Non linear term $\textbf{u}\cdot\nabla \textbf{u}$ 

The non linear term ,it can be written as
$$u\cdot\nabla \textbf{u} =  \frac{\partial u^2}{\partial x} + \frac{\partial uv}{\partial y}$$
The nonlinear term is solved explicitly at time time step n while
computing the predictor velocity in first one. This reduces the
complexity significantly by making the computations easier on the non
linear part.

The grid to be taken for solving the X momentum equation is shown in
figure below. This grid also physically signifies that the $U_{[i,j]}$
is caused due to the change of pressure in that direction.\

In this method , we will be using a **Central difference Scheme** for
solving the non linear part of the velocity . how ever for the
$\frac{\partial uv}{\partial y}$ part we do not have values at the exact
nodes. So in order to calculate them we will interpolate them to the
respective grid points and then perform the calculations. Formulating
the equations based on the scheme discussed above will result in the
following set of relations. Further *upwinding* has been incorporated in
order to avoid the instability that would be caused when the Re$>$2. The
relations used for calculating the non linear term are given below

$$\begin{aligned}
\frac{\partial u^2}{\partial x}  &= \frac{ u_{s [i+1,j]} - u_{s [i-1,j]}} {2h_x} \\
 \frac{\partial uv}{\partial y} &= \frac{\left(\frac{ u_{s [i,j+1]} - u_{s [i,j]}} {2 }\right)\left(\frac{ v_{s [i+1,j]} - u_{s [i,j]}} {2} \right) - \left(\frac{ u_{s [i-1,j]} - u_{s [i+1,j+1]}} {2 } \right)\left(\frac{ v_{s [i,j]} - u_{s [i-1,j]}} {2} \right)}{h_y}\end{aligned}$$

For $v$ , The control volume is chosen in similar way such that the
difference in pressure in y direction at grid boundaries drives the y
component of velocity ($v$) inside the control volume.

### Discretizing $\nabla\cdot \textbf{u}$

We have used **first order up-winding scheme** for calculating the
divergence of velocity while solving for pressure poisson equation (7).
The discretized form of relation for staggered grid is as given below

$$\begin{aligned}
     \nabla\cdot u^*  &= \left( \frac{ u_{s [i,j]} - u_{s [i-1,j]}} { h_x} + \frac{ v_{s [i,j]} - v_{s [i,j-1]}} { h_y} \right) \\.
 \end{aligned}$$
 
 
### Discretizing Pressure Terms

The pressure term in $\nabla^2{P}$ is formulated using the normal
**second order central difference scheme** for double derivative in the
pressure poisson equation. While solving for $\nabla{P}$ in the
corrector equation \[7\] , we have used **first order up-winding
scheme**.

$$\begin{aligned}
\nabla ^2 P  &= \left( \frac{ p_{s [i+1,j]} + p_{s [i-1,j]} - 2p_{s [i,j]}} { h_x ^2} + \frac{ p_{s [i,j+1]} + p_{s [i,j-1]} - 2p_{s [i,j]}} { h_y ^2} \right) \\.
\nabla P  &= \left( \frac{ p_{s [i+1,j]} - p_{s [i,j]} } { h_x } + \frac{ p_{s [i,j+1]} - p_{s [i,j]}} { h_y} \right) 
\end{aligned}$$


Solving Methodologies
=====================

 Predictor Step - Finding $u^*$
------------------------------

For this method , we arrange the equation in the way mentioned below and
we solve them as a **Helmholtz equation** . We can find the solution by
the direct method which involved finding the Eigen values and Eigen
vectors as a functions of sine and cosine functions for a second order
central difference matrix.

$$\begin{aligned}
    \left[ \nabla ^2 - \frac{\text{Re}}{\Delta{t}}h^2   \right] &= h^2\left[ \text{Re}\left( \frac{\partial u^2}{\partial x} + \frac{\partial uv}{\partial y} \right) - \frac{\text{Re}}{\Delta{t}}u^n   \right]
\end{aligned}$$

*Note : While solving the boundary term has to be included in the RHS of
the above equation*

 Solving Pressure Poisson Equation 
----------------------------------

Solving the pressure Poisson equation to find the pressure unknowns is
the Bottle neck in this Computation. The equation states

$$\nabla^{2}P^{n+1}   = \frac{1}{\Delta{t}}\nabla\cdot {\textbf{u}^{\text{*}}}$$


This can be solved as a Helmholtz function , if the Pressure values at
the boundary are already given. however in this problem we are
considering a Neumann Boundary problem where we only know the relation
between pressures at the edges and not the actual values.So we need to
edit alter the matrix as per the boundary conditions mentioned in
section 4.2 and use Iterative solvers to fnd the solution

In this Method , I **tried** using *Bi Conjugate gradient Stabilized
(BiCGSTAB)* method to solve for the pressure terms. Though BiCGSTAB has
higher convergence rate for these kind of matrices compared to other
iterative methods like Gauss Seidel or its own predecessor Conjugate
gradient method ,while Implementing it to solve for pressure poisson , I
could observe that the BiCGSTAB iterations were not converging during
some of the iterations for higher values of time steps or Reynolds
number.Further, the use of iterative solvers takes lot of time due to
the number of iterations it takes to converge compared to using the
direct methods like the Helmholtz function used in the previous section\
**Assumptions Made** : In order to overcome the above mentioned problem
, I made an assumption in the following problem that the pressure at the
end grid points are zero. This logical assumption is made on the fact
that

-   The pressure at end grid points will greatly influence only the
    velocity at the edge cells. And the pressure at interior cells are
    calculated comparatively accurate using the standard methodologies

-   The small magnitude of error produced in the solution due to this
    non zero pressure gradient at the edges can some how be controlled
    by adjusting the corrector step ,where the gradient of pressure is
    added as a corrector to the predicted velocity (7)

**Advantages** :

-   This can significantly reduces lot of iterations required to
    converge , since we are using a direct solver for solving the
    pressure Poisson step.\

**Disadvantages** :

-   As mentioned above , Due to the non zero pressure gradient induced
    at the edges , we need to control that error in to be propagated to
    velocity in the corrector equation by using a **very small time
    step** in the order of $10^{-4}$ to get convergence for higher
    Reynolds number

-   Because of using very small time steps the computation time may also
    increase , Since we have to more iterations to reach to steady state
    solution.



Results
=======

 Contour Plots of 'u' for various Reynolds Number
------------------------------------------------

All the contours obtained blow are computed in a 128\*128 uniform grid
structure. Further all the results are validated using the bench mark
results published by *Ghia Et Al. High-Re solutions for in compressible
flow using the Navier-Stokes equations and a multigrid method* . The
Main comparison is done for $u$ velocity in Y axis along the Mid
line(x=0.5)

###  Time evolution plots for Re = 100 

Time Step : 0.00005 sec\
Exit Tolerance : $10^{-8}$\
Time to solve - 300 minutes

![image](Images/Re_100_1.png)

![image](Images/Re_100_2.png)

![Evolution of $u$ velocity for Re 100](Images/Re_100_3.PNG)

![Evolution of $u$ velocity for Re 100](Images/Re_100_4.png)

###  Time evolution plots for Re = 1000 

Time Step : 0.00005 sec\
Exit Tolerance : $10^-{8}$\
Time to solve - 400 minutes

![image](Images/RE_1000_1.PNG)

![image](Images/RE_1000_2.PNG)

![Evolution of u velocity for Re 1000 ](Images/RE_1000_3.PNG)

![Evolution of u velocity for Re 1000 ](Images/RE_1000_LAST.PNG)

###  Time evolution plots for Re = 2000 

Time Step : 0.00005 sec\
Exit Tolerance : $10^-{8}$\
Time to solve - 440 minutes

![image](Images/RE_2000_1.PNG)

![image](Images/RE_2000_2.PNG)

![Evolution of u velocity for Re 2000 ](Images/RE_2000_3.PNG)

![Evolution of u velocity for Re 2000 ](Images/RE_2000_4.PNG)

###  Time evolution plots for Re = 5000 

Time Step : 0.00001 sec\
Exit Tolerance : $10^-{8}$\
Time to solve - 610 minutes

![image](Images/RE_5000_1.PNG)

![image](Images/RE_5000_2.PNG)

![Evolution of u velocity for Re 5000 ](Images/RE_5000_3.PNG)

![Evolution of u velocity for Re 5000 ](Images/RE_5000_4.PNG)

Observation
-----------

From the above figure ( Figure 4 ) we could clearly observe that the
method gas captured the primary vortex and all the other three secondary
vortices which could only be observed during a high Reynolds number flow

Comparison of Computed Solution vs Benchmark Solution
-----------------------------------------------------

![ Benchmark Comparison for Re 400 ](Images/RE-100.PNG)

![ Benchmark Comparison for Re 400 ](Images/RE-400.PNG)

![ Benchmark Comparison for Re 2000 ](Images/RE-1000.PNG)

![ Benchmark Comparison for Re 2000 ](Images/RE-2000.PNG)

![ Benchmark Comparison for Re 5000 ](Images/RE-3200.PNG)

![ Benchmark Comparison for Re 5000 ](Images/RE-5000.PNG)
