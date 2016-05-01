% Code Notes
% Jinming LYU
% 9 mars 2016

# Fileds List #
- `0 = "pos_ref"` : reference position of the drop (unused)
- `1 = "pos"` : deformed position of the drop
- `2 = "f_caps"` : surface forces
- `3 = "velocity"` : velocity of the drop surface
- `4 = "pos_lim"` : limit value of the position
- `5 = "f_caps_lim"` : limit value of the surface forces
- `6 = "velocity_lim"` : limit value of the velocity
- `7 = "normal_load"` : value of the normal component of the drop forces.
- `8 = tension_1`
- `9 = tension_2`
- `10 = direction_1`
- `11 = direction_2`

Note : the fields 4, 5, 6, which correspond to limit values of fields 1, 2, 3 are never used in the computation, just in the postprocessing.

#程序总体#
对于所有的单个Soft Object，主要的计算都由`SimulateOneObjectInFlow.java`中的*sc.softobject.computeOneTimeStep*调用
（也是不同Soft Object分别运行的开始）。该函数分为RK45格式和梯形积分公式两种。

##Drop##
对于Drop，时间推进公式的主要调用是`computeSurfaceVelocity`方法（位于Drop.java中）。该方法
Solves for velocity at given position, taking into account surface tension and surface viscosity, as necesssary.
具体而言该方法主要有如下两个方面的调用：`computeMembraneForcesByFEM.computeDropForcesConstantTension(m, 2, MassMatrix);`（位于
interfacialForcesComputation内）和`computeVelocityAndStoreInField`（位于velocityComputation中的velocityComputationByBEM内)。
computeVelocityAndStoreInField通过调用computeVelocityAndStoreInFieldWithSingularitySubstractionOPTIAndMT_CPU计算格林函数相关。



# softObjects #
## Drop.java ##
`public class Drop extends SoftObject implements Serializable`

- `public void computeSurfaceVelocity(Mesh m, double[] gvel)`
   Solves for velocity at given position, taking into account surface tension and surface viscosity, as necesssary.

- `public double[] computeOneTimeStepRK45(Mesh m1, double time, double dt)`
  Computes one time step with a Runge-Kutta-Fehlberg method (adaptative time step): using 6 intermediates computation, 
  form two combinations of them, one which is 4 th order accurate, the other one 5 th order. 
  Then, compares the two to estimate the error and updates the time step such that the error is always 
  between $[\epsilon_{min},\epsilon_{max}]$. If the error > $\epsilon_{max}$, then reduces the time step; and if 
  error < $\epsilon_{min}$, then increases the time step. The parameters of RKF45 can be see in [Ref](./ref/RKF45.pdf)

- `public void computeSurfaceVelocity(Mesh m, double[] gvel)`: Solves for velocity at given position, taking into account 
  surface tension and surface viscosity, as necesssary.  Velocity is collocated inside the routine.

# mesh #
## Mesh.java ##
`public class Mesh implements Serializable`

- `public void computeCollocationMatrixNEW()`: Compute the collocation matrix used to convert value at nodes to nodal values.
  Stores it in the mesh instead of in the simulation context (needed for suspensions).

  Inside the `computeCollocationMatrixNEW()` in Mesh.java, it has called a method `computeShapeFunctionAt(double[])` (In package mesh). 
  The nodes distributed in the standard shape are: 

`  2------o-------o   `

  |$\ \ \ \ \ \ \ \ \ \ \ $|$\ \ \ \ \ \ \ \ \ \ \ \ $|   
  |$\ \ \ \ \ \ \ \ \ \ \ $|$\ \ \ \ \ \ \ \ \ \ \ \ $|   
`  5------4-------o   `

  |$\ \ \ \ \ \ \ \ \ \ \ $|$\ \ \ \ \ \ \ \ \ \ \ \ $|   
  |$\ \ \ \ \ \ \ \ \ \ \ $|$\ \ \ \ \ \ \ \ \ \ \ \ $|   
`  0------3-------1   `

 \begin{equation*}
   N_0(\xi,\eta)=-\lambda (1-2\lambda),\ N_1(\xi,\eta)=-\xi (1-2\xi)
 \end{equation*}
 \begin{equation*}
   N_2(\xi,\eta)=-\eta (1-2\eta),\ N_3(\xi,\eta)=4\xi\lambda
 \end{equation*}
 \begin{equation*}
   N_4(\xi,\eta)=4\xi\eta,\ N_5(\xi,\eta)=4\eta\lambda
 \end{equation*}
 where $\lambda=1-\xi-\eta$.

## RemeshingByMembraneForces.java ##
Implementation of a remeshing algorithm based on membrane (i.e. elastic) forces.
The algorithm is as follows :

1. the mesh nodes are moved along the surface by elastic forces (no normal displacement)

2. the intensity of the remeshing can be controlled by a parameter which measures an "elastic" energy of the mesh. 
   The closest to zero the parameter is, the strongest the remeshing.

3. the forces are calculated by finite element computations of the elastic forces due to the deformation of the mesh
   between the current state and a reference position.

4. the reference position is arbitrary : the idea is to use some state where the mesh was correct. At the moment,
   the idea is to use the position at some previous time step as reference position. If the initial mesh was correct,
   and the remeshing is done frequently enough, this should preserve a correct mesh. However, it does not deals with
   possible needs for local mesh "refinements", i.e. nodes should concentrate in more complex regions while flat/simple
   regions could be described with less nodes. At the moment, the idea is to couple this remeshing with a "nodes concentrating" 
   algorithm based e.g. on local curvature.


#simulation#
##SimulateOneObjectInFlow.java##
1. The `run(Mesh m1)` method in `public class SimulateOneObjectInFlow`
    - `if (sc.timeStepAlgorithm.equals("stat")||sc.timeStepAlgorithm.equals("hybrid"))`, the simulation is advanced 
      with a equal time step;
    - `else`: we set a *Initial time step* `dt` and a *Initial time* `time`, `computeOneTimeStep(m1, time, dt)` will update
      the solution and returen the new `time` and `dt`.
      For example, `computeOneTimeStepRK45` will compute one time step with a Runge-Kutta-Fehlberg method (adaptative time step):
      using 6 intermediates computation, form two combinations of them, one which is
      4 th order accurate, the other one 5 th order. Then, compares the two to estimate the error
      and updates the time step such that the error is always between $[\epsilon_{min},\epsilon_{max}]$.
        - If the error is > $\epsilon_{max}$, then reduces the time step
        - If the error is < $\epsilon_{max}$, then increases the time step


#field#
`Field` and `Node` are two abstract class in this package, `node=vetex+field`. The `vertex` defines the mesh (its coordinates,
 the elements related to it, the edges connected to it and the nodes defined on it). `Field` has two subclass `ScalarField` 
 and `VectorField`.

##Field.java##
`public abstract class Field implements Serializable`, Abstract class representing a physical field (velocity, forces, tension) onto a mesh.

