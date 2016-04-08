% Des Questions
% Jinming LYU
% \today

## Le 8 mars 2016 ##
- [PROBLEME PHYSIQUE, P.4]公式$F_{\gamma}=\int_S\gamma dS$是如何将不可压条件($\nabla\cdot\mathbf{V}=0$)考虑进来的？推导过程。

## Le 9 mars 2016 ##
- Under mathsTools package, the class `SparseMatrixColumnMajor` and `SparseMatrixRowMajor` are both extented from 
abstract class `SparseMatrixStorage`. And each cloumn (or row) is TreeMap object, but what is the difference of their
final object? BTW, I think there are a little error about the index in one of the constructor.

## Le 10 mars 2016 ##
- Is all the quantities are dimensionless numbers?
    - `OUI`. The code solve the dimensionless equations.
- As for the *Clean droplet in shear flow* in JFM 2016, when $Ca=1E-3$ (namely the flow strength) and $\gamma=1$, why
  the computation time decrese as I increse $\lambda$? (I've set the maximum compuational time as 50)

## Le 22 mars 2016 ##
- What is the physical meaning of dipole, quadruple, sextuple and octuple? (Sec. 2.2.7 Multipoles of Green’s functions)
  The Physical meaning in aeroacoustic can be found [here](./refQuestion/空气动力性噪声.pdf)

## Le 3 avril 2016 ##
- In *velocityComputationByBEM.java*, line 1858, what dose 

`Elements e = m.getFields().get(outFieldVelocity).getNodes()[i].getVertex().getElements().get(1);` 

really means?

## Le 5 avril 2016 ##
- The functioning of `setOneRing()` method in `LoopElement.java`?
- Why the method `computeGeometricParameters(Mesh m1, PrintWriter file, double time)` in `SimulateOneObjectInFlow.java` 
  has run `i` and `k` loops.
  And also, what dose `get(1)` means in `shape[l]*m1.getElements()[i].getOneRing()[l].getNodes().get(1).getDofValueAt(j)` in
  line 928?

## Le 7 avril 2016 ##
- How is the fourth-order derivatiion is computed? By $f^{(4)}=\sum N_i \alpha_i$ or computed with numerical methods, like FD.
- The direct application of Differential Geometry is in the class `computeMembraneForcesByFEM`? And this class use membrane law
  of soild mechanics to compute the force?
- What dose field 8-11 mean in the method `computeViscoelasticStress` of class `computeMembraneForcesByFEM`?
- In method `computeViscoelasticStress` of class `computeMembraneForcesByFEM`, is $x_0$ the cell's inner integration point?
  The formulation for computing softobject volume and center? And also for the inertia matrix?
- For the method `computeMembraneForcesViscousTension` in class `computeMembraneForcesByFEM`, why `inFieldVelocity`
  is used to as an input variable?
- The method `convertValueAtNodeInNodalValues` in class `SoftObject`.

## Le 8 avril 2016 ##
- How is the remeshing doing? However, it does not deals with 
  possible needs for local mesh "refinements", i.e. nodes should concentrate in more complex regions while flat/simple
  regions could be described with less nodes.
