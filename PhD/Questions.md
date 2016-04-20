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

## Le 9 avril 2016 ##
- `velocityComputation.java`中的`computeSingleLayer`是怎么运行的？
  因为该抽象类中只含有方法声明，但是在其他方法中却调用了该方法。

## Le 20 avril 2016 ##
- 其实对于最原始的Single-Layer积分项，对于我们当前程序所使用的积分方法，计算根本不存在奇异性（所以如何Regularized呢？）。
  这是因为对应Green函数的target point为mesh的vertex，而source point 则离散分布于单元的高斯积分点！
  二者在物理空间的距离实际上永远不可能为零。所以使用最原始的Green函数时，使用当前程序不需要对
  奇异点（单元）进行判定，然后采用如同Pozrikidis使用的Cauchy PV方法。问题是为什么使用Alexander et al. (JCP, 2014)
  的方法确能显著的提高计算精度（本质上不是减去了一个积分恒等式吗？所以改写后的与原始的应该是相等的！！！）。
- For `computeShapeFunctionAt` in `LoopElement.java`, we can find that the shape function has already take the 
  irregualr cell into consideration. Just take a good look at 

  *x0[jj]+= shape[kk]\*e.getOneRing()[kk].getNodes().get(1).getDofValues()[jj];*

  and 
  
  *coordinatesreal[n]+=shape_e[l]\*m.getElements(j).getOneRing(l).getNodes(1).getDofValues(n);*.
