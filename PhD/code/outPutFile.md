% References for Output Files
% Jinming LYU
% 10 mars 2016

# `file=GeometricParameters.txt` #
The output is controlled by the method `computeGeometricParameters`, the output lists are:

1. *time* 
2. *taylorDeformation*: $D=\frac{L-B}{L+B}$, L(B) is the major(minor) axe of the ellipsoid
3. *volume*: $V=\frac{4}{3}\pi R^3=\frac{4}{3}\pi a b c$
4. *aire*: $4\pi R^2$
5. *l1*: $a$
6. *l2*: $b$
7. *l3*: $c$
8. *phi_vn*: $\theta$, inclination angle?
9. *xg*: x coordinate of mass center
10. *yg-sc.wallPosition*: height between the mass center and the wall recorded, for near-wall flow. Or *yg* for
    no wall-flow.
11. *zg*
12. *hb*: minimum height between the soft object and the wall recorded (for near-wall flow)
14. *axe*, `Veig.getColumnCount()*3` parameters  ??

# `file2=Lengths.txt` #
The ouput is controlled by the method `computeRealLengthInExtensional`, the output list are:

1. *time*
2. *L*: x坐标最大绝对值, max(abs($x_i$))
3. *l*: y坐标最大绝对值
4. *w*: z坐标最大绝对值

namely *pos_lim* (物理空间坐标)
