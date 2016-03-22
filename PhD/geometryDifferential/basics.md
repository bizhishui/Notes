% Basic Notes of Geometry Differential
% Jinming LYU
% 13 mars 2016

# Definitions#
- 正则性：(P.80)若曲面的参数曲线关于参数坐标的两个切向量是线性无关的，则称曲面是正则的。
- 第一基本形式：(P.96)曲面上任意点的切向量可以表示成$d\mathbf{r}(u,v)=\mathbf{r}_u(u,v)du+\mathbf{r}_v(u,v)dv$，
  其中$(du,dv)$是切向量$d\mathbf(r)(u,v)$在自然基底${\mathbf{r}_u(u,v),\mathbf{r}_v(u,v)}$下的分量。一般来说，
  ${\mathbf{r}_u(u,v),\mathbf{r}_v(u,v)}$不是单位正交基底，但是，如果知道这个基底的度量系数，则切向量与其自身的内积
  就可以表示成它分量$du,dv$的二次型。命
  \begin{equation*}
    E(u,v)=\mathbf{r}_u(u,v)\cdotp \mathbf{r}_u(u,v),F(u,v)=\mathbf{r}_u(u,v)\cdotp\mathbf{r}_v(u,v),G(u,v)=\mathbf{r}_v(u,v)\cdotp\mathbf{r}_v(u,v)
  \end{equation*}
  为基底${\mathbf{r}_u(u,v),\mathbf{r}_v(u,v)}$的度量系数，称为曲面S的第一基本量。曲面的第一基本形式定义为
  \begin{equation*}
    \begin{aligned}
    I &= d\mathbf{r}(u,v)\cdotp d\mathbf{r}(u,v) \\ 
      &= E(u,v)(du)^2+2F(u,v)dudv+G(u,v)(dv)^2 \\
      &= \begin{bmatrix}
           du & dv
         \end{bmatrix}
         \begin{bmatrix}
           E & F \\ F & G
         \end{bmatrix}
         \begin{bmatrix}
           du \\ dv
         \end{bmatrix}.
    \end{aligned}
  \end{equation*}
  第一基本形式I具有形式不变性，即与正则参数的选取无关。其几何意义是：它是切向量$d\mathbf{r}$的长度平方。
- 内积和夹角(P.99)：如在同点有另一切向量
  \begin{equation*}
     \delta\mathbf{r}(u,v)=\mathbf{r}_u(u,v)\delta u+\mathbf{r}_v(u,v)\delta v,
  \end{equation*}
  则切向量$d\mathbf{r}$和$\delta\mathbf{r}$的内积是
  \begin{equation*}
    \begin{aligned}
      d\mathbf{r}\cdotp\delta\mathbf{r} &= \frac{1}{2}((d\mathbf{r}+\delta\mathbf{r})^2-(d\mathbf{r})^2-(\delta\mathbf{r})^2) \\
                                        &= Edu\delta u + F(du\delta v + dv\delta u) +Gdv\delta v,
    \end{aligned}
  \end{equation*}
  并且
  \begin{equation*}
    \begin{aligned}
      cos\arg(d\mathbf{r},\delta\mathbf{r}) &= \frac{d\mathbf{r}\cdotp\delta\mathbf{r}}{|d\mathbf{r}||\delta\mathbf{r}|} \\
                      &= \frac{Edu\delta u + F(du\delta v + dv\delta u) +Gdv\delta v}
                      {\sqrt{E(du)^2+2Fdudv+G(dv)^2}\sqrt{E(\delta u)^2+2F\delta u\delta v+G(\delta v)^2}}
    \end{aligned}
  \end{equation*}

- 面积元素(P.101):
  \begin{equation*}
    \begin{aligned}
    |(\mathbf{r}_u\Delta u)\times (\mathbf{r}_v\Delta v)| &= |\mathbf{r}_u\times\mathbf{r}_v|\Delta u\Delta v \\
                        &= |\mathbf{r}_u||\mathbf{r}_v|sin\arg{(\mathbf{r}_u,\mathbf{r}_v)}\Delta u\Delta v \\
                        &= \sqrt{EG-F^2}\Delta u\Delta v.
    \end{aligned}
  \end{equation*}
  命 $\delta\sigma = \sqrt{EG-F^2}dudv$为曲面S的面积元素。则曲面S的面积是$A=\int\int_D\sqrt{EG-F^2}dudv$。
  它与参数选取也是无关的。

- 曲面的第二基本形式(P.139):邻近点到一点切平面的有向距离$\delta$可以用来衡量曲面在该点的弯曲程度。设$(u_0,v_0)$的邻近点
  $(u_0+\Delta u,v_0+\Delta v)$，它到$(u_0,v_0)$的切平面的距离是
  \begin{equation*}
    \begin{aligned}
      \delta (\Delta u,\Delta v) &= (\mathbf{r}(u_0+\Delta u,v_0+\Delta v))-\mathbf{r}(u_0,v_0))\cdotp\mathbf{n} \\
                                &= \frac{1}{2}(L(\Delta u)^2+2M\Delta u\Delta v+N(\Delta v)^2)+o((\Delta u)^2+(\Delta v)^2)
    \end{aligned}
  \end{equation*}
  其中
  \begin{equation*}
    \begin{aligned}
      L &= \mathbf{r}_{uu}(u_0,v_0)\cdotp\mathbf{n}(u_0,v_0) = -\mathbf{r}_u\cdotp\mathbf{n}_u, \\
      M &= \mathbf{r}_{uv}(u_0,v_0)\cdotp\mathbf{n}(u_0,v_0) = -\mathbf{r}_u\cdotp\mathbf{n}_v = -\mathbf{r}_v\cdotp\mathbf{n}_u, \\
      N &= \mathbf{r}_{vv}(u_0,v_0)\cdotp\mathbf{n}(u_0,v_0) = -\mathbf{r}_v\cdotp\mathbf{n}_v.
    \end{aligned}
  \end{equation*}
  曲面的二次微分式
  $$\text{\rom{2}}=d^2\mathbf{r}\cdotp\mathbf{n}=-d\mathbf{r}\cdotp d\mathbf{n}=L(du)^2+2Mdudv+N(dv)^2$$
  也与参数变换无关，因此将其称为曲面的第二基本形式。它是有向距离$\delta (\Delta u,\Delta v)$在
  $\sqrt{(\Delta u)^2+(\Delta v)^2}\rightarrow 0$时作为无穷小的主要部份的两倍，即
  $$\text{\rom{2}}\approx 2\delta(du,dv).$$
  
- 曲面法曲率(P.155):\bfseries{正则参数曲面在任意一个固定点，其法曲率必在两个彼此正交的切方向上分别取最大值($\kappa_1$)和最小值($\kappa_2$).}
  正则参数曲面在任意一个固定点，其法曲率取最大值和最小值的方向称为曲面在该点的主方向。
  法曲率$\kappa_n(\theta)$有如下Euler公式:
  $$\kappa_n(\theta)=\kappa_1 cos^2(\theta-\theta_0)+\kappa_2sin^2(\theta-\theta_0).$$
