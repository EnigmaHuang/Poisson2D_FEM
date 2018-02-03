# Finite Element Method for Solving 2D Poisson Equation

Element type: quadrilateral, triangle

Basis (shape) function: bilinear for quadrilateral elements, linear for triangle elements

Boundary condition: Dirichlet (first-type), Robin (generalized Neumann, third-type)

Notice: the Robin BC in this program has the form: $\displaystyle \frac{\partial u}{\partial \vec{n}} + \alpha u = g, (x,y) \in \partial \Omega$ (the coefficients are normalized and merged into $\alpha$ and $g$)