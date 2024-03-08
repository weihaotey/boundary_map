Numerical approximation of 2-D minimal invariant set of set-valued mapping,
representing support of stationary measure of some random dynamical system

1) "Boundary Map Invariant Manifold Method.ipynb"
unstable_mani_brute function approximates the unstable manifolds by taking a number of points on the eigenvector of a periodic point

unstable_mani function approximates the unstable manifolds of periodic points using an adaptive method which produces relatively
equally spaced sample points on the unstable manifolds. This method is based on the thesis https://spiral.imperial.ac.uk/handle/10044/1/100849

The third block defines the deterministic mappings, the corresponding boundary mapping and other functions used in the approximations.
Change the variable func_sym and func_sym_inv to specify deterministic mapping, default mapping is the classical Henon mapping.

The fourth block starts the approximation with the parameters a_sub and b_sub in the deterministic mapping f and epsilon for the magnitude of
the additive noise or radius of the ball in the set-valued mapping B_epsilon(f(x)). If 'dualoff' is 1, then the dual repeller is also
approximated; If 'compact' is 1, then only periodic points where nearby points do not diverge are computed. 'sample_grid' and 'normal_grid'
are the number of sample points to take in the original state space and the normal space respectively. 'coord_lim_grid' is the bounded region
where we search for the periodic points of the boundary mapping.

The fifth block '2d/3d better plot' plots the resulting unstable manifold approximation whether in 3-D (d2 = 0), or its projection in 2d 
(d2 =1). 

2) "Outward_normal_method_for_MIS_clean.ipynb"
Starting from a sample on the boundary of an epsilon ball of a point in the attractor, we compute the subsequent boundary
under the set-valued mapping $F(x) = \overline{B_{\varepsilon}(f(x))}$, i.e. a deterministic mapping with an additive term 
contained in the ball of size epsilon. 

By setting the formula of the original map in (func_sym) and its parameters, the program uses the boundary mapping which 
keep track of the boundary point and its outward pointing normal vector. Then, the fixed (or periodic) points on the boundary 
according to the boundary map are determined by applying the inverse boundary map. The inverse of the original function is
needed in (func_sym_inv).

The subsequent boundary is calculated by using points on the initial epsilon ball $\partial B_{\varepsilon}(x_0)$ of an initial
point $x_0$. Hence, depending on the boundary mapping, two very close point on the initial boundary might be needed to fill
the gaps on the boundary of nth iteration of the set boundary where n is large. Thus, we use mpmath to implement arbritary 
precision of float number.
