# Powers-of-Forms 

This repository accompanies the paper [Unique powers-of-forms decompositions from simple Gram spectrahedra](https://arxiv.org/abs/2305.06860). It contains: 

+ a tentative implementation of Algorithm 1 therein for powers-of-forms decomposition, including subroutines to: 
	- verify whether a form is uniquely Sum-of-Squares representable.
	- compute the Sum-of-Squares support of a given SOS form.
	- decompose power sums of linear forms.
+ data from the numerical experiments of Section 4.2.  

This file describes the content of the repository. Make sure to read `installation-instructions.md`, if you want to try the code yourself. 

## Sums-of-squares representations

Running the method 
```julia
	s = calc_sos_attributes(f)
```
on a Sum of Squares polynomial `f` returns an SOSData object. These contain various informations about f, such as: 
+ `s.sosupp` is the Sum of Squares support of `f`. 
+ `s.dual_nondegenerate`: whether the supporting face of `f` in the SOS cone is exposed. 
+ `s.uniquely_sos_representable`: whether `f` is uniquely SOS representable. 

Calculations are done using semidefinite programming. Therefore, an interior point solver for SDPs is needed. If no solver is provided, the code will set the `default_solver` to `COSMO`. All experiments and the tutorial notebook explicitly set the solver to `Mosek`. Remember to change the corresponding line, if you do not have access to a `Mosek` licence, but note that free academic licences are available on the [Mosek website](https://www.mosek.com/products/academic-licenses/). The `Mosek` licence must be placed into a folder named `mosek` **in your home directory**. 

The SOSData object `s` also stores the psd matrices that were used to calculate these properties. E.g., `s.G` is a matrix in the relative interior of the Gram spectrahedron of f. `s.M_E` is the corresponding solution of the dual SDP from Appendix B. It is a relative interior point of the face $C_f$ of the dual cone, described in Appendix B. Finally, `s.G_boundary` is a matrix on the boundary of the Gram spectrahedron. It may be used to test if f is uniquely SOS representable. Indeed, the rank of `s.G_boundary` is at most the rank of `s.G`, with equality if and only if f is uniquely SOS representable. The vectors `s.λ`, `s.λ_boundary` and `s.μ` contain the eigenvalues of `s.G`, `s.G_boundary` and `s.M_E`, respectively. 

## Power sum decompositions

The function
```julia
	p = pof_decompose(f_2, f_3)
```
implements Algorithm 1 from the [paper](https://arxiv.org/abs/2305.06860), which decomposes weighted power sums 
$$
	f_d = \sum_{i=1}^m \lambda_i q_i^d, \quad d = 2,3,
$$
under certain conditions. The weights $\lambda_i$ are assumed to be positive reals. The algorithm first runs `calc_sos_attributes(f_2)` in order to obtain a basis for the Sum-of-Squares support of $f_2$. Assuming the Sum-of-Squares support has a basis of $k$-forms that are algebraically independent up to degree $3$, it then reduces to a decomposition problem as powers of **linear** forms. The function `positive_weighted_sylvester` is then called to solve for the linear forms. 
The function call `p = pof_decompose(f_2, f_3)` returns a POFDecomposition object `p`, which stores the decomposition and several intermediate quantities: 
+ `p.q` contains the addends of the decomposition.  
+ `p.λ` contains the weights of the decomposition.  

As for the intermediate quantities: 
+ `p.sosdata` is an SOSData object, storing the SOS properties of `f_2`. 
+ `p.u` is the Sum-of-Squares support of $f_2$. 
+ `M.φ_1`, `M.φ_2` and `M.φ_3` are the representing matrices of the evaluation homomorphism $\varphi\colon Y_1,\ldots, Y_N \mapsto (u_1,\ldots,u_N)$ from Section 3, on the graded components. 
+ `p.g_2` and `p.g_3` are quadratic and cubic forms, and the preimages of $f_2$ and $f_3$ under $\varphi$. 
+ `p.ℓ` stores the preimages of $q$ under $\varphi$, after they were computed using `positive_weighted_sylvester`. 



## 1-Waring decomposition

The function 
```julia
	ℓ, λ = positive_weighted_sylvester(g_2, g_3)
```
implements one version of the "simultaneous diagonalization method", see [Leurgans-Ross-Abel](https://epubs.siam.org/doi/10.1137/0614071) and [Kolda](https://www.mathsci.ai/post/jennrich/). The algorithm is described in Theorem 2.3 in the [paper](https://arxiv.org/abs/2305.06860). Assuming the quadratic form $g_2$ and the cubic form $g_3$ have a joint decomposition
$$
	g_d = \sum_{i=1}^m \lambda_i \ell_i^d, \quad d = 2,3,
$$
with **positive** weights and linearly independent linear forms $\ell_i$, then it returns this decomposition, which is the unique rank-$m$ POF decomposition of $g_2, g_3$.

Note that Theorem 2.3 allows arbitrary **nonzero** weights, whereas our implementation requires positive weights. This is an intentional choice: If all weights are positive, then the algorithm from Theorem 2.3 leads to a generalized eigenvalue problem $M_v x = \mu Mx$, where $M$ is a _positive definite_ matrix. These can be solved by the specialized Julia method `eigen(M_v, M)`, which exploits symmetry of the matrices $M_v$ and $M$. If only $g_3$ is known, then it is possible to produce an artificial canditate for $g_2$ with the function `sylvester_produce_g2`. Note that unique recovery of the weights is not possible in this case, but the linearly independent linear forms are uniquely determined up to scalar multiples.   


## Data from Numerical Experiments

The numerical experiments from Section 4.2 in the [paper](https://arxiv.org/abs/2305.06860) are protocolled in the `data` folder. Section 4.2 provides mathematical context for these experiments. Let us write $ m(n) = \lceil \frac{(n+2)(n+1)}{6} \rceil = \Theta(n^2) $, and briefly explain which data files correspond to which experiments: 
+ `data/experiment-1/explicit-family.csv` verifies unique sum-of-squares representability for an explicit family of quadratics of size $m(n)$ in $2n$ variables.
+ `data/experiment-1/explicit-family-nhood.csv` indicates that unique sum-of-squares representability also holds in a neighbourhood of the explicit parameters. 
+ `data/experiment-2/gaussian-tracefree-family.csv` indicates that sums of squares of $m(n)$ random trace-free quadratics in $2n$ variables are uniquely representable, with probability very close to $1$.
+ `data/experiment-2/gaussian-tracefree-family.csv` indicates that sums of squares of $m(n)$ random trace-free quadratics in $n$ variables are NOT uniquely representable, with probability very close to $1$. 
+ `data/experiment-3` conducts a study to estimate the probability when a sum of $m$ squares of random trace-free quadratics in $n$ variables is uniquely representable. The subfolders `m=n-1`, `m=n` etc. correspond to various fixed relations ("regimes") between $m$ and $n$. For each regime, and each $n = 2,\ldots,17$, we sample 20 batches of $m$ quadratics each and verify whether their sum of squares is uniquely representable. E.g., for the regime $m = n+1$ and $n = 12$, the results for the 20 batches is collected in `data/experiment-3/m=n+1/samples-12.csv`. The file `data/experiment-3/m=n+1/probabilities.csv` collects for each $n = 2,\ldots,17$ the probability that $m=n+1$ quadratics in $n$ variables are uniquely representable. There are some cases where $m$ exceeds the (conjectured) generic 2-Waring rank for quartics. We marked these specifically by putting 
a value of $-1$ in the probability table. 

Note that the code files in `experiments` not only protocol the results in `.csv` files, but they also save the generated SOSData objects into `.jld2` files. In principle, this provides proofs of computation, since the `.jld2` files contain, e.g., Gram matrix representations of the generated polynomials. However, since these files totaled to a size of ca. 4.6GB, we decided not to host them anywhere. Please drop me a mail, if you want to access them. 