# InexactProximalOperators

Code associated to the work [`V2`](https://arxiv.org/abs/2006.06041v2), a [`new version`](https://arxiv.org/abs/2006.06041) of this work is available and the associated [`code`](https://github.com/mathbarre/InexactProximalOperators/tree/version-2) can be found in the branch "version-2".

> [1] M. Barré, A. B. Taylor and F. Bach, "Principled Analyses and Design of First-order Methods with Inexact Proximal Operators" arXiv:2006.06041[V2], 2020.

Date:    June 12, 2020

#### Authors

- [**Mathieu Barré**](https://mathbarre.github.io/)
- [**Adrien Taylor**](https://www.di.ens.fr/~ataylor/)
- [**Francis Bach**](https://www.di.ens.fr/~fbach/)

#### Content
- The [Mathematica](https://www.wolfram.com/mathematica/) notebooks in the Proofs folder allow to easily verify some statements present in the paper's proofs.
- The Matlab code in the NumericalExperiments folder allows to reproduce the numerical experiments of the paper on the CUR factorization problem and the deblurring with TV regularization.
- The Matlab code in the PEP folder allow to verify numerically the convergence rates stated in the paper using the PEP framework. It requires [YALMIP](https://yalmip.github.io/), [PESTO](https://github.com/AdrienTaylor/Performance-Estimation-Toolbox), along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).
