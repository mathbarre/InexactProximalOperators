# InexactProximalOperators

Code associated to the work (coming soon), which is a new version of [`reference`](https://arxiv.org/abs/2006.06041).
The code associated to the previous version from 2020 can be found in the branch **old_version_2020**.

> [1] M. Barré, A. B. Taylor and F. Bach, "Principled Analyses and Design of First-order Methods with Inexact Proximal Operators" arXiv:2006.06041, 2020.

Date:    April 20, 2021

#### Authors

- [**Mathieu Barré**](https://mathbarre.github.io/)
- [**Adrien Taylor**](https://www.di.ens.fr/~ataylor/)
- [**Francis Bach**](https://www.di.ens.fr/~fbach/)

#### Content
- The [Mathematica](https://www.wolfram.com/mathematica/) notebooks in the Proofs folder allow to easily verify some statements present in the paper's proofs.
- The Matlab code in the PEP folder allow to verify numerically the convergence rates stated in the paper using the PEP framework. It requires [YALMIP](https://yalmip.github.io/), [PESTO](https://github.com/AdrienTaylor/Performance-Estimation-Toolbox), along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).