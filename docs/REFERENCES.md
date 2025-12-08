> ⚠️ **Current Status:** This bibliography represents a **preliminary literature review**.
> While the core concepts have been identified, these papers are currently under **active study** to deepen the theoretical understanding and refine the numerical implementation.

# Bibliography & References

This project relies on the following academic papers and resources regarding **Projective Integration**, **Cardiac Electrophysiology**, and **High-Performance Computing**.

## Theoretical Foundations (Equation-Free & Projective Integration)

* **[1] Equation-free: The computer-aided analysis of complex multiscale systems**
    * *Authors:* I.G. Kevrekidis, C.W. Gear, et al.
    * *Journal:* Communications in Mathematical Sciences, Vol. 2, No. 4, pp. 715-781, 2004.
    * *Topic:* Foundational paper on the Equation-Free framework.
    * [🔗 Link to Paper](https://doi.org/10.4310/CMS.2004.v2.n4.a1)

* **[2] Equation-free: The computer-aided analysis of complex multiscale systems (Engineering View)**
    * *Authors:* I.G. Kevrekidis, C.W. Gear, G. Hummer.
    * *Journal:* AIChE Journal, Vol. 50, No. 7, pp. 1346-1355, 2004.
    * [🔗 Link to PDF](https://aiichironakano.github.io/cs653/Kevrekidis-EqFree-AIChE04.pdf)

* **[3] Projective methods for stiff differential equations: problems with gaps in their eigenvalue spectrum**
    * *Authors:* C.W. Gear, I.G. Kevrekidis.
    * *Journal:* SIAM Journal on Scientific Computing, Vol. 24, No. 4, pp. 1091-1106, 2003.
    * *Topic:* Mathematical analysis of stability for stiff systems.
    * [🔗 Link to Paper](https://doi.org/10.1137/S106482750240227X)

* **[11] Coarse-grained computations of traveling wave stability**
    * *Authors:* R. Rico-Martinez, I.G. Kevrekidis, et al.
    * *Journal:* Nonlinear Analysis: Theory, Methods & Applications, 2005.
    * *Topic:* Application of PI methods to traveling waves (relevant for cardiac pulses).
    * [🔗 Link to Paper](https://doi.org/10.1016/j.na.2004.11.025)

## Cardiac Modeling (Aliev-Panfilov & Electrophysiology)

* **[4] A simple two-variable model of cardiac excitation**
    * *Authors:* R.R. Aliev, A.V. Panfilov.
    * *Journal:* Chaos, Solitons & Fractals, Vol. 7, No. 3, pp. 293-301, 1996.
    * *Topic:* The reference paper for the model used in this project.
    * [🔗 Link to Paper](https://doi.org/10.1016/0960-0779(95)00089-5)

* **[5] Stability of periodic traveling waves in the Aliev-Panfilov reaction-diffusion system**
    * *Authors:* S. Morgan, et al.
    * *Topic:* Mathematical analysis of wave stability in the AP model.
    * [🔗 Link to Scholar](https://scholar.google.com/scholar?q=Stability+of+periodic+traveling+waves+in+the+Aliev-Panfilov)

* **[6] Mathematical Biology II: Spatial Models and Biomedical Applications**
    * *Author:* J.D. Murray.
    * *Publisher:* Springer-Verlag, 3rd Edition, 2003.
    * *Topic:* The "Bible" of mathematical biology (Chapter 1: Excitable Media).
    * [🔗 Link to Book](https://link.springer.com/book/10.1007/b98869)

* **[7] Alternans and spiral breakup in a human ventricular tissue model**
    * *Authors:* K. H. W. J. Ten Tusscher, A. V. Panfilov.
    * *Journal:* American Journal of Physiology, 2006.
    * *Topic:* Physiological context on fibrillation and spiral breakup.
    * [🔗 Link to Paper](https://doi.org/10.1152/ajpheart.00109.2006)

* **[8] Recent advances in Cardiac Modelling**
    * *Source:* National Institutes of Health (NIH), PMC11252590.
    * [🔗 Link to Article](https://pmc.ncbi.nlm.nih.gov/articles/PMC11252590/)

## Numerical Implementation & HPC

* **[9] Computing the Electrical Activity in the Heart**
    * *Authors:* J. Sundnes, G.T. Lines, X. Cai, B.F. Nielsen, K.A. Mardal, A. Tveito.
    * *Publisher:* Springer, Simula Research Laboratory, 2006.
    * *Topic:* Comprehensive guide on C++ implementation and operator splitting for cardiac PDEs.
    * [🔗 Link to PDF](https://web-backend.simula.no/sites/default/files/publications/Simula.acdc.8.pdf)

* **[10] Simulating the Aliev-Panfilov Model on GPUs**
    * *Authors:* M. Hanslien, et al.
    * *Conference:* Proceedings of the Winter Simulation Conference (WSC), 2012.
    * *Topic:* Parallel optimization strategies.
    * [🔗 Link to Paper](https://doi.org/10.1109/WSC.2012.6465046)

## Future Perspectives (AI & PINNs)

* **[12] Simulation of parametrized cardiac electrophysiology in three dimensions using physics-informed neural networks**
    * *Authors:* R.A. Gomez, J. Stöcker, B. Cansız, M. Kaliske.
    * *Source:* arXiv preprint arXiv:2506.15405, 2025.
    * *Topic:* State-of-the-art combination of Simulation and AI.
    * [🔗 Link to ArXiv](https://arxiv.org/html/2506.15405v1)
