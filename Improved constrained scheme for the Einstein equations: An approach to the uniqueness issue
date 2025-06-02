😃 Improved constrained scheme for the Einstein equations: An approach to the uniqueness issue

\author{
Isabel Cordero-Carrión, ${ }^{1, *}$ Pablo Cerdá-Durán, ${ }^{2, \dagger}$ Harald Dimmelmeier, ${ }^{3, \ddagger}$ José Luis Jaramillo, ${ }^{4,5, \S}$ Jérôme Novak, ${ }^{5,11}$ and Eric Gourgoulhon ${ }^{5,11}$ \ ${ }^{1}$ Departamento de Astronomía y Astrofísica, Universidad de Valencia, C/Dr. Moliner 50, E-46100 Burjassot, Valencia, Spain \ ${ }^{2}$ Max-Planck-Institut für Astrophysik, Karl-Schwarzschild-Strasse 1, D-85741 Garching, Germany \ ${ }^{3}$ Department of Physics, Aristotle University of Thessaloniki, GR-54124 Thessaloniki, Greece \ ${ }^{4}$ Instituto de Astrofísica de Andalucía, CSIC, Apartado Postal 3004, E-18080 Granada, Spain \ ${ }^{5}$ Laboratoire Univers et Théories, Observatoire de Paris, CNRS, Université Paris Diderot, 5 place Jules Janssen, F-92190 Meudon, France
}
(Received 13 September 2008; published 22 January 2009)
\begin{abstract}
Uniqueness problems in the elliptic sector of constrained formulations of Einstein equations have a dramatic effect on the physical validity of some numerical solutions, for instance, when calculating the spacetime of very compact stars or nascent black holes. The fully constrained formulation (FCF) proposed by Bonazzola, Gourgoulhon, Grandclément, and Novak is one of these formulations. It contains, as a particular case, the approximation of the conformal flatness condition (CFC) which, in the last ten years, has been used in many astrophysical applications. The elliptic part of the FCF basically shares the same differential operators as the elliptic equations in the CFC scheme. We present here a reformulation of the elliptic sector of the CFC that has the fundamental property of overcoming the local uniqueness problems. The correct behavior of our new formulation is confirmed by means of a battery of numerical simulations. Finally, we extend these ideas to the FCF, complementing the mathematical analysis carried out in previous studies.
\end{abstract}
DOI: 10.1103/PhysRevD.79.024017
PACS numbers: 04.20.Ex, 04.25.D-, 04.25.Nx, 97.60.-s
\section*{I. INTRODUCTION}
In recent years we have seen the successful application of numerical codes to accurately calculate the spacetimes of compact astrophysical objects like collapsing stellar cores, (proto)neutron stars, and black holes. Most of these codes are based on the $3+1$ formalism of general relativity (see, e.g., [1-3] for reviews). They typically fall into two classes. One approach relies on the free evolution of the $3+1$ Einstein equations, recast in order to cure longterm stability problems. Here the constraint equations are only solved initially, and closely monitored at each time step to control the accuracy of the numerical solution.
Alternatively, formulations based on a constrained evolution, where the constraints are solved in parallel with evolution equations, have proven to be successful as well. Such approaches exhibit the advantage that the solution cannot violate the constraints by definition (within the accuracy of the numerical scheme). In particular, the conformally flat approximation [4,5] (hereafter CFC) of the full Einstein equations, which constitutes a fully constrained formulation, has been shown to yield long-term stable evolutions of such astrophysical scenarios (see, e.g., [6-9]). However, apart from computational challenges,
\footnotetext{
*Isabel.Cordero@uv.es
${ }^{\dagger}$ cerda@mpa-garching.mpg.de
${ }^{\ddagger}$ harrydee@mpa-garching.mpg.de
\$jarama@iaa.es
Jerome.Novak@obspm.fr
${ }^{\text {II }}$ eric.gourgoulhon@obspm.fr
}
arising from the need to frequently solve the elliptic constraint equations, constrained formulations suffer from mathematical nonuniqueness problems when the configuration becomes too compact. In the case of the collapse of a stellar core or a (proto)neutron star to a black hole, such a situation is encountered already before the apparent horizon forms. This issue has, in the past, been prohibitive to successfully applying such formulations in numerical simulations of a wide range of astrophysical problems.
The nonuniqueness of solutions stems from the nonlinearity of the constraint equations and has been studied within the so-called extended conformal thin sandwich (XCTS) [10-12] approach to the initial data problem in general relativity. In Ref. [13] a parabolic branching was numerically found in the solutions to the XCTS equations for perturbations of Minkowski spacetime, providing the first evidence of nonuniqueness in this elliptic system. First analytical studies have been carried out in [14,15], finding support for the genericity of this nonuniqueness behavior. More specifically, the XCTS elliptic system is formed by the Einstein constraint equations in a conformal thin sandwich (CTS) decomposition [10] supplemented with an additional elliptic equation for the lapse function, which follows from the maximal slicing condition. Although no general results on the existence and uniqueness for the XCTS system are available (in contrast to the CTS case and similar elliptic systems encompassing only the constraints; see, e.g., $[1,10,11,16,17]$ ), the analysis in [14] strongly suggests the presence of a wrong sign in a certain term of the lapse equation as the culprit for the loss of
uniqueness, essentially because it spoils the application of a maximum principle to guarantee uniqueness. Moreover, in these circumstances (namely, the existence of a nontrivial kernel for the XCTS elliptic operator) it is shown in [15] that the parabolic behavior found in [13] is indeed generic.
Certain constrained evolution formalisms which incorporate elliptic gauges in their schemes contain elliptic subsystems which share essential points with the XCTS equations. Nonuniqueness in the elliptic subsystem is certainly an issue for the well-posedness of the whole elliptichyperbolic evolution system. In numerical implementations this can depend on the employed numerical scheme, in particular, on its capability to remain close to one of the solutions, at least as long as the solution stays sufficiently far from the branching point. In fact, constrained or partially constrained evolutions have shown to be robust in a variety of contexts (see, e.g., the references in [18] and Sec. 5.2.2 of [19]). However, the problems described above have also emerged, for instance, in the axisymmetric case in [20,21] (see also [22]). The analysis in [18] concludes that the reason behind the failures in these axisymmetric formulations is in fact related to the presence of wrong signs or, more precisely, to the indefinite character of certain nonlinear Helmholtz-like equations present in the scheme (see [18] for details and also for a parallel numerical discussion in terms of a class of relaxation methods for the convergence of the elliptic solvers). Regarding the full three-dimensional case, fully constrained formalisms have been presented in [23-25]. While the work in [23,24] includes an elliptic subsystem closely related to the XCTS equations and therefore suffers potentially from these nonuniqueness problems, the uniqueness properties of the scheme of [25] must yet be studied. In both cases, the full numerical performance still has to be assessed.
The goal of the present work is to discuss a scheme addressing the nonuniqueness issues of XCTS-like elliptic systems in the full three-dimensional case, with astrophysical applications as our main motivation. Having the analysis of the fully constrained formalism (hereafter FCF) of [23,24] as our ultimate aim, we focus on an approximation in the spirit of the CFC approximation by Isenberg, Wilson, and Mathews [4,26]. This methodological choice is justified since the CFC scheme already contains the relevant elliptic system of FCF, but in a setting in which potential additional problematic issues related to the FCF hyperbolic part do not mix with the specific problem we are addressing here. Therefore, we discuss in detail a modification of the CFC scheme (in the presence of matter) where maximum-principle lines of reasoning can be used to infer the uniqueness of the solutions. We investigate numerically the performance of the new CFC scheme and finally indicate the main lines for its generalization to the full Einstein FCF case.
This article is organized as follows. In Sec. II we review the FCF and CFC formalisms, and then discuss the limita-
tions found in the numerical implementations of the latter. In Sec. III we introduce the modification of the CFC scheme, with the aim of solving the uniqueness issues, and we present various numerical tests of the new scheme in Sec. IV. In Sec. V the guidelines for the generalization to the FCF case are discussed, and conclusions are drawn in Sec. VI. In the Appendix we justify a further approximation assumed in Sec. III which is consistent with the CFC setting. Throughout the paper we use the signature $(-,+,+,+)$ for the spacetime metric, and units in which $c=G=M_{\odot}=1$. Greek indices run from 0 to 3, whereas Latin ones run from 1 to 3 only.
\section*{II. THE FULLY CONSTRAINED FORMALISM AND THE CONFORMAL FLATNESS CONDITION}
\section*{A. A brief review of the fully constrained formalism}
Given an asymptotically flat spacetime ( $\mathcal{M}, g_{\mu \nu}$ ) we consider a $3+1$ splitting by spacelike hypersurfaces $\Sigma_{t}$, denoting timelike unit normals to $\Sigma_{t}$ by $n^{\mu}$. The data on each spacelike hypersurface $\Sigma_{t}$ are given by the pair $\left(\gamma_{i j}, K^{i j}\right)$, where $\gamma_{\mu \nu}=g_{\mu \nu}+n_{\mu} n_{\nu}$ is the Riemannian metric induced on $\Sigma_{t}$. We choose the convention $K_{\mu \nu}=$ $-\frac{1}{2} \mathcal{L}{n} \gamma{\mu \nu}$ for the extrinsic curvature. With the lapse function $N$ and the shift vector $\beta^{i}$, the Lorentzian metric $g_{\mu \nu}$ can be expressed in coordinates ( $x^{\mu}$ ) as
$$
\begin{equation}
g_{\mu \nu} d x^{\mu} d x^{\nu}=-N^{2} d t^{2}+\gamma_{i j}\left(d x^{i}+\beta^{i} d t\right)\left(d x^{j}+\beta^{j} d t\right) \tag{1}
\end{equation}
$$
On the other hand, we can write
$$
\begin{equation}
2 N K^{i j}=\partial_{t} \gamma^{i j}+D^{i} \beta^{j}+D^{j} \beta^{i} \tag{2}
\end{equation}
$$
where $D_{i}$ is the Levi-Civita connection associated with $\gamma_{\mu \nu}$ and $\partial_{t} \gamma^{i j}$ represents the Lie derivative with respect to the evolution vector $t^{\mu}:=\left(\partial_{t}\right)^{\mu}=N n^{\mu}+\beta^{\mu}$. As in [23] we introduce a time independent flat metric $f_{i j}$, which satisfies $\mathcal{L}{t} f{i j}=\partial_{t} f_{i j}=0$ and coincides with $\gamma_{i j}$ at spatial infinity. We define $\gamma:=\operatorname{det} \gamma_{i j}$ and $f:=\operatorname{det} f_{i j}$. This fiducial metric permits the use of tensor quantities rather than tensor densities. The next step in the formulation of [23] is the conformal decomposition of the $3+1$ fields. First, a representative $\tilde{\gamma}{i j}$ in the conformal class of $\gamma{i j}$ is chosen, so we can write
$$
\begin{equation}
\gamma_{i j}=\psi^{4} \tilde{\gamma}{i j}, \quad K{i j}=\psi^{\zeta-8} \tilde{A}^{i j}+\frac{1}{3} K \gamma^{i j} \tag{3}
\end{equation}
$$
where $K=\gamma^{i j} K_{i j}$ and $\tilde{\gamma}:=\operatorname{det} \tilde{\gamma}{i j}$, and $\zeta \in \mathbb{R}$. In Ref. [23], the choice $\zeta=4$ was adopted, leading to the following expression of $\tilde{A}^{i j}$ in terms of the lapse $N$ and shift $\beta^{i}$ :
$$
\begin{equation}
\tilde{A}^{i j}=\frac{1}{2 N}\left(\tilde{D}^{i} \beta^{j}+\tilde{D}^{j} \beta^{i}-\frac{2}{3} \tilde{D}{k} \beta^{k} \tilde{\gamma}^{i j}+\partial_{t} \tilde{\gamma}^{i j}\right) \tag{4}
\end{equation}
$$
$\tilde{D}{i}$ being the Levi-Civita connection associated with $\tilde{\gamma}{i j}$. This is in the spirit of the decomposition employed in the
(X)CTS approach to initial data. Regarding the choice of the representative of the conformal metric $\tilde{\gamma}_{i j}$, a unimodular condition $\tilde{\gamma}=f$ was adopted in [23], so that $\psi=$ $(\gamma / f)^{1 / 12}$. The deviation of the conformal metric from the flat fiducial metric is denoted by $h^{i j}$, i.e.
$$
\begin{equation}
h^{i j}:=\tilde{\gamma}^{i j}-f^{i j} . \tag{5}
\end{equation}
$$
Once the $3+1$ conformal decomposition is performed, a choice of gauge is needed in order to properly reformulate the Einstein equations as partial differential equations. The prescriptions in [23] are maximal slicing and the so-called
generalized Dirac gauge,
$$
\begin{equation}
K=0, \quad \mathcal{D}{k} \tilde{\gamma}^{k i}=0 \tag{6}
\end{equation}
$$
where $\mathcal{D}{k}$ stands for the Levi-Civita connection associated with the flat metric $f_{i j}$. The Einstein equations then become a coupled elliptic-hyperbolic system to be solved for the basic variables $h^{i j}, \psi, N$, and $\beta^{i}$ [23].
Expressing the differential operators in terms of the connection of the flat metric, the elliptic part can be written as
$$
\begin{align}
& \Delta \psi=-2 \pi E \psi^{5}-h^{k l} \mathcal{D}{k} \mathcal{D}{l} \psi+\psi \frac{\tilde{R}}{8}-\frac{\psi^{5}}{8(2 N)^{2}} \tilde{\gamma}{i k} \tilde{\gamma}{j l}\left[(L \beta)^{i j}+\frac{\partial h^{i j}}{\partial t}-\mathcal{L}{\boldsymbol{\beta}} h^{i j}-\frac{2}{3} \mathcal{D}{k} \beta^{k} h^{i j}\right] \
& \times\left[(L \beta)^{k l}+\frac{\partial h^{k l}}{\partial t}-\mathcal{L}{\boldsymbol{\beta}} h^{k l}-\frac{2}{3} \mathcal{D}{m} \beta^{m} h^{k l}\right] \tag{7}\
& \Delta(N \psi)=2 N \psi^{5} \pi(E+2 S)+N \psi \frac{\tilde{R}}{8}-h^{k l} \mathcal{D}{k} \mathcal{D}{l}(N \psi)+\frac{7}{32} \frac{\psi^{6}}{(N \psi)} \tilde{\gamma}{i k} \tilde{\gamma}{j l}\left[(L \beta)^{i j}+\frac{\partial h^{i j}}{\partial t}-\mathcal{L}{\boldsymbol{\beta}} h^{i j}-\frac{2}{3} \mathcal{D}{k} \beta^{k} h^{i j}\right] \
& \times\left[(L \beta)^{k l}+\frac{\partial h^{k l}}{\partial t}-\mathcal{L}{\boldsymbol{\beta}} h^{k l}-\frac{2}{3} \mathcal{D}{k} \beta^{k} h^{k l}\right] \tag{8}\
& \Delta \beta^{i}+\frac{1}{3} \mathcal{D}^{i} \mathcal{D}{j} \beta^{j}=16 \pi N \psi^{4} S^{i}-h^{k l} \mathcal{D}{k} \mathcal{D}{l} \beta^{i}-\frac{1}{3} h^{i k} \mathcal{D}{k} \mathcal{D}{l} \beta^{l}+\frac{\psi^{6}}{N} \mathcal{D}{j}\left(\frac{N}{\psi^{6}}\right)\left[(L \beta)^{i j}\right] \
& \quad+\frac{\psi^{6}}{N} \mathcal{D}{j}\left(\frac{N}{\psi^{6}}\right)\left[\frac{\partial h^{i j}}{\partial t}-\mathcal{L}{\boldsymbol{\beta}} h^{i j}-\frac{2}{3} \mathcal{D}{k} \beta^{k} h^{i j}\right]-2 N \Delta{k l}^{i} \tilde{A}^{k l} \tag{9}
\end{align}
$$
where $\Delta$ stands for the flat Laplacian ( $\Delta:=f^{i j} \mathcal{D}{i} \mathcal{D}{j}$ ), and $E, S^{i}$, and $S$ are, respectively, the energy density, momentum density, and trace of the stress tensor, all measured by the observer of four-velocity $n^{\mu}$ (Eulerian observer): in terms of the energy-momentum tensor $T_{\mu \nu}, E:=$ $T_{\mu \nu} n^{\mu} n^{\nu}, S^{i}:=-\gamma^{i \mu} T_{\mu \nu} n^{\nu}$, and $S:=\gamma^{i j} S_{i j}$, with $S_{i j}:=$ $T_{\mu \nu} \gamma_{i}^{\mu} \gamma_{j}^{\nu}$. Furthermore,
$$
\begin{align}
\tilde{R}= & \frac{1}{4} \tilde{\gamma}^{k l} \mathcal{D}{k} h^{m n} \mathcal{D}{l} \tilde{\gamma}{m n}-\frac{1}{2} \tilde{\gamma}^{k l} \mathcal{D}{k} h^{m n} \mathcal{D}{n} \tilde{\gamma}{m l} \tag{10}\
& (L \beta)^{i j}:=\mathcal{D}^{i} \beta^{j}+\mathcal{D}^{j} \beta^{i}-\frac{2}{3} f^{i j} \mathcal{D}{k} \beta^{k} \tag{11}\
& \Delta{i j}^{k}:=\frac{1}{2} \tilde{\gamma}^{k l}\left(\mathcal{D}{i} \tilde{\gamma}{l j}+\mathcal{D}{j} \tilde{\gamma}{i l}-\mathcal{D}{l} \tilde{\gamma}{i j}\right) \tag{12}
\end{align}
$$
Equation (7) follows from the Hamiltonian constraint, whereas Eq. (9) results from the momentum constraint together with the preservation of the Dirac gauge in time. Equation (8) corresponds to the preservation in time of the maximal slicing condition, $\partial K / \partial t=0$. Note that expression (10) for the Ricci scalar of the conformal metric does not involve any second-order derivative of the metric; this property follows from Dirac gauge [23]. The resulting elliptic subsystem coincides with the XCTS system [11], except for the field chosen to solve the maximal slicing equation: Eq. (8) above is to be solved for $N \psi$, whereas in
[11] the conformal lapse $\tilde{N}:=N \psi^{-6}$ is employed instead. This directly affects the value (and, in particular, the sign) of the power of the conformal factor in the nonlinear terms of Eqs. (7) and (8). More generally, one could define a generic rescaling of the lapse, $N=\tilde{N} \psi^{a}$, such that the choice in [11] corresponds to $a=6$, whereas the choice in Eq. (8) above corresponds to $a=-1$ (see [27] for the general equations in the vacuum case). An important remark is the absence of a choice of $a$ such that the factors multiplying $\psi$ and $\tilde{N}$ on the right-hand side of the linearized versions of Eqs. (7) and (8) both present a positive sign. In the presence of matter, terms multiplying the energy density $E$ also contribute to these sign difficulties, though in this case they can be fixed by an appropriate conformal rescaling of the energy density (see below). An additional concern in a generic evolution scenario is the sign of $\tilde{R}$, which is also relevant in the linearized equations. Implications of this issue are discussed in Sec. III.
The Einstein equations in the form of the elliptic equations (7)-(9) and the hyperbolic equation for $h^{i j}$ as given in Ref. [23] are to be solved together with the hydrodynamic equations,
$$
\begin{equation}
\nabla_{\mu}\left(\rho u^{\mu}\right)=0 \tag{13}
\end{equation}
$$
$$
\begin{equation}
\nabla_{\mu} T_{\nu}^{\mu}=0 \tag{14}
\end{equation}
$$
where $\nabla^{\mu}$ is the Levi-Civita connection associated with the metric $g_{\mu \nu}, \rho$ is the rest-mass (baryon mass) density, and $u^{\mu}$ is the four-velocity of the fluid.
\section*{B. The conformal flatness approximation}
If the hyperbolic part of the FCF system is not solved, but rather the condition $h^{i j}=0$ is imposed, the resulting three-metric $\gamma_{i j}$ is conformally flat, and the CFC approximation is recovered. Therefore, the FCF is a natural generalization of the CFC approximation. The latter has been used in many astrophysical applications, like the rotational collapse of cores of massive stars [6,28-30] or supermassive stars [8], the phase-transition-induced collapse of rotating neutron stars to hybrid quark stars [9], and equilibrium models of rotating neutron stars $[31,32]$, as well as for binary neutron star merger [7,26,33,34]. The elliptic subsystem of the FCF, Eqs. (7)-(9), reduces, in the CFC, to
$$
\begin{gather}
\Delta \psi=-2 \pi \psi^{-1}\left[E^{}+\frac{\psi^{6} K_{i j} K^{i j}}{16 \pi}\right] \tag{15}\
\Delta(N \psi)=2 \pi N \psi^{-1}\left[E^{}+2 S^{}+\frac{7 \psi^{6} K^{i j} K_{i j}}{16 \pi}\right] \tag{16}\
\Delta \beta^{i}+\frac{1}{3} \mathcal{D}^{i} \mathcal{D}{j} \beta^{j}=16 \pi N \psi^{-2}\left(S^{}\right)^{i}+2 \psi^{10} K^{i j} \mathcal{D}{j} \frac{N}{\psi^{6}} \tag{17}
\end{gather}
$$
where the following rescaled matter quantities have been introduced, following York [1]:
$$
\begin{align}
E^{} & :=\sqrt{\gamma / f} E=\psi^{6} E \tag{18}\
S^{} & :=\sqrt{\gamma / f} S=\psi^{6} S \tag{19}\
\left(S^{}\right){i} & :=\sqrt{\gamma / f} S{i}=\psi^{6} S_{i} \tag{20}
\end{align*}
$$
Equations (15) and (16) inherit the local nonuniqueness problems already present in the FCF equations. Although the sign problems specifically related to the energy density terms are solved by the conformal rescaling of the components of the energy-momentum tensor and the CFC eliminates the $\tilde{R}$ term, problems related to the $K_{i j} K^{i j}$ term remain in the scalar CFC equations. This is apparent once the extrinsic curvature is expressed in terms of the lapse and the shift.
Conformal rescaling of the hydrodynamical variables is not only relevant for local uniqueness issues. The hydrodynamic equations (13) and (14) can be formulated as a first-order hyperbolic system of conservation equations for the quantities $\left(D^{},\left(S^{}\right){i}, E^{}\right)[35,36]$, where, similarly to Eqs. (18)-(20), $D^{}:=\psi^{6} D, D:=N u^{0} \rho$ being the baryon
mass density as measured by the Eulerian observer. We can thus consider $E^{}$ and $\left(S^{}\right){i}$ as known variables in the computation of the CFC metric. Note that these quantities differ from $E$ and $S_{i}$ by a factor $\psi^{6}$, and hence it is not possible to compute the nonstarred quantities before knowing the value of $\psi$. If the energy-momentum tensor represents a fluid, then the source of Eq. (16) cannot be explicitly expressed in terms of $\left(D^{},\left(S^{}\right){i}, E^{}\right)$, the reason for that being the dependence of $S^{}$ on the pressure $P$. The pressure can only be computed in terms of the "primitive" quantities, e.g., as a function $P(\rho, \epsilon)$ of the rest-mass density and the specific internal energy $\epsilon$. The primitive quantities are, in general, recovered from ( $D, S{i}, E$ ) implicitly by means of an iteration algorithm. So far, two solutions of the problem related to the fact that $S^{*}$ directly contains $P$ have been used in numerical simulations performed with the CFC approximation.
The first approach [26] is to consider $P$, and hence also $S^{*}$, as an implicit function of $\psi$. Then Eqs. (15)-(17) can be solved as a coupled set of nonlinear equations using a fixed-point iteration algorithm. The convergence of the algorithm to the correct solution depends not only on the proximity of the initial seed metric to the solution, but also on the uniqueness of this solution. The latter point is extensively discussed in Sec. III. Furthermore, one problem of this approach is the necessity of performing the recovery of the primitive variables (which is numerically a time-consuming procedure) to compute the pressure during each fixed-point iteration. Because of the uniqueness problem, this approach can only be successfully applied in numerical simulations for, at most, moderately strong gravity (like stellar core collapse to a neutron star or the inspiral and initial merger phase of binary neutron stars), but fails for more compact configurations like the collapse of a stellar core or a neutron star to a black hole. For such scenarios with very strong gravity, one finds convergence of the metric to a physically incorrect solution of the equations or even nonconvergence of the algorithm.
A second approach to the recovery algorithm problem is the attempt to calculate $P$ independently of the CFC equations. This can be achieved by computing the conformal factor by means of the evolution equation
$$
\begin{equation}
\frac{\partial \psi^{\prime}}{\partial t}=\frac{\psi^{\prime}}{6} \mathcal{D}_{k} \beta^{k} \tag{21}
\end{equation}
$$
The conformal factor $\psi^{\prime}$ obtained in this way is analytically identical to the $\psi$ from Eqs. (15)-(17), but here we use a different notation to keep track of the way it is computed. The value of $\psi^{\prime}$ is solely used to evaluate $P$, and the coupled system of Eqs. (15)-(17) is solved to determine $\psi, N$, and $\beta^{i}$. Although this approach allows one to avoid the problem of recovering the primitive variables at each iteration, it also suffers from the convergence problem, and the simulation of configurations with very strong gravity is still not feasible. Furthermore, new com-
plications are introduced by using two differently computed values, $\psi$ and $\psi^{\prime}$, of the same quantity. For some scenarios, like the formation of a black hole from stellar collapse, the numerical values of these two quantities during the evolution of the system start to diverge significantly at some point. We find that this inconsistency cannot be avoided, since any attempt to artificially synchronize both values leads to numerical instabilities.
\section*{III. THE NEW SCHEME IN THE CONFORMALLY FLAT CASE}
\section*{A. Uniqueness of the elliptic equations and convergence of elliptic solvers}
Well-posed, elliptic partial differential systems admit nonunique solutions whenever the associated differential operator has a nontrivial kernel. When discussing sufficient conditions guaranteeing uniqueness, it is illustrative to first consider the case of a scalar elliptic equation. In particular, for the class of scalar elliptic equations for the function $u$ of the form
$$
\begin{equation}
\Delta u+h u^{p}=g, \tag{22}
\end{equation}
$$
where $h$ and $g$ are known functions independent of $u$, a maximum principle can be used to prove local uniqueness of the solutions as long as the sign of the exponent $p$ is different from the sign of the proper function $h$ [1,37-39].
In the CFC case, we are not dealing with a single scalar elliptic equation, but rather with the coupled nonlinear elliptic system (15)-(17). Therefore, assessing whether or not the scalar equations (15) and (16) present good signs for the application of a maximum principle is an important step for understanding the uniqueness properties of the whole system. However, as pointed out in the previous section, the CFC equations for the conformal factor and the lapse possess the wrong signs in the quadratic extrinsic curvature terms (once everything is expressed in terms of the lapse and the shift). This problem can be fixed in Eq. (15) by an appropriate rescaling of the lapse, $N=$ $\tilde{N} \psi^{6}$, but this strategy does not solve the problem for the lapse equation (cf. the discussion on the conformal lapse $\tilde{N}$ in Sec. II A). Therefore, we cannot use the maximum principle to infer local uniqueness of the solutions to the CFC equations. In these conditions of potential nonunique solutions, convergence to an undesirable solution may happen. As mentioned in the Introduction, this pathology has been illustrated using simple analytical examples of scalar equations of the type (22) in [14], as well as in numerical implementations of the vacuum Einstein constraints in the XCTS approach [13] and certain constrained evolution formalisms (see, e.g., [18]).
In the context of the CFC approximation this sign issue has also appeared, in particular, associated with the "recovery algorithm" problem discussed in Sec. II B since it involves the evaluation of the conformal factor. Nonunique solutions of $\psi$, either due to the use of the nonconformally
rescaled $E$ or the quadratic extrinsic curvature term, spoil the convergence of the algorithm when density, and thus compactness, increases. We again emphasize that a possible synchronization of $\psi$ and $\psi^{\prime}$ does not solve the problem in general, since numerical instabilities eventually arise at sufficiently high compactness.
\section*{B. Numerical examples}
The nonuniqueness of solutions has also been observed in the FCF, as described in the following example. Let us consider a vacuum spacetime, with initial data formed by a Gaussian wave packet, as in [23], but with a much higher amplitude, $\chi_{0}=0.9$ instead of $\chi_{0}=10^{-3}$ in [23] (see the latter reference for notations). The integration technique and numerical settings are the same as in [23], but contrary to the results for small amplitudes obtained in that reference, the wave packet does not disperse to infinity and instead starts to collapse. Figure 1 displays the time evolution of the central lapse $N_{\mathrm{c}}$ at $r=0$ and of the system's Arnowitt-Deser-Misner (ADM) mass $M_{\text {ADM }}$, which in the present conformal decomposition can be expressed as

![](https://cdn.mathpix.com/cropped/2025_06_02_4610231bf4be6dd7bda1g-05.jpg?height=1162&width=863&top_left_y=1113&top_left_x=1092)


FIG. 1. Time evolution of the central lapse $N_{\mathrm{c}}$ (top panel) and the ADM mass $M_{\text {ADM }}$ (bottom panel) for a collapsing packet of gravitational waves, using the integration scheme proposed in [23]. The unit of $t$ is given by the initial width of the wave packet.
$$
\begin{align}
M_{\mathrm{ADM}} & =-\frac{1}{2 \pi} \oint_{\infty}\left(\mathcal{D}{i} \psi-\frac{1}{8} \mathcal{D}^{j} \tilde{\gamma}{i j}\right) d \mathcal{A}^{i} \
& =-\frac{1}{2 \pi} \oint_{\infty} \mathcal{D}_{i} \psi d \mathcal{A}^{i} \tag{23}
\end{align}
$$
where the integral is taken over a sphere of radius $r=\infty$ and the second equality follows from the use of Dirac gauge [Eq. (6)].
The very sudden change at $t \simeq 0.4$ in both the central lapse and the ADM mass, which is also present in, e.g., the central conformal factor $\psi_{c}$, originates from the convergence of the elliptic system (7)-(9) to another solution with a different (unphysical) value of the ADM mass. The good conservation of $M_{\text {ADM }}$ and the smooth evolution of $N_{\mathrm{c}}$ for $t \gtrsim 0.4$ indicate that this other solution remains stable until $t \simeq 2$, when high-frequency oscillations appear. These oscillations may be due to the overall inconsistency of the system, destabilizing the whole scheme. On the other hand, the time evolution of $h^{i j}$ does not show any such type of behavior, and $h^{i j}$ exhibits a continuous radial profile at all times. This is numerical evidence that, also for the full Einstein case (i.e. without the approximation), the generalized elliptic equations suffer from a convergence problem similar to the CFC case.
The same subject is also exemplified when one tries to calculate the spacetime metric for an equilibrium neutron star model from the unstable branch using either Eqs. (7)(9) in the FCF case or Eqs. (15)-(17) in the CFC approximation. Even for the simple setup of a polytrope with an adiabatic index $\Gamma=2$ in spherical symmetry, those metric equations yield-when converging at all-a grossly incorrect solution if the matter quantities ( $D, S_{i}, E$ ) in the source terms are held fixed. Both the metric components and the ADM mass can deviate from the physical solution by a few tens of percents, even though that incorrect metric satisfies the asymptotic flatness condition. The reason why programs for constructing rotating relativistic neutron star models, like the KEH code [40], the RNS code [41], or the BGSM code [42], are not obstructed by this nonuniqueness problem is apparently that they all utilize an iteration over both the metric and the hydrodynamic equations simultaneously, thereby allowing the matter quantities to change during the calculation of the metric.
We want to stress here that these nonconvergence issues in the CFC case are not related to the approximation that is made. If one considers this system in the spherical (onedimensional) case, the CFC is no longer an approximation, but is the choice of the so-called isotropic gauge. Even then, the elliptic system (15)-(17) no longer converges to the proper (physical) solution.
\section*{C. The new scheme and its theoretical properties}
Despite the above-mentioned convergence problems, numerically simulating the physical problem of spherical collapse to a black hole in isotropic coordinates has been
successfully studied by Shapiro and Teukolsky in [43]. Because of the spherical symmetry, there exists only one independent component of the extrinsic curvature. It is then possible to compute directly a conformal extrinsic curvature, $\psi^{6} K_{r}^{r}$, from the conserved hydrodynamical variables. The elliptic equation for $\psi$ then decouples from the other elliptic equations by introducing this conformal extrinsic curvature and using the conserved hydrodynamical variables in the source. This source term presents no problem for proving local uniqueness, and the equation for $\psi$ always converges to the physically correct solution. Once the conformal factor, the extrinsic curvature (from the conformal factor and the conformal extrinsic curvature), and the conserved hydrodynamical variables are known, the elliptic equation for $N \psi$ can be solved and, again, the source exhibits no local uniqueness problem. This follows from the fact that the extrinsic curvature is not expressed in terms of the lapse and the shift. This contrasts with the CFC equation (16) where a division by $N^{2}$ occurs in the last term when the extrinsic curvature is expressed in terms of its constituents $N, \psi$, and $\beta^{i}$. In addition, there is no need to use $\psi^{\prime}$. Finally, the elliptic equation for the shift vector can be solved. In summary, no problems of instability or divergence are encountered.
We now generalize this scheme to the CFC case in three dimensions. This involves the use of two different conformal decompositions of the extrinsic curvature: first, two different conformal rescaling and, second, two different decompositions of the traceless part into longitudinal and transverse parts. Adopting maximal slicing, $K=0$, a generic conformal decomposition can be written as
$$
\begin{equation}
K^{i j}=\psi^{\zeta-8}\left(A^{(\zeta)}\right)^{i j}:=\psi^{\zeta-8}\left(\frac{1}{\sigma}(L X)^{i j}+A_{\mathrm{TT}}^{i j}\right) \tag{24}
\end{equation}
$$
where $\zeta$ is a free parameter, $\sigma$ is a free function, $A_{\text {TT }}^{i j}$ is transverse traceless, and $L$ is the conformal Killing operator defined by Eq. (11). We implicitly make use of a flat conformal metric, with respect to which $A_{\text {TT }}^{i j}$ is transverse, although, in principle, it would be more general to use the metric $\tilde{\gamma}^{i j}$ and the conformal Killing operator associated with it, $\tilde{L}$. But such a decomposition would introduce many technical difficulties in our treatment. In particular, it is numerically easier to handle tensors which are divergencefree with respect to the flat metric in the generalization to the FCF. The vector $X^{i}$, on which $L$ is acting, is therefore called the longitudinal part of $\left(A^{(\zeta)}\right)^{i j}$. The first decomposition we use is the one introduced in Eqs. (3) and (4) with the choice $\zeta=4$ and $\sigma=2 N$. This corresponds to a CTS-like decomposition of the traceless part, so that $X^{i}$ is given by the shift vector $\beta^{i}$ and $A_{\mathrm{TT}}^{i j}$ can be expressed in terms of the time derivative of the conformal metric. We denote this traceless part as $\tilde{A}^{i j}:=\left(A^{(4)}\right)^{i j}$. In the CFC approximation this becomes
$$
\begin{equation}
K^{i j}=\psi^{-4} \tilde{A}^{i j}, \quad \tilde{A}^{i j}=\frac{1}{2 N}(L \beta)^{i j} \tag{25}
\end{equation}
$$
The second conformal decomposition,
$$
\begin{equation}
K^{i j}=\psi^{-10} \hat{A}^{i j}, \quad \hat{A}^{i j}=(L X)^{i j}+\hat{A}_{\mathrm{TT}}^{i j}, \tag{26}
\end{equation}
$$
refers to $\zeta=-2$ and $\sigma=1$. It instead corresponds to a conformal transverse traceless (CTT) decomposition of the traceless part of the extrinsic curvature introduced by Lichnerowicz [44]. Notice that we have defined $\hat{A}^{i j}:=$ $\left(A^{(-2)}\right)^{i j}$, not to be confused with $\tilde{A}^{i j}:=\left(A^{(4)}\right)^{i j}$. The relation between $\hat{A}^{i j}$ and $\tilde{A}^{i j}$ is given by
$$
\begin{equation}
\hat{A}^{i j}=\psi^{10} K^{i j}=\psi^{6} \tilde{A}^{i j} \tag{27}
\end{equation}
$$
In terms of $\hat{A}^{i j}$, the CFC momentum constraint can be written as
$$
\begin{equation}
\mathcal{D}{j} \hat{A}^{i j}=8 \pi \psi^{10} S^{i}=8 \pi \psi^{6} f^{i j} S{j}=8 \pi f^{i j} S_{j}^{} \tag{28}
\end{equation*}
$$
Consistency between the CTT-like decomposition (26) and the CTS-like one (25) generically requires a nonvanishing transverse part, $\hat{A}{\text {TT }}^{i j}$ in Eq. (26). However, as it is shown in the Appendix, this $\hat{A}{\text {TT }}^{i j}$ is smaller in amplitude than the nonconformal part $h^{i j}$ of the spatial metric, and $\hat{A}^{i j}$ can be approximated on the CFC approximation level as
$$
\begin{equation}
\hat{A}^{i j} \approx(L X)^{i j}=\mathcal{D}^{i} X^{j}+\mathcal{D}^{j} X^{i}-\frac{2}{3} \mathcal{D}_{k} X^{k} f^{i j} \tag{29}
\end{equation}
$$
From Eqs. (26) and (28), an elliptic equation for the vector $X^{i}$ can be derived,
$$
\begin{equation}
\Delta X^{i}+\frac{1}{3} \mathcal{D}^{i} \mathcal{D}{j} X^{j}=8 \pi f^{i j} S{j}^{} \tag{30}
\end{equation*}
$$
from which $X^{i}$ can be obtained. With this vector field, one can calculate the tensor $\hat{A}^{i j}$ via (29). Notice that in the case of spherical symmetry, $\hat{A}^{r r}=\psi^{10} K^{r r}=\psi^{6} K^{r}{ }_{r}$ is the quantity used by Shapiro and Teukolsky [43].
The elliptic equation for the conformal factor can be rewritten in terms of the conserved hydrodynamical variables and $\hat{A}^{i j}$ :
$$
\begin{equation}
\Delta \psi=-2 \pi \psi^{-1} E^{}-\psi^{-7} \frac{f_{i l} f_{j m} \hat{A}^{l m} \hat{A}^{i j}}{8} \tag{31}
\end{equation*}
$$
This equation can be solved in order to obtain the conformal factor. Once the conformal factor is known, the procedure to implicitly recover the primitive variables from the conserved ones is possible, the pressure $P$ can be computed using the equation of state, and therefore $S^{*}$ is at hand. The elliptic equation for $N \psi$ can be reformulated by means of the conserved hydrodynamical variables, $\hat{A}^{i j}$, and the conformal factor:
From this equation $N \psi$ can then be obtained and, consequently, so can the lapse function $N$. Note that, since $\hat{A}^{i j}$ is already known at this step, no division by $N^{2}$ spoils the good sign for the maximum principle.
Using the relation between the two conformal decompositions of the extrinsic curvature, $\hat{A}^{i j}=\psi^{6} \tilde{A}^{i j}$, Eq. (25) can be expressed as $(L \beta)^{i j}=2 N \psi^{-6} \hat{A}^{i j}$. Taking the divergence, we arrive at an elliptic equation for the shift vector,
$$
\begin{equation}
\Delta \beta^{i}+\frac{1}{3} \mathcal{D}^{i}\left(\mathcal{D}{j} \beta^{j}\right)=\mathcal{D}{j}\left(2 N \psi^{-6} \hat{A}^{i j}\right) \tag{33}
\end{equation}
$$
where the source is completely known. This elliptic equation can be solved in order to obtain the shift vector $\beta^{i}$ consistent with $\partial_{t} \tilde{\gamma}_{i j}=0$, as required by the CFC approximation.
In this recast form of the CFC equations, an extra elliptic vectorial equation for the vector field $X^{i}$ is introduced. However, now the signs of the exponents of $\psi$ and $N$ are compatible with the maximum principle for scalar elliptic equations, and the problem is linearization stable. While this does not guarantee global uniqueness of the solutions, it provides a sufficient result for local uniqueness. This strongly relies on the fact that the system decouples in a hierarchical way, which we summarize here once more:
(1) With the hydrodynamical conserved quantities at hand, solve Eq. (30) for $X^{i}$, and thus for $\hat{A}^{i j}$.
(2) Solve Eq. (31) for $\psi$, where local uniqueness is now guaranteed. Then $S^{*}$ can be calculated consistently.
(3) Solve Eq. (32) for $N \psi$, a linear equation where the maximum principle can be applied and uniqueness and existence follow with appropriate boundary conditions.
(4) As the source of Eq. (33) is then fully known, solve it for $\beta^{i}$.
Note that this scheme is similar to that used by Shibata and Uryū [45] to compute initial data for black hole-neutron star binaries. We will discuss this point further in Sec. VIB.
The new CFC metric equations presented here not only allow us to evolve the hydrodynamical equations and recover the metric variables from the elliptic equations in a consistent way (no auxiliary quantity $\psi^{\prime}$ is needed), but they also permit us to introduce initial perturbations in the hydrodynamical variables (strictly speaking, in the conserved quantities) in a set of previously calculated initial data and directly deliver the correct values for the metric. It is even possible to perturb only the primitive quantities, and consistently resolve for the metric by iterating until the conformal factor $\psi$, which links the primitive to the conserved quantities, converges. We find that such an iteration
method fails for sufficiently strong gravity if the original CFC formulation is used.
\section*{IV. NUMERICAL RESULTS}
We recapitulate that the original CFC formulation exhibits serious convergence problems when dealing with highly compact configurations such as nascent black holes. This weakness of the original formalism is noticeable in the fact that no simulations of rotational collapse to a black hole substantially beyond the formation of the apparent horizon have been performed so far in the CFC. Furthermore, some scenarios which do not involve the formation of a black hole are already feasible with the old formulation only if procedures like using Eq. (21), with all associated problems and inconsistencies, are employed. An example is the migration of a neutron star model from the stable to the unstable branch, which is a standard test for relativistic hydrodynamics codes. In contrast, the new CFC scheme presented in this work solves all the problems that prevented performing such simulations in the past. In order to show the suitability of the new scheme, we present the results of numerical simulations of the migration test and of the rotational collapse to a black hole.
\section*{A. Model setup}
The numerical simulations presented here are performed using the numerical code COCONUT [28,46]. This code solves the evolution of the hydrodynamics equations coupled to the elliptic equations for the spacetime metric in the CFC approximation. Standard high-resolution shock-capturing schemes are used in the hydrodynamic evolution, while spectral methods are employed to solve the metric equations. The code is based on spherical polar coordinates, and for the tests presented here we assume axisymmetry and symmetry with respect to the equatorial plane. Note that the metric equations presented in this paper are covariant. Thus the formalism can be used for any coordinate basis as well as without any symmetry conditions.
The initial models are general relativistic $\Gamma=2$ polytropes in equilibrium with a polytropic constant $K=100$. The models are chosen to be situated on the unstable branch, i.e. $\partial M_{\mathrm{ADM}} / \partial \rho_{\mathrm{c}}<0$, where $\rho_{\mathrm{c}}$ is the central rest-mass density. Therefore, any perturbation of the star induces either a collapse to a black hole or migration to a configuration of the same baryon mass on the stable branch. Table I shows the main features of these initial models. Models D1 to D4 are uniformly rotating models which are identical to those presented in [47]. The model labeled SU is a spherical model, while the model labeled SS is the counterpart model with the same baryon mass but it is located on the stable branch. The equilibrium rotating star models in Dirac gauge (the axisymmetric and stationary limit of the FCF) used here are described in [48],
TABLE I. Initial models used in the migration test and the rotational collapse to a black hole. $\rho_{\mathrm{c}, \mathrm{i}}$ is the initial central restmass density, $\Omega_{\mathrm{i}}$ is the initial angular velocity, $r_{\mathrm{p}, \mathrm{i}} / r_{\mathrm{e}, \mathrm{i}}$ is the initial ratio of the polar to the equatorial coordinate radius, $M_{\text {ADM }}$ is the gravitational ADM mass, and $J$ is the total angular momentum (which is conserved in the CFC during the evolution in the axisymmetric case). Units in which $G=c=M_{\odot}=1$ are used.
\begin{tabular}{lcccccc}
\hline \hline Model & $\rho_{\mathrm{c}, \mathrm{i}}\left[10^{-3}\right]$ & $\Omega_{\mathrm{i}}\left[10^{-2}\right]$ & $r_{\mathrm{p}, \mathrm{i}} / r_{\mathrm{e}, \mathrm{i}}$ & $r_{\mathrm{e}, \mathrm{i}}$ & $M_{\mathrm{ADM}}$ & $J / M_{\mathrm{ADM}}^{2}$ \
\hline SU & 8.000 & 0 & 1.00 & 4.267 & 1.447 & 0 \
SS & 1.346 & 0 & 1.00 & 7.999 & 1.424 & 0 \
D1 & 3.280 & 1.73 & 0.95 & 5.947 & 1.665 & 0.207 \
D2 & 3.189 & 2.88 & 0.85 & 6.336 & 1.727 & 0.362 \
D3 & 3.134 & 3.55 & 0.75 & 6.839 & 1.796 & 0.468 \
D4 & 3.116 & 3.95 & 0.65 & 7.611 & 1.859 & 0.542 \
\hline \hline
\end{tabular}
and are computed using the lorene [49] library. We map the hydrodynamic and metric quantities to the CFC code neglecting the $h^{i j} \sim 10^{-3}$ terms, which are negligible due to their smallness. Alternatively, we compute CFC equilibrium initial models. In this case we find that the differences with respect to the FCF models are small ( $\sim 0.1 \%$ ) for representative metric and hydrodynamic quantities, initially and during the evolution, and therefore we discuss only the FCF initial models here.
The hydrodynamic equations are discretized on the finite difference grid with $n_{r} \times n_{\theta}$ grid points. The radial grid size is $\Delta r_{0}$ for the innermost cell and increases geometrically outwards, while the angular grid is equidistantly spaced. The metric equations are solved on a spectral grid consisting of $n_{\mathrm{d}}-1$ radial domains distributed such as to homogeneously cover the finite difference grid and a compactified exterior domain extending to radial infinity. On the spectral grid we resolve each radial domain with 33 collocation points. The spherical model needs only one angular collocation point, while we use 17 angular points for the rotating models.
We track the location of the apparent horizon by means of a three-dimensional spectral apparent horizon finder, described in detail and tested in [50]. The apparent horizon location is given by a function $\mathcal{H}(r, \theta)$, which is decomposed into a set of spherical harmonics. The coefficients of $\mathcal{H}$ in this basis are computed iteratively, in order to satisfy the condition that the expansion in the outgoing null direction vanishes at the apparent horizon location.
\section*{B. Migration of unstable neutron stars to the stable branch}
The first test we consider is the migration of a neutron star model in equilibrium from the unstable branch to the stable branch, which is a standard but still demanding test for general relativistic hydrodynamic codes, as it involves the dynamic transition between two very compact equilibrium states. This test has been performed in the past in full general relativistic simulations [51]. We start the evolution
with the nonrotating equilibrium model labeled SU. Since it belongs to the unstable branch, any perturbation from exact equilibrium (which can, for instance, be caused by discretization errors) leads either to a collapse or to an expansion to a new equilibrium configuration of the same baryon mass on the stable branch. The corresponding equilibrium configuration with the same baryon mass, model SS, has a smaller ADM mass than the initial system (see Table I). Therefore, to preserve the ADM mass, the final configuration cannot be exactly the equilibrium model given by SS. The energy difference between the SU and SS models should be transformed into kinetic energy, remaining in the final object in the form of pulsations.
In our case the numerical truncation error is sufficient to trigger the migration. Since the final neutron star on the stable branch is larger than the initial model (see Table I), the outer boundary of the finite difference grid is chosen to be 4.5 times the radius of the SS model. We perform two simulations on a finite difference grid with 150 or 300 radial cells and $\Delta r_{0}=0.022$ or 0.012, respectively. We use $n_{\mathrm{d}}=6$ radial domains for the spectral grid. We evolve the system with either a polytropic or an ideal gas equation of state.
Figure 2 shows the time evolution of the central values of the rest-mass density and the lapse. As the star expands, $\rho_{\text {c }}$ decreases while $N_{\text {c }}$ grows until the new stable equilibrium configuration is reached. In the polytropic case, there are no physical mechanisms to damp the strong pulsations, and the final state resembles a star oscillating around the equilibrium configuration until numerical dissipation finally damps the oscillations. This can be seen in the pulsating values of rest-mass density and lapse around the value corresponding to the equilibrium model on the stable branch (solid horizontal line in Fig. 2).
In the ideal gas case, shock waves are formed at every pulsation, and they dissipate kinetic energy into thermal energy, thereby damping the oscillations. As these shock waves reach the surface of the star, a small amount of mass is expelled from the star and matter is ejected outwards into the surrounding artificial low-density atmosphere until it leaves the grid across the outer numerical boundary. We approximately compute the escape velocity as $v_{\mathrm{e}}=$ $\sqrt{2 U} \approx \sqrt{\psi^{2}-1}$, where $U$ is the Newtonian potential. This formula is not exact in general relativity, but it should by sufficiently accurate near the outer numerical boundary where gravity is weaker. We find that the shock waves leaving the computational domain exceed the escape velocity, and therefore the lost mass is gravitationally unbounded. We also check that these results are not affected by changing the resolution or setting the outer boundary twice as far away. As the oscillations are damped, the shock waves become weaker and the mass expelled at each oscillation is smaller. At the end of the simulation the star has lost about $10 \%$ of its initial baryon mass, approaching a state of constant baryon mass. As a conse-

![](https://cdn.mathpix.com/cropped/2025_06_02_4610231bf4be6dd7bda1g-09.jpg?height=1164&width=866&top_left_y=155&top_left_x=1090)


FIG. 2. Time evolution of the central rest-mass density $\rho_{\mathrm{c}}$ (top panel) and the central lapse $N_{\mathrm{c}}$ (bottom panel) for the migration of the unstable neutron star model, SU, to the stable branch, with either a polytropic (solid lines) or an ideal gas (dashed lines) equation of state. The dotted horizontal lines mark the values of $\rho_{\mathrm{c}}$ and $N_{\mathrm{c}}$ for the SS equilibrium configuration from the stable branch with the same baryon mass $M_{\mathrm{b}}$ as the SU model, while the dash-dotted lines are obtained from a series of equilibrium models where mass shedding, like in the migration model with an ideal gas equation of state, is taken into account. In the inset the baryon mass $M_{\mathrm{b}}$ versus the $\rho_{\mathrm{c}}$ relation for this model setup is displayed. The SU (the initial model) and SS (the final state for a polytropic equation of state) models as well as the final state for an ideal gas equation of state are marked. The arrows symbolize the respective migration paths.
quence, the final equilibrium configuration on the stable branch is not the SS model anymore but, rather, the corresponding model from the stable branch with lower baryon mass and central density. In Fig. 2 we plot the central restmass density and lapse of a series of equilibrium models on the stable branch corresponding to the baryon mass remaining in the computational domain at each time. It can be seen that these values deviate with time from the SS model and fit the final state in the hydrodynamical evolution of the star.
As a by-product of this study we draw the reader's attention to the consistency (as it should be) between the
amplitude and the frequency of the oscillations. The period of these oscillations is approximately of the order of the hydrodynamical characteristic time $\tau_{\rho}$, which decreases with density like $\tau_{\rho} \approx \rho^{-1 / 2}$. In the polytropic case, the maxima of the oscillations in $\rho_{\mathrm{c}}$ are systematically higher than in the ideal gas case. Consequently, the characteristic time is shorter than in the ideal gas case, as Fig. 2 shows. A second property worth pointing out is that the low numerical viscosity of our code is responsible for maintaining a nearly constant amplitude of the oscillations (in the polytropic case) during many characteristic times.
Our simulations are consistent with the results from the fully relativistic three-dimensional code in [51]. Similar simulations of this test, with the original, unmodified CFC scheme, lead to a completely incorrect solution with a grossly incorrect ADM mass. When running with the new improved CFC scheme, we obtain $M_{\text {ADM }}=$ $1.451 M_{\odot}$ and initial values for the conformal factor and lapse of $\psi_{\mathrm{c}}=1.561$ and $\alpha_{\mathrm{c}}=0.273$, respectively. On the other hand, with the unmodified conventional CFC scheme, the metric solver already initially converges to a solution with $M_{\mathrm{ADM}}=0.647 M_{\odot} \quad(55 \%), \quad \psi_{\mathrm{c}}=1.221$ (61\%), and $\alpha_{\mathrm{c}}=0.532$ (63\%), where the relative differences from the physically correct solution are given in parentheses.
As presented in [52] the migration test can be successfully simulated using the old CFC scheme, if one resorts to additionally solving the evolution equation (21) for the conformal factor (which would lead to large inconsistencies in scenarios with higher compactness but still yield acceptable results for the standard migration case). Here the superiority of the new, fully consistent CFC scheme, which does not depend on such scenario-dependent amendments, already becomes apparent.
\section*{C. Collapse of unstable neutron stars to a black hole}
As the second test, we present the collapse of a (spherical or rotating) neutron star model to a black hole. Following [47] we trigger the collapse to a black hole by reducing the polytropic constant $K$ by $2 \%$ in the initial D1 to D4 models. Alternatively, in the spherical SU model we increase the rest-mass density by $0.1 \%$, which yields a similar dynamic evolution. However, since the models are initially in equilibrium, the total collapse time depends strongly on the perturbation applied. In these cases, the outer boundary of the finite difference grid is $20 \%$ larger than the star radius. For the spherical SU model, we perform two simulations using 150 or 300 radial cells and $\Delta r_{0} \sim 10^{-3}$ or $10^{-4}$, respectively, to assess the resolution dependence of our simulations. For the rotating D1 to D4 models the grid is made up of $150 \times 20$ and $150 \times$ 40 cells, with the same radial grid spacing as in the spherical model. We choose $n_{\mathrm{d}}=8$ radial domains for the spectral grid. As in [47] we use a polytropic equation of state in the evolution.
The top panel of Fig. 3 shows the evolution of the restmass density and lapse at the center. Since for the maximal slicing condition the singularity cannot be reached in a finite time, $N_{\mathrm{c}}$ rapidly approaches zero once the apparent horizon has formed. In parallel, $\rho_{\mathrm{c}}$ grows, which results in a decrease in the numerical time step due to the Courant condition applied to the innermost grid cell. We terminate the evolution as the central regions of the collapsing star inside the apparent horizon become increasingly badly resolved on the regular grid, and thus numerical errors grow. We check in the SU model that by refining the radial resolution we are able to follow the collapse to even higher densities. Therefore, the only limitation to perform a stable evolution after the apparent horizon formation is the nu-

![](https://cdn.mathpix.com/cropped/2025_06_02_4610231bf4be6dd7bda1g-10.jpg?height=1242&width=871&top_left_y=805&top_left_x=1088)


FIG. 3. Collapse to a black hole for the spherical model SU, and the rotating models D1 and D4. The top panel shows the time evolution of the central lapse $N_{\mathrm{c}}$ (thin lines) and the central restmass density $\rho_{\mathrm{c}}$ relative to the initial value $\rho_{\mathrm{c}, 0}$ (thick lines). The bottom panel shows the time evolution of the apparent horizon radius $r_{\mathrm{AH}, \mathrm{e}}$ in the equatorial plane (thin lines) and the rest mass $M_{\text {outside AH }}$ remaining outside the apparent horizon relative to the total rest mass $M$ (thick lines). The dashed vertical lines mark the time when the apparent horizon first appears. If the axes of the lower panel were exchanged, the resulting plot would resemble the typical spacetime diagram of a star collapsing to a black hole.
merical resolution used. Note, however, that the spatial gauge condition is fixed in the CFC, and thus we are not able to utilize the common method of exploiting the gauge freedom for the radial component of the shift vector in order to effectively increase the central resolution.
In the bottom panel of Fig. 3 we display the time evolution of the apparent horizon radius. As expected, the apparent horizon appears at a finite radius and already encompasses a significant fraction of the total mass of the star ( $\sim 70 \%-80 \%$ ) at that time. Afterwards, its radius grows as the surrounding matter falls inside beyond the horizon. The fraction of the rest mass remaining outside the horizon is also plotted in the figure. In the rotating case the apparent horizon is slightly nonspherical. The ratio of the polar to the equatorial proper circumferential radius of the apparent horizon at the end of the simulation is $R_{\mathrm{p}} / R_{\mathrm{e}}=0.998-0.978$ for models D 1 to D 4 , where $R_{\mathrm{e}}:=$ $\int_{0}^{2 \pi} \sqrt{g_{\varphi \varphi}} d \varphi /(2 \pi)$ and $R_{\mathrm{p}}:=\int_{0}^{\pi} \sqrt{g_{\theta \theta}} d \theta / \pi$.
Since we cannot reasonably determine the location of the event horizon, as this would require the evolution of spacetime until the black hole has become practically stationary, we utilize the apparent horizon radius to estimate the mass of the newly formed black hole. Following the prescription in [47] we use the expression $M_{\mathrm{BH}}=R_{\mathrm{e}} / 2$. Note that this formula is only strictly valid for a stationary Kerr black hole. In our case, however, first of all, some (albeit small) amount of matter is still outside the horizon and the black hole is still dynamically evolving, and second, the metric of a Kerr black hole is not conformally flat [53]. Still, according to [47] this approximation (excluding the effects of the CFC) introduces an error in the mass estimate of only $\sim 2 \%$. For the spherical model the estimated value for $M_{\mathrm{BH}}$ at the end of the simulation agrees within $0.5 \%$ with the ADM mass $M_{\text {ADM }}$ of the initial model, while in the rotating D1 to D3 models the error is $\leq 4 \%$. In all these cases the above formula overestimates the black hole mass. Because of its rapid rotation and the resulting strong centrifugal forces, in model D4 the collapse deviates significantly from sphericity, leading to a strongly oblate form of the density stratification. Consequently, we still find a non-negligible amount of matter outside the apparent horizon at the end of the simulation (about $12 \%$ of the total rest mass). Therefore the value for $M_{\mathrm{BH}}$ is $8.2 \%$ smaller than $M_{\mathrm{ADM}}$. In Fig. 4 we present the distribution of the rest-mass density and the location of the apparent horizon at the end of the simulation for this particular model. Since the time evolution is limited by our chosen, still computationally affordable, grid resolution in the central region, we are not able to evolve this model to times when a disk forms as in [47]. Nevertheless, all other quantities qualitatively agree with the results in that work, although we refrain from performing a more detailed comparison due to the respective differences in the gauge of the two formulations used in [47] and in this study, respectively.

![](https://cdn.mathpix.com/cropped/2025_06_02_4610231bf4be6dd7bda1g-11.jpg?height=844&width=860&top_left_y=158&top_left_x=1093)


FIG. 4. Isocontours of the rest-mass density for model D4 after the apparent horizon first appears at $t=129.9$. The dashed line shows the location of the apparent horizon.
In the near future we plan to carry out an exhaustive analysis of the scenario of a collapse to a black hole by comparing, on one hand, the CFC formulation with the FCF (see Sec. V) and, on the other hand, the FCF with other (free evolution) formulations. The difficulties induced by the use of different gauges can be overcome by using gauge-invariant quantities for comparison and analyzing their behavior as a function of proper time.
\section*{V. GENERALIZATION TO THE FULLY CONSTRAINED FORMALISM}
The ideas presented in Sec. III can be generalized to the FCF approach of the full Einstein equations described in Sec. II A.
As shown in [24], the hyperbolic part of the FCF can be split into a first-order system. The reformulation of the CFC equations presented in Sec. III relies on the rescaled extrinsic curvature $\hat{A}^{i j}$ given by Eq. (27). Consequently, we write the FCF hyperbolic part as a first-order system in ( $h^{i j}, \hat{A}^{i j}$ ), instead of a first-order system in ( $h^{i j}, \partial h^{i j} / \partial t$ ) as in [24], arriving at
$$
\begin{align}
\frac{\partial h^{i j}}{\partial t}= & 2 N \psi^{-6} \hat{A}^{i j}+\beta^{k} w_{k}^{i j}-\tilde{\gamma}^{i k} \mathcal{D}{k} \beta^{j}-\tilde{\gamma}^{k j} \mathcal{D}{k} \beta^{i} \
& +\frac{2}{3} \tilde{\gamma}^{i j} \mathcal{D}{k} \beta^{k} \tag{34}
\end{align}
$$
$$
\begin{align}
\frac{\partial \hat{A}^{i j}}{\partial t}= & -\mathcal{D}{k}\left(-\frac{N \psi^{2}}{2} \tilde{\gamma}^{k l} w_{l}^{i j}-\beta^{k} \hat{A}^{i j}\right)-\hat{A}^{k j} \mathcal{D}{k} \beta^{i}-\hat{A}^{i k} \mathcal{D}{k} \beta^{j}+\frac{2}{3} \hat{A}^{i j} \mathcal{D}{k} \beta^{k}+2 N \psi^{-6} \tilde{\gamma}{k l} \hat{A}^{i k} \hat{A}^{j l} \
& -8 \pi N \psi^{6}\left(\psi^{4} S^{i j}-\frac{S \tilde{\gamma}^{i j}}{3}\right)+N\left(\psi^{2} \tilde{R}{}^{i j}+8 \tilde{\gamma}^{i k} \tilde{\gamma}^{j l} \mathcal{D}{k} \psi \mathcal{D}{l} \psi\right)+4 \psi\left(\tilde{\gamma}^{i k} \tilde{\gamma}^{j l} \mathcal{D}{k} \psi D_{l} N+\tilde{\gamma}^{i k} \tilde{\gamma}^{j l} \mathcal{D}{k} N \mathcal{D}{l} \psi\right) \
& -\frac{1}{3}\left[N\left(\psi^{2} \tilde{R}+8 \tilde{\gamma}^{k l} \mathcal{D}{k} \psi \mathcal{D}{l} \psi\right)+8 \psi \tilde{\gamma}^{k l} \mathcal{D}{k} \psi \mathcal{D}{l} N\right] \tilde{\gamma}^{i j}-\frac{1}{2}\left(\tilde{\gamma}^{i k} w_{k}^{l j}+\tilde{\gamma}^{k j} w_{k}^{i l}\right) \mathcal{D}{l}\left(N \psi^{2}\right)-\tilde{\gamma}^{i k} \tilde{\gamma}^{j l} \mathcal{D}{k} \mathcal{D}{l}\left(N \psi^{2}\right) \
& +\frac{1}{3} \tilde{\gamma}^{i j} \tilde{\gamma}^{k l} \mathcal{D}{k} \mathcal{D}{l}\left(N \psi^{2}\right) \tag{35}
\end{align}
$$
where
$$
\begin{gather}
w{k}^{i j}:=\mathcal{D}{k} h^{i j}, \tag{36}\
\tilde{R}{}^{i j}:=\frac{1}{2}\left[-w_{l}^{i k} w_{k}^{j l}-\tilde{\gamma}{k l} \tilde{\gamma}^{m n} w{m}^{i k} w_{n}^{j l}+\tilde{\gamma}{n l} w{k}^{m n}\left(\tilde{\gamma}^{i k} w_{m}^{j l}+\tilde{\gamma}^{j k} w_{m}^{i l}\right)\right]+\frac{1}{4} \tilde{\gamma}^{i k} \tilde{\gamma}^{j l} w_{k}^{m n} \mathcal{D}{l} \tilde{\gamma}{m n} . \tag{37}
\end{gather}
$$
The system is closed by adding the equation
$$
\begin{equation}
\frac{\partial w_{k}^{i j}}{\partial t}-\mathcal{D}{k}\left(\beta^{l} w{l}^{i j}+2 N \psi^{-6} \hat{A}^{i j}\right)=-w_{k}^{i l} \mathcal{D}{l} \beta^{j}-\tilde{\gamma}^{i l} \mathcal{D}{k} \mathcal{D}{l} \beta^{j}-w{k}^{l j} \mathcal{D}{l} \beta^{i}-\tilde{\gamma}^{l j} \mathcal{D}{k} \mathcal{D}{l} \beta^{i}+\frac{2}{3} \tilde{\gamma}^{i j} \mathcal{D}{k} \mathcal{D}{l} \beta^{l}+\frac{2}{3} w{k}^{i j} \mathcal{D}{l} \beta^{l} \tag{38}
\end{equation}
$$
which is derived from applying partial derivatives with respect to $t$ in the definition of $w{k}^{i j}$. Moreover, the system observes the constraint of the Dirac gauge, $w_{i}^{i j}=0$ [Eq. (6)], and for the determinant of the conformal metric, we obtain $\tilde{\gamma}=f$. The first-order system given by Eqs. (34)-(38) has the same properties regarding hyperbolicity and existence of fluxes as the one in [24]. It has the advantage over the second-order system for $h^{i j}$ proposed in Ref. [23] of getting rid of partial derivatives with respect to $t$ of the lapse $N$, the shift $\beta^{i}$, or the conformal factor $\psi$.
The elliptic part of the FCF can be rewritten, using the tensor $\hat{A}^{i j}$, as
$$
\begin{align}
& \tilde{\gamma}^{k l} \mathcal{D}{k} \mathcal{D}{l} \psi=-2 \pi \psi^{-1} E^{}-\frac{\tilde{\gamma}{i l} \tilde{\gamma}{j m} \hat{A}^{l m} \hat{A}^{i j}}{8 \psi^{7}}+\frac{\psi \tilde{R}}{8} \tag{39}\
& \tilde{\gamma}^{k l} \mathcal{D}{k} \mathcal{D}{l}(N \psi)= {\left2 \pi \psi^{-2}\left(E^{}+2 S^{}\right)\right.} \
&+\left.\left(\frac{7 \tilde{\gamma}{i l} \tilde{\gamma}{j m} \hat{A}^{l m} \hat{A}^{i j}}{8 \psi^{8}}+\frac{\tilde{R}}{8}\right)\right \tag{40}\
& \tilde{\gamma}^{k l} \mathcal{D}{k} \mathcal{D}{l} \beta^{i}+\frac{1}{3} \tilde{\gamma}^{i k} \mathcal{D}{k} \mathcal{D}{l} \beta^{l} \
&=16 \pi N \psi^{-6} \tilde{\gamma}^{i j}\left(S^{}\right){j}+\hat{A}^{i j} \mathcal{D}{j}\left(2 N \psi^{-6}\right)-2 N \psi^{-6} \Delta_{k l}^{i} \hat{A}^{k l} \tag{41}
\end{align}
$$
The strategy to evolve the two symmetric tensors $h^{i j}$ and $\hat{A}^{i j}$ relies on a decomposition of these tensors in longitudinal and transverse traceless parts. The longitudinal parts (divergences with respect to the flat metric) are either known a priori or are determined by the elliptic equations. More specifically, the divergence of $h^{i j}$ vanishes according
to the Dirac gauge, whereas the divergence of $\hat{A}^{i j}$ is determined by the momentum constraint (42)-see below. Consequently, focus is placed on the transverse traceless parts of these tensors. The latter are described in a purespin tensor harmonic decomposition, as discussed in a previous article [24]. In particular, each transverse traceless tensor is fully expressed in terms of two scalar potentials (named $A$ and $\tilde{B}$ in [24]) that are evolved according to evolution equations obtained from the transverse traceless parts of Eqs. (34) and (35) for $h^{i j}$ and $\hat{A}^{i j}$, respectively, by consistently applying the decomposition in [24]. Once the scalar potentials on the next time slice are determined, the tensors $h^{i j}$ and $\hat{A}_{\text {TT }}^{i j}$ can be reconstructed completely, satisfying the divergence-free conditions. This fully fixes $h^{i j}$, whereas in the case of $\hat{A}^{i j}$ the longitudinal part is computed in a very similar way to the CFC case, i.e. by determining the vector $X^{i}$ from the momentum constraint as described hereafter.
From Eq. (26), the momentum constraint can be written as
$$
\begin{equation}
\mathcal{D}{j} \hat{A}^{i j}=8 \pi \tilde{\gamma}^{i j}\left(S^{}\right){j}-\Delta_{k l}^{i} \hat{A}^{k l} \tag{42}
\end{equation}
$$
which is equivalent to the following elliptic equation for $X^{i}$ :
$$
\begin{align}
& \mathcal{D}{j} \mathcal{D}^{j} X^{i}+\frac{1}{3} \mathcal{D}^{i} \mathcal{D}{k} X^{k}+\tilde{\gamma}^{i m}\left(\mathcal{D}{k} \tilde{\gamma}{m l}-\frac{\mathcal{D}{m} \tilde{\gamma}{k l}}{2}\right) \
& \quad \times\left(\mathcal{D}^{k} X^{l}+\mathcal{D}^{l} X^{k}-\frac{2}{3} f^{k l} \mathcal{D}{p} X^{p}\right) \
& =8 \pi \tilde{\gamma}^{i j}\left(S^{}\right){j}-\tilde{\gamma}^{i m}\left(\mathcal{D}{k} \tilde{\gamma}{m l}-\frac{\mathcal{D}{m} \tilde{\gamma}{k l}}{2}\right) \hat{A}_{\mathrm{TT}}^{k l} \tag{43}
\end{align}
$$
This elliptic equation for the vector $X^{i}$ is linear. Since $h^{i j}$ and $\hat{A}{\text {TT }}^{i j}$ have been calculated previously, we can solve the elliptic equation (43) to obtain the vector $X^{i}$. With this method, the Dirac gauge and the momentum constraint are guaranteed to be satisfied. Then, $\hat{A}^{i j}$ is reconstructed from $\hat{A}{\mathrm{TT}}^{i j}$ and $X^{i}$ on the new time slice.
At this point, since the tensors $h^{i j}$ and $\hat{A}^{i j}$ are known, we can follow exactly the same scheme as in the CFC case to solve in a hierarchical way the elliptic equations. First the conformal factor is obtained from Eq. (39), then the lapse function is obtained from Eq. (40), and finally the shift vector is acquired from Eq. (41). These equations are decoupled in the order mentioned. No sign problems are exhibited in the scalar elliptic equation, and therefore the maximum principle can be applied. A minor concern is associated with the sign of the term $\tilde{R}$ in Eq. (39), but unique solutions also exist for negative conformal Ricci scalars (closely related to $\tilde{R}$ ). Note that, contrary to the CFC case, here no (additional) approximation has been made: it is simply a new scheme to write down the FCF, where the elliptic part is better behaved from the point of view of local uniqueness. Numerical simulations with this FCF scheme will be presented in a future publication.
\section*{VI. DISCUSSION}
\section*{A. Summary}
We have presented an approach to the solution of the uniqueness issues appearing in certain constrained formulations of Einstein equations. We have illustrated the problem and its solution through a detailed analytical and numerical study of a waveless approximation that retains all the involved essential features.
More specifically, we have reformulated XCTS-like elliptic systems appearing in constrained evolution schemes of the Einstein equations, like the FCF of $[23,24]$, as well as in the CFC approximation [4,5]. Such systems require the simultaneous solution of the constraints, in particular, the momentum constraint for the shift, together with a maximal slicing condition for the lapse. The resulting elliptic system presents potential local nonuniqueness problems, and numerical implementations have indeed encountered such obstacles. The original CFC formulation has not been able to cope with these problems, as it suffers from convergence of the system to unphysical solutions or from nonconvergence in high density regimes. We have suggested that these problems are not due to the approximative nature of the CFC, since FCF in the variant of [23,24], which is a natural generalization of the CFC to the nonconformally flat case, also suffers from the same problems. In order to address these issues, first focusing on the simpler CFC case, we have considered the conformal rescaling of the traceless part of the extrinsic curvature, resulting in the expression for $\hat{A}^{i j}$ in Eq. (27), which is a rescaling that is different from the respective ones em-
ployed in the FCF and the CFC approximation, but coincides with the one in the XCTS approach of $[10,11]$. This is motivated by the work of Shapiro and Teukolsky [43], who simulated the collapse of a neutron star model using such a reformulation of the CFC metric equations (however, restricted to spherical symmetry in their case) and apparently did not encounter any of the problems described above. Extending their approach to three dimensions, we have decomposed $\hat{A}^{i j}$ into longitudinal and transverse parts as in the CTT formulation of the constraint equations (29). The divergence (i.e. the longitudinal part) of this tensor is determined by the momentum constraints, Eqs. (28) in the CFC case, just as in the CTT formulation. In the CFC scheme, we have neglected the transverse part of this tensor, as the order of its error is higher than the one arising from the CFC approximation itself. In the nonapproximate FCF case, the transverse part of $\hat{A}^{i j}$ is determined by an evolution equation. Once the conformal extrinsic curvature is obtained, it can be employed in the Hamiltonian equation to calculate the conformal factor $\psi$. The lapse is then fixed through the maximal slicing condition, and the resulting equation allows the application of a maximumprinciple uniqueness argument. Finally, the shift is found through the kinematical relationship defining the extrinsic curvature, leading to Eq. (33).
By performing a variety of tests, we have provided evidence that the problem of convergence to an unphysical solution of the metric equations (or even complete nonconvergence) in the original formulation of the CFC scheme is fully cured by our new reformulation. Not only can numerical results in the original CFC scheme (in the, at most, moderately gravitationally compact regime where that system still yields physically correct solutions) be reproduced by the new formulation but, more importantly, the new numerical results presented here exhibit the proper numerical and physical behavior even for highly compact configurations. For the first time, it has been possible to successfully perform both the migration test and the collapse of a neutron star to a black hole in the CFC case in a consistent way. Our new formulation thus facilitates simulations in the high density regime of those scenarios where the CFC is still a reasonably fair approximation, that is, for systems which are not too far from sphericity, like stellar gravitational collapse.
\section*{B. Comparison with previous works}
As compared to the original CFC formulation by Isenberg [4] and Mathews and Wilson [5], the scheme presented here is augmented by an additional vector elliptic equation for $X^{i}$, while the elliptic character of the system of metric equations is preserved. The new scheme reformulates the CFC approximation in a CTT shape (one scalar and one vector elliptic equation), and then solves for the lapse and the shift (one additional scalar and one vector elliptic equation). In contrast, the original CFC scheme
employed a (X)CTS approach where, together with two scalar elliptic equations, only one vector elliptic equation was present. In contrast to the original scheme, the elliptic system in the new formulation not only corrects the problem of local uniqueness in the scalar elliptic equations, but also introduces a hierarchical structure that decouples the system in one direction.
In the context of the conformally flat approximation, the same "augmented CFC" scheme as that discussed here has been introduced already by Saijo [8] to compute the gravitational collapse of differentially rotating supermassive stars. However, in this work the inconsistency between Eqs. (25) and (29), i.e. setting to zero the transverse traceless part of $\hat{A}^{i j}$, has not been pointed out. On the contrary, we have analyzed this inconsistency in detail (cf. the Appendix) and have shown that it leads to an error of the same order as that of the CFC approximation. In addition, we have shown here that the introduction of the vector potential $X^{i}$ is the key ingredient for solving the nonuniqueness issue.
The same scheme, but without the conformal rescaling of the matter quantities, has also been used recently by Shibata and Uryū [45] in the context of computing initial data. As in [8], the inconsistency resulting from setting to zero the transverse traceless part of $\hat{A}^{i j}$ and the uniqueness issue are not discussed in their work. We emphasize that these studies [8,45] do not discuss the extension of the new scheme to the nonconformally flat case, as done here.
Let us also mention that the augmented CFC scheme presented here can be regarded as a hybrid mixture of some of the waveless approximation theories (WAT) proposed by Isenberg [4]. In fact, the CFC approximation using the two choices $\tilde{\gamma}{i j}=f{i j}$ and $\partial_{t} \tilde{\gamma}{i j}=0$ [as employed in Eq. (33)] corresponds to WAT-I. On the other hand, the approximation $\hat{A}{\text {TT }}^{i j}=0$ used in Eq. (29) is in the spirit of the vanishing transverse traceless part of the extrinsic curvature in the (coupled) version, WAT-II (although WAT-II refers to the physical extrinsic curvature, whereas here we have dealt with the conformal one). As mentioned above, both assumptions are consistent at the considered level of approximation, as shown in the Appendix.
Regarding the complete constrained evolution of the Einstein equations, we have generalized the ideas presented here for the CFC case to the elliptic part of the FCF. In previous studies [23,24], the hyperbolic part of the Einstein equations resulted in a wave-type equation for the tensor $h^{i j}$, representing the deviation of the three-metric from conformal flatness. With the introduction of $\hat{A}^{i j}$ we have recovered here a first-order evolution system, analogous to the standard Hamiltonian $3+1$ system, in which we have, however, retained only the divergence-free terms. Thus, for both $h^{i j}$ and $\hat{A}^{i j}$, the transverse (divergence-free) parts are evolved by this system, while the longitudinal parts are fixed either by the gauge (for $h^{i j}$ ), or by the
momentum constraint (for $\hat{A}^{i j}$ ). Numerical results for this case will be presented in future studies.
We finally comment on the recent work by Rinne [18], where uniqueness problems appearing in certain constrained and partially constrained schemes for vacuum axisymmetric Einstein equations [20,54] are addressed. As in the present case, uniqueness issues related to the Hamiltonian constraint equation are solved by adopting an appropriate rescaling of the extrinsic curvature. On the other hand, problems associated with the slicing condition are tracked to the substitution in that equation of the extrinsic curvature by its kinematical expression in terms of the (shift and the) lapse. The latter spoils the uniqueness properties by reversing the sign of the relevant term in the slicing equation. This problem is solved by enlarging the elliptic system with an additional vector so as to reexpress the relevant components of the extrinsic curvature without resorting to the lapse. The resulting elliptic system also presents a hierarchical structure. Although the spirit of such an approach is close to the one presented here, the specific manner of introducing the additional vector variable in [18] critically relies on the two dimensionality of the axisymmetric problem (specifically, on a choice of a particular gauge and on the fact that vectors and rank-two traceless symmetric tensors have the same number of components in two dimensions, a property lost in three dimensions). On the contrary, the introduction of the vector $X^{i}$ through the CTT decomposition (29) is properly devised to work in three dimensions. Relevant discussions in the three-dimensional context can be found in Sec. 3.4 of [18] (where the relation between nonuniqueness problems in XCTS and axisymmetric constrained evolution schemes is discussed) and in the threedimensional constrained evolution scheme presented by Moncrief et al. in [25].
\section*{ACKNOWLEDGMENTS}
It is a pleasure to thank J. M. Ibáñez, M. Saijo, and K. Uryū for many fruitful discussions. I. C.-C. acknowledges support from the Spanish Ministerio de Educación y Ciencia (AP2005-2857), and H. D. is supported by the Marie Curie Intra-European Fellowship within the 6th European Community Framework Programme (IEF 040464). This work was supported by the Ministerio de Educación y Ciencia through Grant No. AYA2007-67626-C03-01, by the Deutsche Forschungsgemeinschaft through the Transregional Collaborative Research Center Grant No. SFB/TR 7 "Gravitational Wave Astronomy," by the French ANR Grant No. 06-2-134423 "Méthodes mathématiques pour la relativité générale," by the FrenchSpanish bilateral research Grant No. HF2005-0115, and by the European Network of Theoretical Astroparticle Physics ENTApP ILIAS/N6 under Contract No. RII3-CT-2004-506222.
\section*{APPENDIX: CONSISTENCY OF THE APPROXIMATION}
In the derivation of the new formalism, we make use of the fact that $(L X)^{i j} \approx \hat{A}^{i j}$ in the CFC. We show next that this assumption is completely consistent at the accuracy level of the CFC approximation. In the first place, we need to estimate the error of the CFC approximation itself. By definition, the CFC three-metric deviates linearly with $h^{i j}$ from the (exact) FCF case. It can easily be shown from the FCF equations (39)-(41) that the metric quantities behave as
$$
\begin{align}
& \psi=\psi_{\mathrm{CFC}}+\mathcal{O}(h), \tag{A1}\
& N=N_{\mathrm{CFC}}+\mathcal{O}(h), \tag{A2}\
& \beta^{i}=\beta_{\mathrm{CFC}}^{i}+\mathcal{O}(h) . \tag{A3}
\end{align}
$$
Therefore $h^{i j}$ can be used as an estimator for the error of the CFC approximation.
Two limits in which the CFC is exact will be considered. First, in spherical symmetry the CFC metric system is an exact reformulation of the Einstein equations since $h^{i j}=0$ in the FCF metric. If the system is close to spherical symmetry (i.e. spheroidal), and if we are able to define a quasispherical surface of the system (e.g., the surface of a star or the apparent horizon of a black hole), then the equatorial and polar circumferential proper radius, $R_{\mathrm{e}}$ and $R_{\mathrm{p}}$, can be computed, and we can define the ellipticity of the system as
$$
\begin{equation}
e^{2}:=1-R_{\mathrm{p}}^{2} / R_{\mathrm{e}}^{2} \tag{A4}
\end{equation}
$$
Close to sphericity, $e^{2}$ scales linearly with $h^{i j}$, and we can ensure that the error of the CFC is $h^{i j} \sim \mathcal{O}\left(e^{2}\right)$. The second limit to consider is if a post-Newtonian expansion of the gravitational sources is possible, i.e. if the post-Newtonian parameter $\max \left(v^{2} / c^{2}, G M / L c^{2}\right)<1$, where $v, M$, and $L$ are the typical velocity, mass, and length of the system, respectively. In this case the CFC metric behaves like the first post-Newtonian approximation [55,56], i.e.
$$
\begin{gather}
\psi=\psi_{\mathrm{CFC}}+\mathcal{O}\left(1 / c^{4}\right), \tag{A5}\
N=N_{\mathrm{CFC}}+\mathcal{O}\left(1 / c^{4}\right), \tag{A6}\
c \beta^{i}=c \beta_{\mathrm{CFC}}^{i}+\mathcal{O}\left(1 / c^{4}\right) . \tag{A7}
\end{gather}
$$
Note that, for clarity, we explicitly retain powers of the speed of light $c$ as factors in the equations throughout this appendix. In the case that both limits are valid, i.e. close to sphericity and in the post-Newtonian expansion, the nonconformally-flat part of the three-metric behaves like $h^{i j} \sim \mathcal{O}\left(e^{2} / c^{4}\right)$. The next step is to compute the behavior of the CFC metric if we assume $(L X)^{i j} \approx \hat{A}^{i j}$, considering the two limiting cases introduced above.
In the spherically symmetric case the relation $(L X)^{i j}=$ $\hat{A}^{i j}$ is trivially fulfilled. Therefore the behavior for a quasispherical configuration is also $h^{i j} \sim \mathcal{O}\left(e^{2}\right)$ even if $\hat{A}_{\text {TT }}^{i j}=0$ is assumed. This limit in the approximation is very important, since it is independent of the strength of the gravitational field. For example, it allows us to evolve black holes, with the only condition being that $h^{i j}$ should be small, i.e. close to the sphericity.
To check the approximation in the post-Newtonian limit, we need to compare $\beta_{\mathrm{CFC}}^{i}$ and $X^{i}$. This can be done by means of the post-Newtonian expansion of the sources of Eqs. (17) and (30), respectively,
$$
\begin{gather}
\Delta \beta_{\mathrm{CFC}}^{i}+\frac{1}{3} \mathcal{D}^{i} \mathcal{D}{j} \beta{\mathrm{CFC}}^{j}=16 \pi S^{ i}+\mathcal{O}\left(1 / c^{7}\right) \tag{A8}\
\Delta X^{i}+\frac{1}{3} \mathcal{D}^{i} \mathcal{D}_{j} X^{j}=8 \pi S^{* i}+\mathcal{O}\left(1 / c^{7}\right) \tag{A9}
\end{gather*}
$$
From the comparison of Eqs. (A8) and (A9) we obtain that
$$
\begin{equation}
c^{3} \frac{\beta_{\mathrm{CFC}}^{i}}{2}=c^{3} X^{i}+\mathcal{O}\left(1 / c^{2}\right) \tag{A10}
\end{equation}
$$
Thus $\hat{A}^{i j}$ can be computed in terms of $X^{i}$ as
$$
\begin{equation}
c^{4} \hat{A}^{i j}=\frac{\psi_{\mathrm{CFC}}^{6}}{2 N_{\mathrm{CFC}}} c^{4}\left(L \beta_{\mathrm{CFC}}\right)^{i j}=c^{4}(L X)^{i j}+\mathcal{O}\left(1 / c^{2}\right) \tag{A11}
\end{equation}
$$
where we make use of $\psi_{\mathrm{CFC}}^{6} / N_{\mathrm{CFC}}=1+\mathcal{O}\left(1 / c^{2}\right)$. The effect of using $(L X)^{i j}$ instead of $\hat{A}^{i j}$ in the calculation of the CFC metric can be seen in the expressions
$$
\begin{align}
\psi_{\mathrm{CFC}} & =\Delta_{\mathrm{s}}^{-1} \mathcal{S}{(\psi)}\left(N{\mathrm{CFC}}, \psi_{\mathrm{CFC}}, \hat{A}^{i j}\right) \
& =\Delta_{\mathrm{s}}^{-1} \mathcal{S}{(\psi)}\left(N{\mathrm{CFC}}, \psi_{\mathrm{CFC}},(L X)^{i j}\right)+\mathcal{O}\left(1 / c^{8}\right) \tag{A12}
\end{align}
$$
$$
\begin{align}
N_{\mathrm{CFC}} & =\psi_{\mathrm{CFC}}^{-1} \Delta_{\mathrm{s}}^{-1} \mathcal{S}{(N \psi)}\left(N{\mathrm{CFC}}, \psi_{\mathrm{CFC}}, \hat{A}^{i j}\right) \
& =\psi_{\mathrm{CFC}}^{-1} \Delta_{\mathrm{s}}^{-1} \mathcal{S}{(N \psi)}\left(N{\mathrm{CFC}}, \psi_{\mathrm{CFC}},(L X)^{i j}\right)+\mathcal{O}\left(1 / c^{8}\right), \tag{A13}
\end{align}
$$
$$
\begin{align}
c \beta_{\mathrm{CFC}}^{i} & =c \Delta_{\mathrm{V}}^{-1} \mathcal{S}{(\beta)}\left(N{\mathrm{CFC}}, \psi_{\mathrm{CFC}}, \hat{A}^{i j}\right) \
& =c \Delta_{\mathrm{V}}^{-1} \mathcal{S}{(\beta)}\left(N{\mathrm{CFC}}, \psi_{\mathrm{CFC}},(L X)^{i j}\right)+\mathcal{O}\left(1 / c^{6}\right) \tag{A14}
\end{align}
$$
where $\mathcal{S}{(\psi)}, \mathcal{S}{(N \psi)}$, and $\mathcal{S}{(\beta)}$ are the sources of Eqs. (31)(33), and $\Delta{\mathrm{s}}^{-1}$ and $\Delta_{\mathrm{v}}^{-1}$ are just the inverse operators appearing in the right-hand sides of these equations (for the scalars $\psi$ and $N \psi$, and for the vector $\beta^{i}$, respectively). When comparing Eqs. (A12)-(A14) with Eqs. (A5)-(A7), it becomes obvious that in all cases the error introduced by making the approximation $(L X)^{i j} \approx \hat{A}^{i j}$ is smaller than the error of the CFC approximation itself.
As an illustration of the above properties, we study the influence of the $\hat{A}_{\text {TT }}^{i j}$ term in Eq. (29) when computing rotating neutron star models with a polytropic $\Gamma=2$ equa-

![](https://cdn.mathpix.com/cropped/2025_06_02_4610231bf4be6dd7bda1g-16.jpg?height=1208&width=874&top_left_y=155&top_left_x=173)


FIG. 5. Consistency of the approximation for rotating neutron star models. In the top panel max $\left|\hat{A}{\text {TT }}^{i j} / \hat{A}^{i j}\right|$ for the FCF (solid line) and the CFC (dashed line) as well as the maximum deviation from conformal flatness max $\left|h^{i j}\right|$ for the FCF (dashdotted line) are plotted against the ellipticity $e$. The bottom panel shows the absolute difference $\left|N{\mathrm{c}, \mathrm{CFC}}-N_{\mathrm{c}}\right|$ in the central value of the lapse between the CFC and FCF (solid line) and the absolute difference $\left|N_{\mathrm{c}, \mathrm{CFC}}-N_{\mathrm{c}, \mathrm{CFC}}\right|$ between the regular CFC and the CFC neglecting $\hat{A}{\text {TT }}^{i j}$ in Eq. (29) (dashed line). The Kepler limit is marked by vertical dotted lines, while the slanted dotted lines represent the order of accuracy with respect to powers of $e$.
tion of state. This model setup contains the initial models used in Sec. IV. They assume axial symmetry and stationarity, in combination with rigid rotation. We build a sequence of rotating polytropes with increasing rotational frequencies, while keeping the central enthalpy fixed,
which produces models of increasing masses from $M=$ $1.33 M{\odot}$ (no rotation) to $M=1.57 M_{\odot}$ (the Kepler limit; see below). For all these models, we use three gravitational field schemes: the exact Einstein equations using the stationary ansatz in the FCF, and the two approximate ones, regular CFC and CFC, neglecting the term $\hat{A}{\text {TT }}^{i j}$ in Eq. (29). The results are displayed on a logarithmic scale in Fig. 5. In the top panel we show the maximal amplitudes of $\hat{A}{\text {TT }}^{i j}$ (relatively, to $\hat{A}^{i j}$ ) in both the FCF and the regular CFC, as functions of the ellipticity $e$ defined in Eq. (A4). This quantity is physically and numerically limited by the minimal rotational period at the so-called mass-shedding limit (or Kepler limit), when centrifugal forces exactly balance gravitational and pressure forces at the star's equator. In the FCF case we plot the maximal amplitude of $h^{i j}$. This quantity is dimensionless and represents the deviation of the three-metric from conformal flatness, which can be interpreted as the relative error one makes in the metric when using the CFC instead of the FCF. Note that this error in computing $\hat{A}^{i j}$ by discarding the $\hat{A}{\text {TT }}^{i j}$ term in the CFC approximation is roughly of the same magnitude as the error on the metric in the CFC case. All these quantities decrease like $\mathcal{O}\left(e^{2}\right)$ as expected, except for stars rotating close to the Kepler limit. Indeed, the development in powers $e$ is equivalent to a slow-rotation approximation (see, e.g., [57]) by perturbing spherically symmetric configurations, and, when comparing these slow-rotation results to numerical "exact" ones for rigidly rotating stars (see, e.g., [58] in the two-fluids case), one sees that they usually agree extremely well, excepted very close to the Kepler limit, where this "perturbed spherical symmetry" approach is no longer valid. Finally, because $\hat{A}^{i j}$ appears as a quadratic source term in the Poisson-like equations (15) and (16), the overall errors on the lapse $N$ or the conformal factor $\psi$ are even smaller, as shown in the bottom panel of Fig. 5. In the case of the central value $N{\mathrm{c}}$ of the lapse, the error due to the CFC approximation is maximal at the Kepler limit and $\lesssim 10^{-4}$ for the studied sequence. The error which is then due to neglecting $\hat{A}_{\text {TT }}^{i j}$ within the CFC scheme amounts to $\lesssim 10^{-6}$ and decreases faster than the error due to the CFC approximation, namely, as $\mathcal{O}\left(e^{4}\right)$, again except near the Kepler limit. Our tests thus show that for stationary rotating neutron star models this additional approximation induces an error which falls within the overall CFC approximation.
[1] J.W. York, Sources of Gravitational Radiation (Cambridge University Press, Cambridge, England, 1979), p. 83.
[2] M. Alcubierre, Introduction to $3+1$ Numerical Relativity (Oxford University Press, Oxford, England, 2008).
[3] E. Gourgoulhon, arXiv:gr-qc/0703035.
[4] J. A. Isenberg, Int. J. Mod. Phys. D 17, 265 (2008).
[5] J. R. Wilson and G. J. Mathews, in Frontiers in Numerical Relativity, edited by C.R. Evans, L. S. Finn, and D. W. Hobill (Cambridge University Press, Cambridge, England, 1989), p. 306.
[6] H. Dimmelmeier, J. A. Font, and E. Müller, Astron. Astrophys. 393, 523 (2002).
[7] R. Oechslin, S. Rosswog, and F.-K. Thielemann, Phys. Rev. D 65, 103005 (2002).
[8] M. Saijo, Astrophys. J. 615, 866 (2004).
[9] E. B. Abdikamalov, H. Dimmelmeier, L. Rezzolla, and J. C. Miller (unpublished).
[10] J. W. York Jr., Phys. Rev. Lett. 82, 1350 (1999).
[11] H. P. Pfeiffer and J. W. York Jr., Phys. Rev. D 67, 044022 (2003).
[12] H. P. Pfeiffer, J. Hyperbol. Diff. Equat. 2, 497 (2005).
[13] H. P. Pfeiffer and J. W. York Jr., Phys. Rev. Lett. 95, 091101 (2005).
[14] T. W. Baumgarte, N. O'Murchadha, and H. P. Pfeiffer, Phys. Rev. D 75, 044009 (2007).
[15] D. Walsh, Classical Quantum Gravity 24, 1911 (2007).
[16] N. O'Murchadha and J. W. York Jr., J. Math. Phys. (N.Y.) 14, 1551 (1973).
[17] N. O'Murchadha and J. W. York Jr., Phys. Rev. D 10, 428 (1974).
[18] O. Rinne, Classical Quantum Gravity 25, 135009 (2008).
[19] J. L. Jaramillo, J. A. Valiente Kroon, and E. Gourgoulhon, Classical Quantum Gravity 25, 093001 (2008).
[20] M. W. Choptuik, E. W. Hirschmann, S. L. Liebling, and F. Pretorius, Classical Quantum Gravity 20, 1857 (2003).
[21] O. Rinne, Ph.D. thesis, University of Cambridge, 2005 [arXiv:gr-qc/0601064].
[22] O. Rinne and J. M. Stewart, Classical Quantum Gravity 22, 1143 (2005).
[23] S. Bonazzola, E. Gourgoulhon, P. Grandclément, and J. Novak, Phys. Rev. D 70, 104007 (2004).
[24] I. Cordero-Carrión, J. M. Ibáñez, E. Gourgoulhon, J. L. Jaramillo, and J. Novak, Phys. Rev. D 77, 084007 (2008).
[25] V. Moncrief, L. Buchman, H. P. Pfeiffer, O. Rinne, and O. Sarbach, in From Geometry to Numerics, Proceedings of Institut Henri Poincaré, Paris, 2006, www.luth.obspm.fr/ IHP06/workshops/geomnum/; V. Moncrief (private communication).
[26] J. R. Wilson, G. J. Mathews, and P. Marronetti, Phys. Rev. D 54, 1317 (1996).
[27] J. L. Jaramillo, M. Ansorg, and F. Limousin, Phys. Rev. D 75, 024019 (2007).
[28] H. Dimmelmeier, J. Novak, J. A. Font, J. M. Ibáñez, and E. Müller, Phys. Rev. D 71, 064023 (2005).
[29] C. D. Ott, H. Dimmelmeier, A. Marek, H.-T. Janka, I. Hawke, B. Zink, and E. Schnetter, Phys. Rev. Lett. 98, 261101 (2007).
[30] P. Cerdá-Durán, J. A. Font, and H. Dimmelmeier, Astron. Astrophys. 474, 169 (2007).
[31] G. B. Cook, S. L. Shapiro, and S. A. Teukolsky, Phys. Rev. D 53, 5533 (1996).
[32] H. Dimmelmeier, N. Stergioulas, and J. A. Font, Mon. Not. R. Astron. Soc. 368, 1609 (2006).
[33] R. Oechslin, H.-T. Janka, and A. Marek, Astron. Astrophys. 467, 395 (2007).
[34] J. A. Faber, P. Grandclément, and F. A. Rasio, Phys. Rev. D 69, 124036 (2004).
[35] F. Banyuls, J. A. Font, J. M. Ibáñez, J. M. Martí, and J. A. Miralles, Astrophys. J. 476, 221 (1997).
[36] J. A. Font, J. Phys. Conf. Ser. 91, 012002 (2007).
[37] M. E. Taylor, Partial Differential Equations III. Nonlinear Equations, Applied Mathematical Sciences Vol. 117 (Springer-Verlag, New York, 1997, corrected reprint of the 1996 original).
[38] L.C. Evans, Partial Differential Equations (American Mathematical Society, Providence, 1998).
[39] M. Protter and H. Weinberger, Maximum Principles in Differential Equations (Prentice-Hall, Englewood Cliffs, 1967).
[40] H. Komatsu, Y. Eriguchi, and I. Hachisu, Mon. Not. R. Astron. Soc. 237, 355 (1989).
[41] N. Stergioulas and J. L. Friedman, Astrophys. J. 444, 306 (1995).
[42] S. Bonazzola, E. Gourgoulhon, M. Salgado, and J. A. Marck, Astron. Astrophys. 278, 421 (1993).
[43] S. L. Shapiro and S. A. Teukolsky, Astrophys. J. 235, 199 (1980).
[44] A. Lichnerowicz, J. Math. Pures Appl. 23, 37 (1944); reprinted in Choix d'œuvres Mathématiques (Hermann, Paris, 1982), p. 4.
[45] M. Shibata and K. Uryū, Phys. Rev. D 74, 121503(R) (2006); Classical Quantum Gravity 24, S125 (2007).
[46] H. Dimmelmeier, J. A. Font, and E. Müller, Astron. Astrophys. 388, 917 (2002).
[47] L. Baiotti, I. Hawke, P. J. Montero, F. Löffler, L. Rezzolla, N. Stergioulas, J. A. Font, and E. Seidel, Phys. Rev. D 71, 024035 (2005).
[48] L.-M. Lin and J. Novak, Classical Quantum Gravity 23, 4545 (2006).
[49] http://www.lorene.obspm.fr/.
[50] L.-M. Lin and J. Novak, Classical Quantum Gravity 24, 2665 (2007).
[51] J. A. Font, T. Goodale, S. Iyer, M. Miller, L. Rezzolla, E. Seidel, N. Stergioulas, W.-M. Suen, and M. Tobias, Phys. Rev. D 65, 084024 (2002).
[52] A. Marek, H. Dimmelmeier, H.-T. Janka, E. Müller, and R. Buras, Astron. Astrophys. 445, 273 (2006).
[53] A. Garat and R. H. Price, Phys. Rev. D 61, 124011 (2000).
[54] D. Garfinkle and G. C. Duncan, Phys. Rev. D 63, 044011 (2001).
[55] P. Cerdá-Durán, G. Faye, H. Dimmelmeier, J. A. Font, J. M. Ibáñez, E. Müller, and G. Schäfer, Astron. Astrophys. 439, 1033 (2005).
[56] W. Kley and G. Schäfer, Phys. Rev. D 60, 027501 (1999).
[57] J. B. Hartle and K. S. Thorne, Astrophys. J. 153, 807 (1968).
[58] R. Prix, J. Novak, and G. L. Comer, Phys. Rev. D 71, 043005 (2005).