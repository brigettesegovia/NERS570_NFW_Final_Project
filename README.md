# NERS 570 Final Project: Reproducing NFW Profiles Using N-body Simulations

In this repository you will find N-body simulation code used to reproduce a Navarro-Frenk-White density profile. We initialize our simulations using a fine grid particle mesh code, a uniform in space but random in velocity distribution of particles, and evolve these until it meets a criteria of convergence. We adapt the leap-frog integration scheme. 

To fulfill the requirements of this project we implement the following:

1.  __Stationary Iterative Methods__ for solving the discretized Poisson gravitational potential equation
2.  __Parallelization__ of our potential and force equations
3.  __Benchmarking__ the parallel implementation against the serial baseline
4.  __Performing detailed profiling__ to identify its main computational bottlenecks.
5.  TBD



_______
#### Background on NFW Profiles
Navarro, Frenk, and White (NFW) established one of the most influential models for describing the structure of dark matter halos. Using cosmological N-body simulations, they demonstrated that the spherically averaged density profiles of dark matter halos could be described by a universal functional form, independent of halo mass and cosmological parameters, with variations characterized by a set of scaling parameters. This profile describes how the density of dark matter varies with radius, decreasing from the central regions of a halo toward its outskirts [Navarro et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJ...490..493N/abstract). 
The Navarro–Frenk–White density profile is

$$
\rho(r) = \frac{\rho_s}{(r/r_s)(1+r/r_s)^2}
$$

Here, $r_s$ is the scale radius and $\rho_s$ is the characteristic density.

The NFW density profile provides a framework for characterizing the structure and gravitational properties of galaxies and galaxy clusters. In subsequent work, Navarro, Frenk, and White further investigated the structural properties of dark matter halos, including their characteristic densities, virial radii, and circular velocity profiles. The enclosed mass derived from the NFW density profile can be used to calculate the circular velocity as a function of radius, providing a theoretical framework for modeling galaxy rotation curves.
