# dendrite_growth

Three-dimensional random deposition simulation source code inspired by Aryanfar et. al., J. Phys. Chem. Lett. 5, 1721−1726 (2014).
Published work: Inhyuk Jang, Arun Yethiraj*, Effect of diffusion constant on the morphology of dendrite growth in lithium metal batteries, J. Chem. Phys. 154, 234705 (2021). 

- Language: Fortran 90/2008
- Supports OpenMP parallelization

- Compilation commands
  Compile (gfortran): gfortran -O3 -fopenmp -march=native -ffast-math mc.f90 -o mc.x
  Compile (ifx): ifx -O3 -qopenmp -xHost mc.f90 -o mc.x

You need cond.dat for setting simulation parameteres. Other parameters are determined based on the semi-empirical values (see the published work by I. Jang or Aryanfar et. al. above), so not recommended to modify unless you simulate deposition systems other than Lithium ion deposition. 
