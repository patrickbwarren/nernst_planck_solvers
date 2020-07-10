# Nernst-Planck solvers

##Solvers for Nernst-Planck equations

See "Non-Faradaic electric currents in the Nernst-Planck equations and
non-local diffusiophoresis of suspended colloids in crossed salt
gradients", P. B. Warren, [Phys. Rev. Lett. 124, 248004 (2020)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.248004)

If you use this stuff please cite:
```
@article{PhysRevLett.124.248004,
  title = {Non-Faradaic Electric Currents in the Nernst-Planck Equations and Nonlocal Diffusiophoresis of Suspended Colloids in Crossed Salt Gradients},
  author = {Warren, Patrick B.},
  journal = {Phys. Rev. Lett.},
  volume = {124},
  issue = {24},
  pages = {248004},
  year = {2020},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.124.248004},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.124.248004}
}
```

Basic test case, after `make` do:
```
./crossgrad_creator.py --ionx=H1Cl1 --iony=K1Br1 --nside=100 --nsave=8 --nstep=800 --steady --dt=0.025 testcase --run
./crossgrad_viewer.py --rc=testcase.rc --show
```

Further examples in the shell scripts.

## Inventory

`Makefile` - default target builds drivers.

`crossgrad_mod.f90` - core numerical routines in FORTRAN 90.

`crossgrad_solver.f90` - FORTRAN 90 driver code, reads run config `.rc` file and runs a complete calculation.

`crossgrad_creator.py` - python script to generate the run config and initial converged solution to the inhomogeneous Poisson equation.

`crossgrad_viewer.py` - reads FORTRAN data files and generates images: to view on the screen use `--show`

`crossgrad_resize.py` - for convenience, resizes the grid using interpolation.

`crossgrad_example.sh` - runs a complete calculation including
generating an image stack (`.png` files) and a movie `.mp4` file.

`crossgrad_more_examples.sh` - more examples (beware some of these take a very long time to run).

### Copying

Code in this repository is licensed under GPL v2:

This program is free software: you can redistribute it and / or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see
<http://www.gnu.org/licenses/>.

### Copyright

Copyright &copy; (2020) Patrick B Warren <patrick.warren@stfc.ac.uk>
