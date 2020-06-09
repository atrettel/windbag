Windbag direct numerical simulation code
========================================

A computational fluid dynamics code for the direct numerical simulation of
transitional and turbulent fluid flows.


Usage
-----

    $ mpirun -n 2 ./windbag [INPUT_FILE]


Planned features and goals
--------------------------

- simple, clear, economical code written in modern Fortran

- block structured grid

- high-order finite difference method (FDM)

- compressible Navier-Stokes equations for the governing equations

- only MPI for parallelization (no additional external libraries)

- simple Makefile for compilation

- no preprocessor use

- minimal number of hard-coded variables and global variables

- single command line argument (the namelist-based input file)

- input file serves as detailed documentation of the case

- architecture based on series of operations acting on a single set of fields

- dependency graph to calculate additional fields

- output files natively readable by Paraview

- open formats for post-processed output files (CSV)

- arbitrary equations of state and material properties (nothing hard-coded)

- structured, readable, and "auditable" log files

- detailed documentation covering verification, validation, and uncertainty
  quantification

- designed for DNS of flows with simple geometries but still possible to
  simulate flows with more complex geometries


-------------------------------------------------------------------------------

Copyright © 2019-2020 Andrew Trettel

Windbag is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Windbag is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Windbag.  If not, see <https://www.gnu.org/licenses/>.
