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

- only MPI for parallelization

- simple Makefile for compilation

- no preprocessor commands or `common` blocks

- minimal number of hard-coded variables and external libraries

- single command line argument (the namelist-based input file)

- input file serves as detailed documentation of the case

- architecture based on series of operations acting on a single set of fields

- output files natively readable by Paraview

- open formats for post-processed output files (CSV)

- arbitrary equations of state and material properties (nothing hard-coded).

- structured, readable log files (redirected from `stdout`).

- detailed documentation covering verification, validation, and uncertainty
  quantification.


-------------------------------------------------------------------------------

Copyright © 2019-2020 Andrew Trettel

