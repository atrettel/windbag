Windbag source directory
========================


Compilation
-----------

To compile, just type `make`.  The program does not have any external
dependencies other than MPI.


Description of source files
---------------------------

- `wb_base.f90`

    - Contains functions and subroutines for the `WB_Subdomain` data structure,
      and all subroutines needed to modify field data.

- `wb_exit.f90`

    - Contains error codes and subroutine to forcibly abort the program.

- `wb_grids.f90`

    - Contains functions to create algebraic structured grids.

- `wb_representation.f90`

    - Contains functions, subroutines, and global parameters related to the
      representation of data in memory.

    - This module is the only module that contains global variables, though
      their values should never change.

- `wb_text.f90`

    - Contains subroutines to create text tables and other human-readable
      output.

- `wb_variables.f90`

    - Contains subroutines to create variable lists with dependency graphs.

- `windbag.f90`

    - Contains the main loop of the program.


-------------------------------------------------------------------------------

Copyright © 2020 Andrew Trettel

Windbag is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Windbag is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Windbag.  If not, see <https://www.gnu.org/licenses/>.
