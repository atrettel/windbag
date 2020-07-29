Windbag direct numerical simulation code
========================================

A computational fluid dynamics code for the direct numerical simulation of
transitional and turbulent fluid flows.

NOTE: This project is a work in progress and does not work yet.  I have
implemented much of the backend that supports various operations, but I have
yet to implement any of the numerical methods or post-processing operations.
Basic grid generation in parallel does work, though.


Usage
-----

    $ mpirun -n 2 ./windbag [INPUT_FILE]


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
