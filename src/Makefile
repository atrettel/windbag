# Copyright (C) 2020 Andrew Trettel
# 
# Windbag is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# Windbag is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# Windbag.  If not, see <https://www.gnu.org/licenses/>.

FC=mpifort
FFLAGS=-std=f2008 -nocpp -fimplicit-none -Wall -O2

project=windbag

OBJ := $(patsubst wb%.f90,wb%.o,$(wildcard wb*.f90))

$(project): $(project).f90 $(OBJ)
	$(FC) $(FFLAGS) $^ -o $@

wb%.o: wb%.f90
	$(FC) $(FFLAGS) -c $^ -o $@

.PHONY: clean
clean:
	-rm -fv $(project)
	-rm -fv *.mod
	-rm -fv *.o
