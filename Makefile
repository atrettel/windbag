# Copyright (C) 2019 Andrew Trettel.  All rights reserved.

FC=mpifort

flags=-std=f2008 -pedantic -pedantic-errors -nocpp -ffree-form \
-fimplicit-none -Wtabs -fmax-errors=1 -Werror -Wfatal-errors -Wall \
-Wextra -Wconversion-extra

windbag: windbag.f90 wbbase.o
	@$(FC) $(flags) $^ -o $@
	@echo "FC $<"

wbbase.o: wbbase.f90
	@$(FC) $(flags) -c $^
	@echo "FC $<"

.PHONY: clean
clean:
	-rm -frv *.mod
	-rm -frv *.o
	-rm -frv windbag
