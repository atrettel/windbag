# Copyright (C) 2019 Andrew Trettel.  All rights reserved.

flags=-std=f2008 -pedantic -pedantic-errors -nocpp -ffree-form \
-fimplicit-none -Wtabs -fmax-errors=1 -Werror -Wfatal-errors -Wall \
-Wextra -Wconversion-extra

windbag: windbag.f90 wbbase.o
	mpifort $(flags) $^ -o $@

wbbase.o: wbbase.f90
	mpifort $(flags) -c $^

.PHONY: clean
clean:
	-rm -frv *.mod
	-rm -frv *.o
	-rm -frv windbag
