# Copyright (C) 2019 Andrew Trettel.  All rights reserved.

flags=-std=f2003 -pedantic -pedantic-errors -nocpp -ffree-form \
-fimplicit-none -Wtabs -fmax-errors=1 -Werror -Wfatal-errors -Wall \
-Wextra -Wconversion-extra

windbag: windbag.f90
	gfortran $(flags) $^ -o $@

.PHONY: clean
clean:
	-rm -frv windbag
