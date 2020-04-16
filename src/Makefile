# Copyright (C) 2020 Andrew Trettel.

FC=mpifort

project=windbag

obj_files := $(patsubst wb%.f90,wb%.o,$(wildcard wb*.f90))

$(project): $(project).f90 $(obj_files)
	$(FC) $^ -o $@

wb%.o: wb%.f90
	$(FC) -c $^ -o $@

.PHONY: clean
clean:
	-rm -v $(project)
