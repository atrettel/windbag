Windbag design
==============


Objective
----------

Create a simple, clear, modular, and economical direct numerical simulation
code built from the ground up in modern Fortran.


Background
----------

I have used several direct numerical simulation codes in the past, but I never
have written one from scratch.  Other DNS codes exist, but I have found they
tend to concentrate on a different set of goals than the ones that matter to me
personally (see the list of goals below).  To that end, I realized that I need
to create my own DNS code to satisfy these goals.


Goals
-----

- Arbitrary-order statistics from post-processing

    - The code MUST output all statistics for all fields up to a particular
      order (likely 5th-order statistics).  The code SHOULD output both
      unweighted and density-weighted statistics, though the code MAY only
      output unweighted statistics provided the order is high enough to recover
      the density-weighted statistics at a lower order.

- Block structured grid

    - The code MUST use a block structured grid.

    - Originally I considered the much simpler design of using a single
      structured grid.  A single structure grid works for most scenarios,
      especially for periodic boxes, channels, and flat plates, but for cases
      like backward-facing steps or external flows it may prove impractical.

    - In light of this problem, I picked a simple way to handle block
      structured grids.  Each process belongs to a single block, and block
      interfaces are treated simply as process interfaces.  This allows for
      different blocks to be treated like different processes (provided
      additional constraints are in place).  This tradeoff allows for both
      additional flexibility while keeping the overall design simple.  The
      tradeoff is that it proves impractical for real-world geometries (complex
      geometries), since it is difficult to satisfy all of the needed
      constraints with a large number of blocks.

- Modularity

    - The code SHOULD set variables and procedures to private whenever
      possible.

    - The code SHOULD separate use and implementation in data structures (data
      abstraction).

    - The code MUST allow for arbitrary equations of state and constitutive
      relations.

    - The code MUST use data structures (and not just arrays) to organize data,
      and the code MUST avoid global variables.  Most of the previous Fortran
      DNS codes I have seen use global variables extensively.  Side effects
      could become a major issue.

- Numerical methods

    - The code MUST use the compressible Navier-Stokes equations for the
      governing equations.  The code MUST use conservative variables for the
      state variables.

    - The code SHOULD use a high-order finite difference method (primarily).
      It SHOULD be possible to set the order arbitrarily, and it SHOULD be
      possible to change numerical methods if desired (modularity).  A command
      design pattern satisfies this goal.

- Openness, accessibility, and reproducibility
      
    - Main goal: reproducibility by a third party.  Peer reviewers MUST be able
      to examine and use the source code before publication and readers MUST be
      able to examine and use the source code after publication (if they would
      like to).

    - The code SHOULD avoid external dependencies (especially proprietary
      dependencies) that prevent other people from being able to compile and
      review the code and reproduce its results.  The code SHOULD use no
      external dependencies other than MPI.

    - The code MUST use open formats like binary and XDMF XML files for the
      outputted field data.  The code MUST ensure that the output is readable
      by open source visualization programs like Paraview.  The code MUST use
      open formats like CSV for post-processed output.  These steps ensure that
      any results remain accessible in the long-term.

    - The code MUST a single, human-readable, input file for each case.  The
      input file SHOULD use Fortran namelists (for convenience).

- Software architecture

    - The code MUST Fortran rather than C or C++.

        - Fortran is antiquated, but it offers several advantages for numerical
          work, including namelists, simpler memory management, ability to
          change precision easily, and arrays as first class citizens.  The
          tradeoff is that C and C++ handle complex data structures much more
          easily.  This project involves arrays primarily, so Fortran seems
          more appropriate here.

    - The code MUST use a simple Makefile for compilation.  Compilation MUST
      NOT require the user to specify any environment variables.

    - The code MUST NOT use any preprocessor.

    - The same executable SHOULD be used for all cases.  There SHOULD be no
      need to recompile the code at all for different "cases" (types of
      simulations).  The same executable SHOULD be used for pre-processing,
      time advancement, and post-processing (excluding visualization).

    - The code MUST accept a single command line argument, the input file.  The
      filename for the input file MUST NOT be hard coded.  The code MUST use
      case names (labels) to organize the data rather than separate
      directories.  These practices allow for better data organization,
      especially when different cases merely change a single parameter like the
      Mach or Reynolds number or some aspect of the grid.  These practices also
      allow for the same executable to be used for all cases.

    - The code MUST only use MPI for parallelization.

    - The code SHOULD use descriptive variable names (even if they are long).


Non-goals
---------

These are beyond the scope of the project, so I will not address them in the
design, nor will the project attempt to satisfy them.  This is somewhat a list
of tradeoffs in the proposed design.

- Additional file formats

    - Binary and CSV are both archival and sufficient.  Other file formats
      occasionally require additional libraries, and this also violates the
      goal of having no external dependencies other than MPI.

    - That said, I am not against eventually using other open formats like HDF5
      with the code, but I find using binary sufficient for the time being,
      Binary has all of the properties the code needs without having to use any
      additional libraries.

- Real-world geometries (complex geometries)

    - The project is used for research purposes only and not for design
      purposes.  As such, I can accept severe limitations in the kind of
      geometries that can be represented provided that they assist in other
      aspects (especially in terms of using different numerical methods).

- Speed

    - While speed certainly is always a concern, I would accept a slower speed
      as a tradeoff for easier-to-maintain code.

    - A consequence of this is that the project will not concentrate on using
      the latest and greatest technologies (GPUs in particular).  This tradeoff
      keeps the code simple by letting MPI handle parallelism.

- Test bed for numerical methods development

    - The code should only use established numerical methods published in
      journals or commonly taught at universities.  I do not intend to use the
      code to test out new numerical methods or to develop new numerical
      methods in particular.

    - However, the command design pattern allows for different numerical
      methods to be easily swapped in or out depending on the circumstances.


Design details
--------------

The high-level architecture of Windbag generally follows the *command design
pattern*.  The program performs a series of operations on a subdomain data
structure.  The subdomain data structure contains all of the field data and
state information, but it contains no information about the numerical methods
(operations that change the data) or post-processing operations (operations
that use the data).  Each command then operates in sequence and either modifies
the subdomain data structure or performs some calculation based on it.

The sequence of operations:

1. Initialize subdomain data structure.

2. Run all pre-processing commands.

    - Pre-processing commands include grid generation and setting the initial
      conditions.

3. Run post-processing commands on the initial subdomain.

4. Run the numerical methods command.

    - Each input file should specify only one numerical method to use.

5. Update additional fields based on dependency graph.

    - The variable list data structure uses a dependency graph, currently
      implemented as an adjacency matrix, to determine the order of evaluation
      needed to update any additional instantaneous fields.  This allows the
      user to select which fields they want, and then the data structure
      figures out what the dependencies are (how to calculate them).  The end
      result is that the fields are sequentially put in the order they need to
      be evaluated in, so the field numbers are meaningful and not arbitrary.

    - Circular dependencies should be impossible given that the dependency
      graph is hard-coded and I ensured that there is only one way for each
      field to be calculated.

6. Save fields to file system (if needed).

7. Run post-processing commands.

    - Post-processing includes averaging and taking probes at various points.
      There could be multiple post-processing commands.

8. Repeat steps 4 to 7 until completion.

    - Note that this entire iterative process could be represented as a queue,
      provided that after each iteration from steps 4 to 7 it enqueues new
      commands.


Additional notes and discussion
-------------------------------

- The command design pattern would be easier to implement in a functional
  programming language.  However, in Fortran this may be possible if the code
  uses a standardized interface and passes subroutines as arguments.  The
  precise implementation is still unclear to me.

- The command design pattern leaves a log and makes it clear what the program
  is doing at one moment or another.  This aids in repeatability.

- I have considered using the subdomain data structure as the output of the
  averaging operation, but this presents several issues.  The primary advantage
  is that it would re-use existing code (economy or minimalism).  The
  disadvantage is that the variable lists do not match in each case.  The
  instantaneous case has a limited number of fields, all instantaneous in time,
  while the averaged case has a much larger number of fields.  There also is
  the issue of the number of dimensions of the averaged field.  If a simulation
  should be statistically homogeneous in one direction, then its averaged field
  has only a single point in that direction.  This motivates having a separate
  domain decomposition for that (otherwise there could be multiple processes in
  the homogeneous direction, all sharing the same point).

    - Instead, the more obvious route would be to have each block leader manage
      the post-processing for each block.  Reductions from each worker would
      end up at their local block leader.  For repeatability, each timestep
      should output the spatially-averaged fields.  To time-average, the code
      should then average over all of these fields.

- Most of the DNS cases I want to run are 3D, but I also wanted to allow the
  code to be 2D if possible.  To that end, I have made the field array of the
  subdomain data structure to always have 3 dimensions, but the number of
  dimensions is still variable (the third dimension would only be 1 point
  wide).  This works for now but I need to look into this more.


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
