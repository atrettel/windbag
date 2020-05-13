! Copyright (C) 2020 Andrew Trettel
! 
! Windbag is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
! 
! Windbag is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with
! Windbag.  If not, see <https://www.gnu.org/licenses/>.
program windbag
   use mpi_f08
   use wbbase

   implicit none
   character(len=STRING_LENGTH) :: filename
   integer :: ierr

   call check_input_file( filename )

   write (*,"(A)") "windbag: input file exists"

   call find_mpi_fp

   print *, "MPI_REAL  = ", MPI_REAL
   print *, "MPI_REAL4 = ", MPI_REAL4
   print *, "MPI_REAL8 = ", MPI_REAL8
   print *, "MPI_FP    = ", MPI_FP

   call mpi_finalize( ierr )
end program windbag
