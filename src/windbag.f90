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
   use wb_representation
   use wb_base

   implicit none
   character(len=STRING_LENGTH) :: filename
   integer :: ierr
   type(WB_Subdomain) :: sd

   call mpi_init( ierr )
   call find_mpi_precisions
   call find_input_file( filename )
   call wb_subdomain_construct( sd, filename )
   call print_initial_information( sd )
   call wb_subdomain_destroy( sd )
   call mpi_finalize( ierr )
end program windbag
