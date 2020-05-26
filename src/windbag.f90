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
   use wb_representation, only : STRING_LENGTH
   use wb_base

   implicit none
   character(len=STRING_LENGTH) :: filename
   integer :: ierr
   type(WB_State) :: s

   call mpi_init( ierr )
   call check_input_file( filename )
   call initialize_state( s, filename )
   call print_initial_information( s )
   call deallocate_state( s )
   call mpi_finalize( ierr )
end program windbag
