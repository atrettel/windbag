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
   integer :: ierr, i_rank
   type(WB_State) :: s

   call check_input_file( filename )
   call initialize_state( s, filename )

   do i_rank = 0, s%world_size-1
      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( s%world_rank .eq. i_rank ) then
         write (*,"(A, I4)") "world_rank = ", s%world_rank
         write (*,"(A, A)") "case_name = ", s%case_name
         write (*,"(A, I4)") "nb = ", s%nb
         write (*,"(A, I4)") "ng = ", s%ng
      end if
   end do

   call mpi_finalize( ierr )
end program windbag
