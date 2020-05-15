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
   integer :: ib, id, ierr, world_rank
   type(WB_State) :: s

   call mpi_init( ierr )
   call check_input_file( filename )
   call initialize_state( s, filename )

   if ( s%world_rank .eq. WORLD_MASTER ) then
      do ib = 1, s%nb
         write (*,"(A)") "----------------------------------------"
         write (*,"(A, I1)") "Block ", ib
         write (*,"(A, I1)") "block_size = ", s%blocks(ib)%block_size
         write (*,"(A)", advance="no") "np         = "
         do id = 1, ND
            write(*,"(I3, A)", advance="no") s%blocks(ib)%np(id), ", "
         end do
         write (*,"(A)") ""
         write (*,"(A)", advance="no") "nx         = "
         do id = 1, ND
            write(*,"(I4, A)", advance="no") s%blocks(ib)%nx(id), ", "
         end do
         write (*,"(A)") ""
         write (*,"(A, L)")  "reorder    = ", s%blocks(ib)%reorder
      end do
      write (*,"(A)") "----------------------------------------"
      do world_rank = 0, s%world_size-1
         write (*,"(A, I3, A, I2)", advance="no") "Process ", world_rank, &
            ": ", s%processes(world_rank)%ib
         do id = 1, ND
            write(*,"(I3, A)", advance="no") &
               s%processes(world_rank)%block_coords(id), ", "
         end do
         do id = 1, ND
            write(*,"(I4, A)", advance="no") &
               s%processes(world_rank)%nx(id), ", "
         end do
         write (*,"(A)") ""
      end do
      write (*,"(A)") "----------------------------------------"
   end if

   call deallocate_state( s )
   call mpi_finalize( ierr )
end program windbag
