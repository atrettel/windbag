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
   use iso_fortran_env
   use mpi_f08
   use wbbase

   implicit none
   character(len=64) :: filename
   integer :: argc, filename_length, ierr, mpi_rank
   logical :: file_exists

   call mpi_init( ierr )
   call mpi_comm_rank( mpi_comm_world, mpi_rank, ierr )
   if ( mpi_rank .eq. MPI_MASTER ) then
      argc = command_argument_count()
      if ( argc .eq. 0 ) then
         write (*,"(A)") "Usage: windbag [INPUT_FILE]"
         call mpi_abort( mpi_comm_world, EXIT_SUCCESS, ierr )
      end if

      call get_command_argument( 1, filename, filename_length, ierr )
      inquire( file=filename, exist=file_exists, iostat=ierr )
      if ( file_exists .eqv. .false. ) then
         write (*,"(A)") "windbag: input file does not exist"
         call mpi_abort( mpi_comm_world, EXIT_FAILURE, ierr )
      end if
   end if

   call mpi_barrier( mpi_comm_world, ierr )
   write (*,"(A, I4, A)") "windbag: input file exists (mpi_rank = ", mpi_rank, ")"
   call mpi_finalize( ierr )
end program windbag
