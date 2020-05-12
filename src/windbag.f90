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
   integer :: argc, filename_length, ierr, world_rank
   logical :: file_exists
   integer :: mpi_float_size

   call mpi_init( ierr )
   call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr )
   if ( world_rank .eq. WORLD_MASTER ) then
      argc = command_argument_count()
      if ( argc .eq. 0 ) then
         write (*,"(A)") "Usage: windbag [INPUT_FILE]"
         call mpi_abort( MPI_COMM_WORLD, MPI_SUCCESS, ierr )
      end if

      call get_command_argument( 1, filename, filename_length, ierr )
      inquire( file=filename, exist=file_exists, iostat=ierr )
      if ( file_exists .eqv. .false. ) then
         write (*,"(A)") "windbag: input file does not exist"
         call mpi_abort( MPI_COMM_WORLD, MPI_ERR_NO_SUCH_FILE, ierr )
      end if
   end if

   call mpi_barrier( MPI_COMM_WORLD, ierr )
   write (*,"(A, I4, A)") "windbag: input file exists (world_rank = ", &
      world_rank, ")"

   call mpi_sizeof( 1.0_FP, mpi_float_size, ierr )
   call mpi_type_match_size( MPI_TYPECLASS_REAL, mpi_float_size, MPI_FP, ierr )

   if ( world_rank .eq. WORLD_MASTER ) then
      print *, MPI_REAL
      print *, MPI_REAL4
      print *, MPI_REAL8
      print *, MPI_FP
   end if

   call mpi_finalize( ierr )
end program windbag
