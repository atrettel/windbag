! Copyright (C) 2020 Andrew Trettel.  All rights reserved.
program windbag
   use iso_fortran_env
   use mpi

   implicit none
   character(len=64) :: filename
   integer :: argc, filename_length, ierr, rank
   logical :: file_exists

   call mpi_init( ierr )
   call mpi_comm_rank( mpi_comm_world, rank, ierr )
   if ( rank .eq. 0 ) then
      argc = command_argument_count()
      if ( argc .eq. 0 ) then
         write (*,"(A)") "Usage: windbag [INPUT_FILE]"
         call mpi_abort( mpi_comm_world, 0, ierr )
      end if

      call get_command_argument( 1, filename, filename_length, ierr )
      inquire( file=filename, exist=file_exists, iostat=ierr )
      if ( file_exists .eqv. .false. ) then
         write (*,"(A)") "windbag: input file does not exist"
         call mpi_abort( mpi_comm_world, 1, ierr )
      end if
   end if

   call mpi_barrier( mpi_comm_world, ierr )
   write (*,"(A, I4, A)") "windbag: input file exists (rank = ", rank, ")"
   call mpi_finalize( ierr )
end program windbag
