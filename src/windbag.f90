! Copyright (C) 2020 Andrew Trettel
program windbag
   use iso_fortran_env
   use mpi
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
