! Copyright (C) 2019 Andrew Trettel.  All rights reserved.

subroutine stop_windbag( message )
   use mpi

   implicit none
   integer(4) :: process_rank, root_process_rank=0, error_status
   character(len=*), intent(in) :: message

   call mpi_comm_rank( mpi_comm_world, process_rank, error_status )

   if ( process_rank .eq. root_process_rank ) then
      print *, "windbag: ", message
   endif

   call mpi_finalize( error_status )
   stop
end subroutine

program windbag
   use mpi

   implicit none
   integer(4) :: number_of_arguments, process_rank, root_process_rank=0, error_status
   character(len=64) :: input_file_name
   logical :: input_file_exists

   call mpi_init( error_status )

   call mpi_comm_rank( mpi_comm_world, process_rank, error_status )

   ! Check if there are any command line arguments.
   if ( process_rank .eq. root_process_rank ) then
      number_of_arguments = command_argument_count()
   end if

   call mpi_bcast( number_of_arguments, 1, mpi_integer, root_process_rank, &
                   mpi_comm_world, error_status )

   if ( number_of_arguments .eq. 0 ) then
      call stop_windbag( "no argument given" )
   end if

   ! Check if the first command line argument exists.
   if ( process_rank .eq. root_process_rank ) then
      call get_command_argument( 1, input_file_name )
      inquire( file=input_file_name, exist=input_file_exists )
   end if

   call mpi_bcast( input_file_name, 64, mpi_char, root_process_rank, &
                   mpi_comm_world, error_status )
   call mpi_bcast( input_file_exists, 1, mpi_logical, root_process_rank, &
                   mpi_comm_world, error_status )

   if ( input_file_exists .eqv. .false. ) then
      call stop_windbag( "input file does not exist" )
   end if

   call stop_windbag( "input file exists" )
end program windbag
