! Copyright (C) 2019 Andrew Trettel.  All rights reserved.
program windbag
   use wbbase

   implicit none
   integer(IP) :: number_of_arguments, current_process_rank, error_status
   character(len=64) :: input_file_name
   logical :: input_file_exists
   type(WB_Field_Data) :: field_data

   call mpi_init( error_status )

   call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )

   ! Check if there are any command line arguments.
   if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
      number_of_arguments = command_argument_count()
   end if

   call mpi_bcast( number_of_arguments, 1, MPI_IP, ROOT_PROCESS_RANK, &
                   mpi_comm_world, error_status )

   if ( number_of_arguments .eq. 0 ) then
      call stop_windbag( "no argument given" )
   end if

   ! Check if the first command line argument exists.
   if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
      call get_command_argument( 1, input_file_name )
      inquire( file=input_file_name, exist=input_file_exists )
   end if

   call mpi_bcast( input_file_name, 64, mpi_char, ROOT_PROCESS_RANK, &
                   mpi_comm_world, error_status )
   call mpi_bcast( input_file_exists, 1, mpi_logical, ROOT_PROCESS_RANK, &
                   mpi_comm_world, error_status )

   if ( input_file_exists .eqv. .false. ) then
      call stop_windbag( "input file does not exist" )
   end if

   field_data = WB_Field_Data( 128, 128, 128 )

   if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
      print *, field_data%n_proc
      print *, field_data%i_proc
      print *, field_data%nx_global
      print *, field_data%ny_global
      print *, field_data%nz_global
   end if

   call stop_windbag( "input file exists" )
end program windbag
