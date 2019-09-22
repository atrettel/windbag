! Copyright (C) 2019 Andrew Trettel.  All rights reserved.
module wbbase
   use iso_fortran_env
   use mpi

   implicit none

   ! FP = float precision
   integer(4), parameter ::     FP =     real64
   integer(4), parameter :: MPI_FP = mpi_real8

   ! IP = integer precision
   integer(4), parameter ::     IP =     int32
   integer(4), parameter :: MPI_IP = mpi_integer4

   integer(IP), parameter :: STRING_LENGTH = 64

   integer(IP), parameter :: ROOT_PROCESS_RANK   = 0
   integer(IP), parameter :: ROOT_PROCESS_NUMBER = 1

   type WB_Field_Data
      integer(IP), public :: nx_global, ny_global, nz_global
      integer(IP), public :: nx_local,  ny_local,  nz_local
      integer(IP), public :: n_proc, n_proc_x, n_proc_y, n_proc_z
      integer(IP), public :: i_proc, i_proc_x, i_proc_y, i_proc_z
   end type WB_Field_Data

   interface WB_Field_Data
      module procedure init_WB_Field_Data
   end interface WB_Field_Data
contains
   subroutine boot_windbag( input_file_name )
      integer(IP) :: number_of_arguments, current_process_rank, error_status
      character(len=STRING_LENGTH), intent(out) :: input_file_name
      logical :: input_file_exists

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
   end subroutine boot_windbag

   function init_WB_Field_Data( nx_global, ny_global, nz_global ) &
   result( field_data )
      integer(IP), intent(in) :: nx_global, ny_global, nz_global
      integer(IP) :: n_proc, current_process_rank, error_status
      type(WB_Field_Data) :: field_data

      call mpi_comm_size( mpi_comm_world, n_proc, error_status )
      call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )

      field_data%n_proc    = n_proc
      field_data%i_proc    = current_process_rank + 1
      field_data%nx_global = nx_global
      field_data%ny_global = ny_global
      field_data%nz_global = nz_global
   end function init_WB_Field_Data

   subroutine stop_windbag( message )
      integer(IP) :: current_process_rank, error_status
      character(len=*), intent(in) :: message
   
      call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )
   
      if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
         write (*, "(A, A)") "windbag: ", message
      endif
   
      call mpi_finalize( error_status )
      stop
   end subroutine stop_windbag
end module wbbase
