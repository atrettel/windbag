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

   integer(IP), parameter ::       EXIT_FAILURE = 1_IP
   integer(IP), parameter ::       EXIT_SUCCESS = 0_IP
   integer(IP), parameter ::     STATUS_SUCCESS = 0_IP
   integer(IP), parameter :: MPI_STATUS_SUCCESS = 0_IP

   integer(IP), parameter :: STRING_LENGTH = 64_IP

   integer(IP), parameter :: ROOT_PROCESS_RANK   = 0_IP
   integer(IP), parameter :: ROOT_PROCESS_NUMBER = 1_IP

   integer(IP), parameter :: INPUT_FILE_UNIT = 100_IP

   type WB_Field_Data
      integer(IP), public :: number_of_components, number_of_ghost_points
      integer(IP), public :: nx_global, ny_global, nz_global
      integer(IP), public :: nx_local,  ny_local,  nz_local
      integer(IP), public :: n_proc, n_proc_x, n_proc_y, n_proc_z
      integer(IP), public :: i_proc, i_proc_x, i_proc_y, i_proc_z

      real(FP), public :: time
   end type WB_Field_Data

   interface WB_Field_Data
      module procedure init_WB_Field_Data
   end interface WB_Field_Data
contains
   subroutine boot_program( input_file_name )
      integer(IP) :: number_of_arguments, input_file_name_length
      integer(IP) :: current_process_rank, error_status
      character(len=STRING_LENGTH), intent(out) :: input_file_name
      logical :: input_file_exists

      call mpi_init( error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( "error starting MPI during boot process", &
            EXIT_FAILURE )
      end if

      call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( "error getting process rank during boot process", &
            EXIT_FAILURE )
      end if

      ! Check if there are any command line arguments.
      if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
         number_of_arguments = command_argument_count()
      end if

      call mpi_bcast( number_of_arguments, 1_IP, MPI_IP, ROOT_PROCESS_RANK, &
                      mpi_comm_world, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error broadcasting number of command line arguments", &
            EXIT_FAILURE )
      end if

      if ( number_of_arguments .eq. MPI_STATUS_SUCCESS ) then
         call stop_program( "no command line argument given", EXIT_SUCCESS )
      end if

      ! Check if the first command line argument exists.
      if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
         call get_command_argument( 1_IP, input_file_name, &
            input_file_name_length, error_status )
         if ( error_status .ne. STATUS_SUCCESS ) then
            call stop_program( "error reading first command line argument", & 
               EXIT_FAILURE )
         end if

         inquire( file=input_file_name, exist=input_file_exists, &
            iostat=error_status )
         if ( error_status .ne. STATUS_SUCCESS ) then
            call stop_program( "error inquiring if input file exists", & 
               EXIT_FAILURE )
         end if
      end if

      call mpi_bcast( input_file_name, STRING_LENGTH, mpi_char, &
         ROOT_PROCESS_RANK, mpi_comm_world, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( "error broadcasting input file name", & 
            EXIT_FAILURE )
      end if

      call mpi_bcast( input_file_exists, 1, mpi_logical, ROOT_PROCESS_RANK, &
                      mpi_comm_world, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error broadcasting whether the input file exists", &
            EXIT_FAILURE )
      end if

      if ( input_file_exists .eqv. .false. ) then
         call stop_program( "input file does not exist", EXIT_FAILURE )
      end if
   end subroutine boot_program

   function count_integer_factors( n ) result( n_factors )
      integer(IP), intent(in) :: n
      integer(IP) :: i, n_factors

      n_factors = 0_IP
      do i = 1_IP, n
         if ( mod( n, i ) .eq. 0_IP ) then
            n_factors = n_factors + 1
         end if
      end do
   end function count_integer_factors

   subroutine integer_factorization( n, factors )
      integer(IP), intent(in) :: n
      integer(IP), dimension(:), allocatable, intent(out) :: factors
      integer(IP) :: i, i_factor, n_factors

      n_factors = count_integer_factors( n )

      allocate( factors(n_factors) )

      i_factor = 1
      do i = 1_IP, n
         if ( mod( n, i ) .eq. 0_IP ) then
            factors(i_factor) = i
            i_factor = i_factor + 1
         end if
      end do
   end subroutine integer_factorization

   function init_WB_Field_Data( number_of_components, number_of_ghost_points, &
      nx_global, ny_global, nz_global ) result( field_data )
      integer(IP), intent(in) :: number_of_components, number_of_ghost_points
      integer(IP), intent(in) :: nx_global, ny_global, nz_global
      integer(IP) :: n_proc, current_process_rank, error_status
      type(WB_Field_Data) :: field_data

      call mpi_comm_size( mpi_comm_world, n_proc, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error getting number of processes", &
            EXIT_FAILURE )
      end if

      call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error getting process rank before creating field data object", &
            EXIT_FAILURE )
      end if

      field_data%n_proc                 = n_proc
      field_data%i_proc                 = current_process_rank + 1
      field_data%number_of_components   = number_of_components
      field_data%number_of_ghost_points = number_of_ghost_points
      field_data%nx_global              = nx_global
      field_data%ny_global              = ny_global
      field_data%nz_global              = nz_global

      field_data%time = 0.0_FP
   end function init_WB_Field_Data

   subroutine read_general_namelist( input_file_name, casename, &
      number_of_components, number_of_ghost_points, &
      nx_global, ny_global, nz_global )
      character(len=STRING_LENGTH), intent(in) :: input_file_name
      character(len=STRING_LENGTH), intent(out) :: casename
      integer(IP), intent(out) :: number_of_components, number_of_ghost_points
      integer(IP), intent(out) :: nx_global, ny_global, nz_global
      integer(IP) :: current_process_rank, error_status
      namelist /general/ casename, number_of_components, &
         number_of_ghost_points, nx_global, ny_global, nz_global

      ! Default values, in effect
      casename = ""
      number_of_components = 1_IP
      number_of_ghost_points = 3_IP
      nx_global = 0_IP
      ny_global = 0_IP
      nz_global = 0_IP

      call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error getting process rank before reading general namelist", &
            EXIT_FAILURE )
      end if

      if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
         open(                    &
            unit=INPUT_FILE_UNIT, &
            file=input_file_name, &
            form="formatted",     &
            delim="quote",        &
            action="read",        &
            iostat=error_status   &
         )
         if ( error_status .ne. STATUS_SUCCESS ) then
            call stop_program( &
               "error opening input file to read general namelist", &
               EXIT_FAILURE )
         end if

         read( unit=INPUT_FILE_UNIT, nml=general, iostat=error_status )
         if ( error_status .ne. STATUS_SUCCESS ) then
            call stop_program( "error reading general namelist", &
               EXIT_FAILURE )
         end if

         close( unit=INPUT_FILE_UNIT, iostat=error_status )
         if ( error_status .ne. STATUS_SUCCESS ) then
            call stop_program( &
               "error closing input file after reading general namelist", &
               EXIT_FAILURE )
         end if
      end if

      call mpi_bcast( number_of_components, 1_IP, MPI_IP, ROOT_PROCESS_RANK, &
                      mpi_comm_world, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error broadcasting number_of_components", &
            EXIT_FAILURE )
      end if

      call mpi_bcast( number_of_ghost_points, 1_IP, MPI_IP, &
         ROOT_PROCESS_RANK, mpi_comm_world, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error broadcasting number_of_ghost_points", &
            EXIT_FAILURE )
      end if

      call mpi_bcast( nx_global, 1_IP, MPI_IP, ROOT_PROCESS_RANK, &
                      mpi_comm_world, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error broadcasting nx_global", &
            EXIT_FAILURE )
      end if

      call mpi_bcast( ny_global, 1_IP, MPI_IP, ROOT_PROCESS_RANK, &
                      mpi_comm_world, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error broadcasting ny_global", &
            EXIT_FAILURE )
      end if

      call mpi_bcast( nz_global, 1_IP, MPI_IP, ROOT_PROCESS_RANK, &
                      mpi_comm_world, error_status )
      if ( error_status .ne. MPI_STATUS_SUCCESS ) then
         call stop_program( &
            "error broadcasting nz_global", &
            EXIT_FAILURE )
      end if

      if ( number_of_components .lt. 1 ) then
         call stop_program( "number_of_components < 1", &
            EXIT_FAILURE )
      end if 

      if ( number_of_ghost_points .lt. 1 ) then
         call stop_program( "number_of_ghost_points < 1", &
            EXIT_FAILURE )
      end if 

      if ( nx_global * ny_global * nz_global .eq. 0_IP ) then
         call stop_program( "grid has no points in at least one direction", &
            EXIT_FAILURE )
      end if
   end subroutine read_general_namelist

   subroutine stop_program( message, exit_status )
      integer(IP) :: current_process_rank, error_status
      character(len=*), intent(in) :: message
      integer(IP), intent(in) :: exit_status
   
      call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )
   
      if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
         if ( exit_status .eq. EXIT_SUCCESS ) then
            write (output_unit, "(A, A)") "windbag: ", message
         else
            write (error_unit, "(A, A)") "windbag: ", message
         end if
      end if
   
      call mpi_finalize( error_status )

      if ( exit_status .eq. EXIT_SUCCESS ) then
         stop
      else
         error stop
      end if
   end subroutine stop_program
end module wbbase
