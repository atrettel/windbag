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

   integer(IP), parameter :: ROOT_PROCESS_RANK = 0

   type WB_Field_Data
      integer(IP), public :: n_proc, nx, ny, nz
   end type WB_Field_Data
contains
   subroutine stop_windbag( message )
      integer(IP) :: current_process_rank, error_status
      character(len=*), intent(in) :: message
   
      call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )
   
      if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
         write (*, "(A, A)") "windbag: ", message
      endif
   
      call mpi_finalize( error_status )
      stop
   end subroutine
end module wbbase
