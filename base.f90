! Copyright (C) 2019 Andrew Trettel.  All rights reserved.

module base
   use iso_fortran_env
   implicit none

   ! FP = float precision
   !integer, parameter :: FP = kind(1.0d0)
   integer, parameter :: FP = real64

   ! IP = integer precision
   integer, parameter :: IP = int32

   integer, parameter :: ROOT_PROCESS_RANK = 0
contains
   subroutine stop_windbag( message )
      use mpi
   
      implicit none
      integer(IP) :: process_rank, error_status
      character(len=*), intent(in) :: message
   
      call mpi_comm_rank( mpi_comm_world, process_rank, error_status )
   
      if ( process_rank .eq. ROOT_PROCESS_RANK ) then
         write (*, "(A, A)") "windbag: ", message
      endif
   
      call mpi_finalize( error_status )
      stop
   end subroutine
end module base
