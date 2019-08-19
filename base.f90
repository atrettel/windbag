! Copyright (C) 2019 Andrew Trettel.  All rights reserved.

module base
   use iso_fortran_env
   use mpi

   ! FP = float precision
   integer(4), parameter :: FP = real64

   ! IP = integer precision
   integer(4), parameter :: IP = int32

   integer(IP), parameter :: ROOT_PROCESS_RANK = 0
contains
   subroutine stop_windbag( message )
   
      implicit none
      integer(IP) :: current_process_rank, error_status
      character(len=*), intent(in) :: message
   
      call mpi_comm_rank( mpi_comm_world, current_process_rank, error_status )
   
      if ( current_process_rank .eq. ROOT_PROCESS_RANK ) then
         write (*, "(A, A)") "windbag: ", message
      endif
   
      call mpi_finalize( error_status )
      stop
   end subroutine
end module base
