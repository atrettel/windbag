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
   contains
      procedure, private :: decompose_domain => WB_Field_Data_decompose_domain
   end type WB_Field_Data

   interface WB_Field_Data
      module procedure init_WB_Field_Data
   end interface WB_Field_Data
contains
   function init_WB_Field_Data( n_proc, nx, ny, nz ) result( field_data )
      integer(IP), intent(in) :: n_proc, nx, ny, nz
      type(WB_Field_Data) :: field_data

      call field_data%decompose_domain( n_proc, nx, ny, nz )
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

   subroutine WB_Field_Data_decompose_domain( self, n_proc, nx, ny, nz )
      class(WB_Field_Data), intent(inout) :: self
      integer(IP), intent(in) :: n_proc, nx, ny, nz

      self%n_proc = n_proc
      self%nx     = nx
      self%ny     = ny
      self%nz     = nz
   end subroutine WB_Field_Data_decompose_domain
end module wbbase
