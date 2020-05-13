! Copyright (C) 2020 Andrew Trettel
! 
! Windbag is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
! 
! Windbag is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with
! Windbag.  If not, see <https://www.gnu.org/licenses/>.
module wbbase
   use iso_fortran_env, only : real64
   use mpi_f08
   implicit none

   private

   public WB_State, check_input_file, initialize_state

   integer, public, parameter ::            FP = real64
   integer, public, parameter ::            ND = 3
   integer, public, parameter ::  WORLD_MASTER = 0
   integer, public, parameter :: STRING_LENGTH = 64

   type(MPI_Datatype), public, save :: MPI_FP

   type WB_State
      character(len=STRING_LENGTH) :: case_name
      integer :: block_rank, block_size
      integer :: world_rank, world_size
      integer :: ib, nb
      integer :: ng
   end type WB_State
contains
   subroutine check_input_file( filename )
      character(len=STRING_LENGTH), intent(out)  :: filename
      integer :: argc, filename_length, ierr, world_rank
      logical :: file_exists
      call mpi_init( ierr )
      call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr )
      call find_mpi_fp
      if ( world_rank .eq. WORLD_MASTER ) then
         argc = command_argument_count()
         if ( argc .eq. 0 ) then
            write (*,"(A)") "Usage: windbag [INPUT_FILE]"
            call mpi_abort( MPI_COMM_WORLD, MPI_SUCCESS, ierr )
         end if
         call get_command_argument( 1, filename, filename_length, ierr )
         inquire( file=filename, exist=file_exists, iostat=ierr )
         if ( file_exists .eqv. .false. ) then
            write (*,"(A)") "windbag: input file does not exist"
            call mpi_abort( MPI_COMM_WORLD, MPI_ERR_NO_SUCH_FILE, ierr )
         end if
      end if
      call mpi_barrier( MPI_COMM_WORLD, ierr )
   end subroutine check_input_file

   subroutine find_mpi_fp
      integer :: mpi_float_size, ierr
      call mpi_sizeof( 1.0_FP, mpi_float_size, ierr )
      call mpi_type_match_size( MPI_TYPECLASS_REAL, mpi_float_size, MPI_FP, &
         ierr )
   end subroutine find_mpi_fp

   subroutine initialize_state( s, filename )
      character(len=STRING_LENGTH) :: filename
      integer :: ierr
      type(WB_State), intent(out) :: s

      call mpi_comm_rank( MPI_COMM_WORLD, s%world_rank, ierr )
      call mpi_comm_size( MPI_COMM_WORLD, s%world_size, ierr )
      call read_general_namelist( s, filename )
   end subroutine initialize_state

   subroutine read_general_namelist( s, filename )
      character(len=STRING_LENGTH), intent(in) :: filename
      character(len=STRING_LENGTH) :: case_name
      integer :: ierr, file_unit, nb, ng
      type(WB_State), intent(inout) :: s
      namelist /general/ case_name, nb, ng

      if ( s%world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         read( unit=file_unit, nml=general )

         s%case_name = trim(case_name)
         s%nb = nb
         s%ng = ng
      end if

      call mpi_bcast( s%case_name, len(s%case_name), MPI_CHARACTER, &
         WORLD_MASTER, MPI_COMM_WORLD, ierr )

      call mpi_bcast( s%nb, 1, MPI_INTEGER, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )

      call mpi_bcast( s%ng, 1, MPI_INTEGER, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
   end subroutine read_general_namelist
end module wbbase