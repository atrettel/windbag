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
module wb_exit
   use iso_fortran_env, only : error_unit
   use mpi_f08, only : MPI_COMM_WORLD, mpi_abort
   use wb_representation
   implicit none

   private

   public wb_abort

   integer(MP), public, parameter :: EXIT_SUCCESS =  0_MP
   integer(MP), public, parameter :: EXIT_FAILURE =  1_MP
   integer(MP), public, parameter :: EXIT_USAGE   = 64_MP
   integer(MP), public, parameter :: EXIT_DATAERR = 65_MP
   integer(MP), public, parameter :: EXIT_NOINPUT = 66_MP
contains
   subroutine wb_abort( message, exit_code, ints, floats )
      character(len=*), intent(in) :: message
      integer(MP), intent(in) :: exit_code
      integer(SP) :: i
      integer(MP) :: ierr
      integer(SP), dimension(:), optional, intent(in) :: ints
      real(FP), dimension(:), optional, intent(in) :: floats

      write (error_unit, "(A, A, A)") "fatal error", ": ", message
      if ( present(ints) ) then
         do i = 1_SP, size(ints)
            write (error_unit, "(A, I1, A, I8)") "N", i, " = ", ints(i)
         end do
      end if
      if ( present(floats) ) then
         do i = 1_SP, size(floats)
            write (error_unit, "(A, I1, A, ES9.2)") "F", i, " = ", floats(i)
         end do
      end if
      call mpi_abort( MPI_COMM_WORLD, exit_code, ierr )
   end subroutine wb_abort
end module wb_exit
