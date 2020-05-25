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
module wb_representation
   use iso_fortran_env, only : int32, int64, real32, real64
   use mpi_f08, only : mpi_sizeof, mpi_type_match_size, MPI_Datatype, &
      MPI_TYPECLASS_INTEGER, MPI_TYPECLASS_REAL
   implicit none

   private

   public find_mpi_precisions

   integer, public, parameter       ::     FP = real64
   type(MPI_Datatype), public, save :: MPI_FP

   integer, public, parameter       ::     SP = int64
   type(MPI_Datatype), public, save :: MPI_SP

   integer, public, parameter       ::     MP = int32
   type(MPI_Datatype), public, save :: MPI_MP

   ! Big-endian architectures put the most significant byte first, and
   ! little-endian architectures put the least significant byte first.  Binary
   ! output requires knowing the architecture's endianness.  This statement
   ! casts an integer for 1 as a single character (the X is arbitrary) and then
   ! gets the ASCII code for the first character (the result only contains the
   ! leading bits of the integer).  If the result is 1, then the architecture
   ! is little-endian, since the least significant byte came first.  If the
   ! result is 0, then the architecture is big-endian.  Little-endian
   ! architectures are more common nowadays.
   logical, public, parameter :: ARCH_IS_BIG_ENDIAN = &
      ichar( transfer( 1_SP, "X" ) ) .eq. 0
contains
   subroutine find_mpi_precisions
      integer(MP) :: mpi_size, ierr

      call mpi_sizeof( 1.0_FP, mpi_size, ierr )
      call mpi_type_match_size( MPI_TYPECLASS_REAL, mpi_size, MPI_FP, ierr )

      call mpi_sizeof( 1_SP, mpi_size, ierr )
      call mpi_type_match_size( MPI_TYPECLASS_INTEGER, mpi_size, MPI_SP, ierr )

      call mpi_sizeof( 1_MP, mpi_size, ierr )
      call mpi_type_match_size( MPI_TYPECLASS_INTEGER, mpi_size, MPI_MP, ierr )
   end subroutine find_mpi_precisions
end module wb_representation
