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
   use mpi_f08
   implicit none

   private

   public find_mpi_precisions,        &
          identify_integer_precision, &
          identify_real_precision,    &
          real_precision_in_bytes

   integer, public, parameter       ::     FP = real64
   type(MPI_Datatype), public, save :: MPI_FP

   integer, public, parameter       ::     SP = int64
   type(MPI_Datatype), public, save :: MPI_SP

   integer, public, parameter       ::     MP = int32
   type(MPI_Datatype), public, save :: MPI_MP

   integer(SP), public, parameter :: STRING_LENGTH = 64_SP

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

   interface identify_integer_precision
      module procedure identify_integer_precision_fortran, &
         identify_integer_precision_mpi
   end interface identify_integer_precision

   interface identify_real_precision
      module procedure identify_real_precision_fortran, &
         identify_real_precision_mpi
   end interface identify_real_precision
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

   function identify_integer_precision_fortran( vp ) &
   result( string )
      integer, intent(in) :: vp
      character(len=STRING_LENGTH) :: string

      if ( vp .eq. int64 ) then
         write (string,"(A)") "`int64`"
      else if ( vp .eq. int32 ) then
         write (string,"(A)") "`int32`"
      else
         write (string,"(A)") "(unknown)"
      end if
   end function identify_integer_precision_fortran

   function identify_integer_precision_mpi( datatype ) &
   result( string )
      type(MPI_Datatype), intent(in) :: datatype
      character(len=STRING_LENGTH) :: string

      if ( datatype .eq. MPI_INTEGER8 ) then
         write (string,"(A)") "`MPI_INTEGER8`"
      else if ( MPI_MP .eq. MPI_INTEGER4 ) then
         write (string,"(A)") "`MPI_INTEGER4`"
      else if ( MPI_MP .eq. MPI_INTEGER ) then
         write (string,"(A)") "`MPI_INTEGER`"
      else
         write (string,"(A)") "(unknown)"
      end if
   end function identify_integer_precision_mpi

   function identify_real_precision_fortran( vp ) &
   result( string )
      integer, intent(in) :: vp
      character(len=STRING_LENGTH) :: string

      if ( vp .eq. real64 ) then
         write (string,"(A)") "`real64`"
      else if ( vp .eq. real32 ) then
         write (string,"(A)") "`real32`"
      else
         write (string,"(A)") "(unknown)"
      end if
   end function identify_real_precision_fortran

   function identify_real_precision_mpi( datatype ) &
   result( string )
      type(MPI_Datatype), intent(in) :: datatype
      character(len=STRING_LENGTH) :: string

      if ( datatype .eq. MPI_REAL8 ) then
         write (string,"(A)") "`MPI_REAL8`"
      else if ( MPI_MP .eq. MPI_REAL4 ) then
         write (string,"(A)") "`MPI_REAL4`"
      else if ( MPI_MP .eq. MPI_REAL ) then
         write (string,"(A)") "`MPI_REAL`"
      else
         write (string,"(A)") "(unknown)"
      end if
   end function identify_real_precision_mpi

   function real_precision_in_bytes( vp ) &
   result( bytes )
      integer, intent(in) :: vp
      integer(SP) :: bytes

      if ( vp .eq. real64 ) then
         bytes = 8_SP
      else if ( vp .eq. real32 ) then
         bytes = 4_SP
      else
         bytes = 0_SP
      end if
   end function
end module wb_representation
