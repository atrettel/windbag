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
module wb_text
   implicit none

   integer, public, parameter ::      UNALIGNED = 0
   integer, public, parameter ::   LEFT_ALIGNED = 1
   integer, public, parameter ::  RIGHT_ALIGNED = 2
   integer, public, parameter :: CENTER_ALIGNED = 3
contains
   subroutine write_string_table_entry( f, entry, width, end_row )
      integer, intent(in) :: f, width
      character(len=*), intent(in) :: entry
      logical, optional, intent(in) :: end_row
      character(len=64) :: write_fmt

      write (write_fmt,"(A, I2.1, A)") "(A, A", width, ", A)"
      write (f, write_fmt, advance="no") "| ", entry, " "
      if ( present(end_row) .and. end_row .eqv. .true. ) then
         write (f, "(A)", advance="yes") "|"
      end if
   end subroutine write_string_table_entry

   subroutine write_integer_table_entry( f, entry, width, end_row )
      integer, intent(in) :: f, entry, width
      logical, optional, intent(in) :: end_row
      character(len=64) :: write_fmt

      write (write_fmt,"(A, I2.1, A)") "(A, I", width, ", A)"
      write (f, write_fmt, advance="no") "| ", entry, " "
      if ( present(end_row) .and. end_row .eqv. .true. ) then
         write (f, "(A)", advance="yes") "|"
      end if
   end subroutine write_integer_table_entry

   subroutine write_rule_table_entry( f, width, alignment, end_row )
      integer, intent(in) :: f, width
      integer, optional, intent(in) :: alignment
      logical, optional, intent(in) :: end_row
      logical, dimension(2) :: aligned = (/ .false., .false. /)
      integer :: n_dash, i_dash

      n_dash = width
      if ( present(alignment) ) then
         if      ( alignment .eq.   LEFT_ALIGNED ) then
            aligned = (/ .true., .false. /)
         else if ( alignment .eq.  RIGHT_ALIGNED ) then
            aligned = (/ .false., .true. /)
         else if ( alignment .eq. CENTER_ALIGNED ) then
            aligned = (/ .true., .true. /)
         end if
      end if
      n_dash = n_dash - count(aligned)

      write (f, "(A)", advance="no") "| "
      if ( aligned(1) ) then
         write (f, "(A)", advance="no") ":"
      end if
      do i_dash = 1, n_dash
         write (f, "(A)", advance="no") "-"
      end do
      if ( aligned(2) ) then
         write (f, "(A)", advance="no") ":"
      end if
      write (f, "(A)", advance="no") " "
      if ( present(end_row) .and. end_row .eqv. .true. ) then
         write (f, "(A)", advance="yes") "|"
      end if
   end subroutine write_rule_table_entry
end module wb_text
