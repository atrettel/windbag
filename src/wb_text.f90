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
   use wb_representation
   implicit none

   private

   public write_blank_line, write_log_heading, write_table_entry, &
      write_table_rule_entry

   integer(SP), public, parameter ::      UNALIGNED = 0_SP
   integer(SP), public, parameter ::   LEFT_ALIGNED = 1_SP
   integer(SP), public, parameter ::  RIGHT_ALIGNED = 2_SP
   integer(SP), public, parameter :: CENTER_ALIGNED = 3_SP

   integer(SP), public, parameter ::       ANSWER_COLUMN_WIDTH = 10_SP
   integer(SP), public, parameter :: BLOCK_NUMBER_COLUMN_WIDTH =  7_SP
   integer(SP), public, parameter ::       COORDS_COLUMN_WIDTH =  6_SP
   integer(SP), public, parameter ::    DATA_TYPE_COLUMN_WIDTH = 20_SP
   integer(SP), public, parameter ::     HOSTNAME_COLUMN_WIDTH = 16_SP
   integer(SP), public, parameter ::         NAME_COLUMN_WIDTH = 40_SP
   integer(SP), public, parameter ::           NP_COLUMN_WIDTH =  7_SP
   integer(SP), public, parameter ::           NX_COLUMN_WIDTH =  7_SP
   integer(SP), public, parameter ::       POINTS_COLUMN_WIDTH = 16_SP
   integer(SP), public, parameter ::     PROPERTY_COLUMN_WIDTH = 30_SP
   integer(SP), public, parameter ::     QUESTION_COLUMN_WIDTH = 30_SP
   integer(SP), public, parameter ::         RANK_COLUMN_WIDTH = 12_SP
   integer(SP), public, parameter ::         SIZE_COLUMN_WIDTH = 12_SP
   integer(SP), public, parameter ::        VALUE_COLUMN_WIDTH = 20_SP
   integer(SP), public, parameter ::     VARIABLE_COLUMN_WIDTH = 15_SP

   interface write_table_entry
      module procedure write_table_entry_character, &
         write_table_entry_integer, write_table_entry_logical
   end interface write_table_entry
contains
   subroutine write_log_heading( f, title, level )
      integer, intent(in) :: f
      character(len=*), intent(in) :: title
      integer(SP), optional, intent(in) :: level
      integer(SP) :: i_hash, n_hash

      if ( present(level) ) then
         n_hash = level
      else
         n_hash = 1_SP
      end if

      do i_hash = 1_SP, n_hash
         write (f, "(A)", advance="no") "#"
      end do
      write (f, "(A, A)") " ", trim(title)
      write (f, "(A)" ) ""
   end subroutine write_log_heading

   subroutine write_blank_line( f, n )
      integer, intent(in) :: f
      integer(SP), optional, intent(in) :: n
      integer(SP) :: i_line, total_lines

      if ( present(n) ) then
         total_lines = n
      else
         total_lines = 1_SP
      end if

      do i_line = 1_SP, total_lines
         write (f, "(A)" ) ""
      end do
   end subroutine write_blank_line

   subroutine write_table_entry_character( f, entry, width, end_row )
      integer, intent(in) :: f
      integer(SP), intent(in) :: width
      character(len=*), intent(in) :: entry
      logical, optional, intent(in) :: end_row
      character(len=STRING_LENGTH) :: write_fmt

      write (write_fmt,"(A, I2.1, A)") "(A, A", width, ", A)"
      write (f, write_fmt, advance="no") "| ", trim(entry), " "
      if ( present(end_row) .and. end_row .eqv. .true. ) then
         write (f, "(A)", advance="yes") "|"
      end if
   end subroutine write_table_entry_character

   subroutine write_table_entry_integer( f, entry, width, end_row )
      integer, intent(in) :: f
      integer(SP), intent(in) :: entry, width
      logical, optional, intent(in) :: end_row
      character(len=STRING_LENGTH) :: write_fmt

      write (write_fmt,"(A, I2.1, A)") "(A, I", width, ", A)"
      write (f, write_fmt, advance="no") "| ", entry, " "
      if ( present(end_row) .and. end_row .eqv. .true. ) then
         write (f, "(A)", advance="yes") "|"
      end if
   end subroutine write_table_entry_integer

   subroutine write_table_entry_logical( f, entry, width, end_row )
      integer, intent(in) :: f
      integer(SP), intent(in) :: width
      logical, intent(in) :: entry
      logical, optional, intent(in) :: end_row
      character(len=STRING_LENGTH) :: write_fmt

      if ( width .lt. 5_SP ) then
         write (write_fmt,"(A, I2.1, A)") "(A, L", width, ", A)"
         write (f, write_fmt, advance="no") "| ", entry, " "
         if ( present(end_row) .and. end_row .eqv. .true. ) then
            write (f, "(A)", advance="yes") "|"
         end if
      else
         write (write_fmt,"(A, I2.1, A)") "(A, A", width, ", A)"
         if ( entry .eqv. .true. ) then
            write (f, write_fmt, advance="no") "| ", "True", " "
         else
            write (f, write_fmt, advance="no") "| ", "False", " "
         end if
         if ( present(end_row) .and. end_row .eqv. .true. ) then
            write (f, "(A)", advance="yes") "|"
         end if
      end if
   end subroutine write_table_entry_logical

   subroutine write_table_rule_entry( f, width, alignment, end_row )
      integer, intent(in) :: f
      integer(SP), intent(in) :: width
      integer(SP), optional, intent(in) :: alignment
      logical, optional, intent(in) :: end_row
      integer(SP) :: n_dash, i_dash
      logical, dimension(2) :: aligned

      aligned = (/ .false., .false. /)
      if ( present(alignment) ) then
         if      ( alignment .eq.   LEFT_ALIGNED ) then
            aligned = (/ .true., .false. /)
         else if ( alignment .eq.  RIGHT_ALIGNED ) then
            aligned = (/ .false., .true. /)
         else if ( alignment .eq. CENTER_ALIGNED ) then
            aligned = (/ .true., .true. /)
         end if
      end if
      n_dash = width - int(count(aligned),SP)

      write (f, "(A)", advance="no") "| "
      if ( aligned(1_SP) ) then
         write (f, "(A)", advance="no") ":"
      end if
      do i_dash = 1_SP, n_dash
         write (f, "(A)", advance="no") "-"
      end do
      if ( aligned(2_SP) ) then
         write (f, "(A)", advance="no") ":"
      end if
      write (f, "(A)", advance="no") " "
      if ( present(end_row) .and. end_row .eqv. .true. ) then
         write (f, "(A)", advance="yes") "|"
      end if
   end subroutine write_table_rule_entry
end module wb_text
