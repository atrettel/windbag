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
module wb_grids
   use wb_representation
   implicit none

   private

   public uniform_grid
contains
   function uniform_grid( xc, xp0, lx ) &
   result( xp )
      real(FP), intent(in) :: xc, xp0, lx
      real(FP) :: xp

      xp = xp0 + xc * lx
   end function uniform_grid
end module wb_grids
