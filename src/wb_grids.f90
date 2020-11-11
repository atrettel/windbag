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

   integer(FP), public, parameter :: DEFAULT_LENGTH     = 1.0_FP
   integer(FP), public, parameter :: DEFAULT_LOCATION   = 0.0_FP
   integer(FP), public, parameter :: DEFAULT_ORIGIN     = 0.0_FP
   integer(FP), public, parameter :: DEFAULT_UNIFORMITY = 1.0_FP

   public uniform_grid, polynomial_stretched_grid
contains
   function uniform_grid( xc, xp0, lx ) &
   result( xp )
      real(FP), intent(in) :: xc, xp0, lx
      real(FP) :: xp

      xp = xp0 + xc * lx
   end function uniform_grid

   function polynomial_stretched_grid( xc, xp0, lx, xc_dmin, xc_dmax, R ) &
   result( xp )
      real(FP), intent(in) :: xc, xp0, lx, xc_dmin, xc_dmax, R
      real(FP) :: xp, k, xpf, B, quad_a, quad_b, quad_c

      xpf = xp0 + lx
      k   = lx * R

      quad_a = 1.0_FP / 5.0_FP &
             - xc_dmax + ( 4.0_FP / 3.0_FP ) * xc_dmax**2.0_FP             &
             + ( 4.0_FP / 3.0_FP ) * xc_dmin * xc_dmax                     &
             - ( 2.0_FP / 3.0_FP ) * xc_dmin**2.0_FP                       &
             + 2.0_FP * xc_dmin * xc_dmax * ( xc_dmin - 2.0_FP * xc_dmax ) &
             + xc_dmin**2.0_FP * ( xc_dmin - 2.0_FP * xc_dmax )**2.0_FP

      quad_b = k**0.5_FP * ( -2.0_FP / 3.0_FP + 2.0_FP * xc_dmax &
             + 2.0_FP * xc_dmin * ( xc_dmin - 2.0_FP * xc_dmax ) )

      quad_c = k + xp0 - xpf

      B = ( -quad_b + ( quad_b**2.0_FP - 4.0_FP * quad_a * quad_c )**0.5_FP ) &
        / ( 2.0_FP * quad_a )


      xp = xc**5.0_FP * ( ( 1.0_FP / 5.0_FP ) * B**2.0_FP )                   &
         + xc**4.0_FP * ( -B**2.0_FP * xc_dmax )                              &
         + xc**3.0_FP * ( B**2.0_FP * ( ( 4.0_FP / 3.0_FP ) * xc_dmax**2.0_FP &
                        + ( 4.0_FP / 3.0_FP ) * xc_dmin * xc_dmax             &
                        - ( 2.0_FP / 3.0_FP ) * xc_dmin**2.0_FP )             &
                        - ( 2.0_FP / 3.0_FP ) * k**0.5_FP * B )               &
         + xc**2.0_FP * ( 2.0_FP * B**2.0_FP * xc_dmin * xc_dmax *            &
                          ( xc_dmin - 2.0_FP * xc_dmax )                      &
                        + 2.0_FP * k**0.5_FP * B * xc_dmax )                  &
         + xc**1.0_FP * ( B**2.0_FP * xc_dmin**2.0_FP *                       &
                          ( xc_dmin - 2.0_FP * xc_dmax )**2.0_FP              &
                        + 2.0_FP * k**0.5_FP * B * xc_dmin *                  &
                          ( xc_dmin - 2.0_FP * xc_dmax )                      &
                        + k )                                                 &
         + xp0
   end function polynomial_stretched_grid
end module wb_grids
