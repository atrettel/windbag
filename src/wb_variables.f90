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
module wb_variables
   use wb_representation
   implicit none

   private

   public WB_Variable_List, wb_variable_list_construct, &
      wb_variable_list_destroy

   integer(SP), public, parameter :: MAX_NUMBER_OF_VARIABLES =  128_SP
   integer(SP), public, parameter :: VACANT_VARIABLE_NUMBER  =   -1_SP

   character(len=*), public, parameter :: DEFAULT_VARIABLE_NAME = "Variable"

   type, public :: WB_Variable_List
      private
      integer(SP) :: next_variable_number
      logical, dimension(:), allocatable :: is_a_required_variable
      logical, dimension(:,:), allocatable :: adjacency_matrix
      integer(SP), dimension(:), allocatable :: order_of_evaluation
      character(len=STRING_LENGTH), dimension(:), allocatable :: variable_names
   end type WB_Variable_List
contains
   subroutine wb_variable_list_construct( vl )
      type(WB_Variable_List) :: vl

      allocate( vl%is_a_required_variable(MAX_NUMBER_OF_VARIABLES), &
                      vl%adjacency_matrix(MAX_NUMBER_OF_VARIABLES,  &
                                          MAX_NUMBER_OF_VARIABLES), &
                   vl%order_of_evaluation(MAX_NUMBER_OF_VARIABLES), &
                        vl%variable_names(MAX_NUMBER_OF_VARIABLES) )

      vl%next_variable_number      = 1_SP
      vl%is_a_required_variable(:) = .false.
      vl%adjacency_matrix(:,:)     = .false.
      vl%order_of_evaluation(:)    = VACANT_VARIABLE_NUMBER
      vl%variable_names(:)         = DEFAULT_VARIABLE_NAME
   end subroutine wb_variable_list_construct

   subroutine wb_variable_list_destroy( vl )
      type(WB_Variable_List) :: vl

      deallocate( vl%is_a_required_variable, &
                  vl%adjacency_matrix,       &
                  vl%order_of_evaluation,    &
                  vl%variable_names )
   end subroutine wb_variable_list_destroy
end module wb_variables
