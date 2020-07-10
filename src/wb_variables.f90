! Copyright (C) 2020 Andrew Trettel
! 
! Windbag is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!  ! Windbag is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with
! Windbag.  If not, see <https://www.gnu.org/licenses/>.
module wb_variables
   use wb_representation
   use wb_text
   implicit none

   private

   public WB_Variable_List, wb_variable_list_construct, &
      wb_variable_list_destroy, wb_variable_list_required_number, &
      wb_variable_list_min_required_number, &
      write_graphviz_file, write_variable_list_information, &
      wb_variable_list_add, wb_variable_list_mark_as_required, &
      wb_variable_list_set_as_minimum, wb_variable_list_add_dependency, &
      wb_variable_list_require, phase_rule

   integer(SP), public, parameter :: NUMBER_OF_PHASES        =  1_SP
   integer(SP), public, parameter :: UNUSED_VARIABLE_NUMBER  = -1_SP

   character(len=*), public, parameter :: DEFAULT_VARIABLE_NAME = "Variable"

   type, public :: WB_Variable_List
      private
      integer(SP) :: number_of_variables, max_number_of_variables, &
         min_number_of_required_variables
      logical, dimension(:), allocatable :: is_a_required_variable
      logical, dimension(:,:), allocatable :: adjacency_matrix
      integer(SP), dimension(:), allocatable :: order_of_evaluation
      character(len=STRING_LENGTH), dimension(:), allocatable :: variable_names
   end type WB_Variable_List

   interface wb_variable_list_add
      module procedure wb_variable_list_add_variable, &
         wb_variable_list_add_vector, wb_variable_list_add_tensor
   end interface
contains
   pure function phase_rule( number_of_components, number_of_phases ) &
      result( degrees_of_freedom )
      integer(SP), intent(in) :: number_of_components, number_of_phases
      integer(SP) :: degrees_of_freedom

      degrees_of_freedom = number_of_components - number_of_phases + 2_SP
   end function phase_rule

   subroutine wb_variable_list_add_dependency( vl, source_number, &
      target_number )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: source_number, target_number
      
      vl%adjacency_matrix(source_number,target_number) = .true.
   end subroutine wb_variable_list_add_dependency

   subroutine wb_variable_list_add_variable( vl, variable_name, is_required, &
      variable_number )
      type(WB_Variable_List), intent(inout) :: vl
      character(len=*), intent(in) :: variable_name
      logical, intent(in) :: is_required
      integer(SP), intent(out) :: variable_number

      variable_number = UNUSED_VARIABLE_NUMBER

      if ( wb_variable_list_number(vl) .lt. wb_variable_list_max_number(vl) ) then
         variable_number = wb_variable_list_number(vl) + 1_SP
         vl%variable_names(variable_number) = trim(variable_name)
         if ( is_required ) then
            call wb_variable_list_mark_as_required( vl, variable_number )
         end if
         vl%number_of_variables = wb_variable_list_number(vl) + 1_SP
      end if
   end subroutine wb_variable_list_add_variable

   subroutine wb_variable_list_add_tensor( vl, variable_base_name, m, n, &
      is_required, variable_numbers )
      type(WB_Variable_List), intent(inout) :: vl
      character(len=*), intent(in) :: variable_base_name
      integer(SP), intent(in) :: m, n
      logical, intent(in) :: is_required
      integer(SP), dimension(:,:), allocatable, intent(inout) :: &
         variable_numbers
      integer(SP) :: i, j, variable_number
      character(len=STRING_LENGTH) :: variable_name

      variable_numbers(:,:) = UNUSED_VARIABLE_NUMBER
      do i = 1, m
         do j = 1, n
            write (variable_name,"(A, A, I0, A, I0)") &
               trim(variable_base_name), " ", i, ",", j
            call wb_variable_list_add_variable( vl, &
               variable_name, is_required, variable_number )
            variable_numbers(i,j) = variable_number
         end do
      end do
   end subroutine wb_variable_list_add_tensor

   subroutine wb_variable_list_add_vector( vl, variable_base_name, n, &
      is_required, variable_numbers )
      type(WB_Variable_List), intent(inout) :: vl
      character(len=*), intent(in) :: variable_base_name
      integer(SP), intent(in) :: n
      logical, intent(in) :: is_required
      integer(SP), dimension(:), allocatable, intent(inout) :: &
         variable_numbers
      integer(SP) :: i, variable_number
      character(len=STRING_LENGTH) :: variable_name

      variable_numbers(:) = UNUSED_VARIABLE_NUMBER
      do i = 1, n
         write (variable_name,"(A, A, I0)") trim(variable_base_name), " ", i
         call wb_variable_list_add_variable( vl, &
            variable_name, is_required, variable_number )
         variable_numbers(i) = variable_number
      end do
   end subroutine wb_variable_list_add_vector

   subroutine wb_variable_list_construct( vl, max_number_of_variables )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: max_number_of_variables

      allocate( vl%is_a_required_variable(max_number_of_variables), &
                      vl%adjacency_matrix(max_number_of_variables,  &
                                          max_number_of_variables), &
                   vl%order_of_evaluation(max_number_of_variables), &
                        vl%variable_names(max_number_of_variables) )

      vl%number_of_variables       = 0_SP
      vl%max_number_of_variables   = max_number_of_variables
      vl%is_a_required_variable(:) = .false.
      vl%adjacency_matrix(:,:)     = .false.
      vl%order_of_evaluation(:)    =  UNUSED_VARIABLE_NUMBER
      vl%variable_names(:)         = DEFAULT_VARIABLE_NAME
   end subroutine wb_variable_list_construct

   subroutine wb_variable_list_destroy( vl )
      type(WB_Variable_List), intent(inout) :: vl

      deallocate( vl%is_a_required_variable, &
                  vl%adjacency_matrix,       &
                  vl%order_of_evaluation,    &
                  vl%variable_names )
   end subroutine wb_variable_list_destroy

   function wb_variable_list_is_dependent( vl, source_number, target_number ) &
      result( is_dependent )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: source_number, target_number
      logical :: is_dependent
      
      is_dependent = vl%adjacency_matrix(source_number,target_number)
   end function wb_variable_list_is_dependent

   function wb_variable_list_is_required( vl, variable_number ) &
      result( is_required )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: variable_number
      logical :: is_required

      is_required = vl%is_a_required_variable(variable_number)
   end function wb_variable_list_is_required

   function wb_variable_list_is_unrequired( vl, variable_number ) &
      result( is_unrequired )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: variable_number
      logical :: is_unrequired

      is_unrequired = wb_variable_list_is_required(vl,variable_number) &
         .eqv. .false.
   end function wb_variable_list_is_unrequired

   subroutine wb_variable_list_mark_as_required( vl, variable_number )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: variable_number

      if ( wb_variable_list_is_unrequired( vl, variable_number ) ) then
         vl%is_a_required_variable(variable_number) = .true.
         vl%order_of_evaluation(wb_variable_list_required_number(vl)) = &
            variable_number
      end if
   end subroutine wb_variable_list_mark_as_required

   function wb_variable_list_max_number( vl ) result( max_number_of_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: max_number_of_variables

      max_number_of_variables = vl%max_number_of_variables
   end function wb_variable_list_max_number

   function wb_variable_list_min_required_number( vl ) &
      result( min_number_of_required_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: min_number_of_required_variables

      min_number_of_required_variables = &
         vl%min_number_of_required_variables
   end function wb_variable_list_min_required_number

   function wb_variable_list_number( vl ) result( number_of_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: number_of_variables

      number_of_variables = vl%number_of_variables
   end function wb_variable_list_number

   function wb_variable_list_ordered_variable( vl, field_number ) &
      result( variable_number )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: field_number
      integer(SP) :: variable_number

      variable_number = vl%order_of_evaluation(field_number)
   end function wb_variable_list_ordered_variable

   recursive subroutine wb_variable_list_require( vl, target_number )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: target_number
      integer(SP) :: source_number

      if ( wb_variable_list_is_unrequired( vl, target_number ) ) then
         do source_number = 1, wb_variable_list_number(vl)
            if ( wb_variable_list_is_dependent( vl, &
                 source_number, target_number ) .and. &
                 wb_variable_list_is_unrequired( vl, source_number ) ) then
               call wb_variable_list_require( vl, source_number )
            end if
         end do

         call wb_variable_list_mark_as_required( vl, target_number )
      end if
   end subroutine wb_variable_list_require

   function wb_variable_list_required_number( vl ) &
      result( number_of_required_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: number_of_required_variables

      number_of_required_variables = int(count(vl%is_a_required_variable),SP)
   end function wb_variable_list_required_number

   subroutine wb_variable_list_set_as_minimum( vl )
      type(WB_Variable_List), intent(inout) :: vl

      vl%min_number_of_required_variables = &
         wb_variable_list_required_number(vl)
   end subroutine wb_variable_list_set_as_minimum

   subroutine wb_variable_list_variable_name( vl, variable_number, &
      variable_name )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: variable_number
      character(len=*), intent(inout) :: variable_name

      variable_name = vl%variable_names(variable_number)
   end subroutine wb_variable_list_variable_name

   subroutine write_graphviz_file( f, vl )
      integer, intent(in) :: f
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: source_number, target_number

      write (f, "(A)") "digraph variables {"
      write (f, "(A)") "concentrate=true"
      write (f, "(A)") "rankdir=LR"
      write (f, "(A)") 'ranksep="1.0"'
      write (f, "(A)") 'node [shape=box, fixedsize=true, width="2.0", height="0.5", margin="0.5"]'

      do source_number = 1, wb_variable_list_number(vl)
         do target_number = 1, wb_variable_list_number(vl)
            if ( wb_variable_list_is_dependent( vl, source_number, &
                                                    target_number ) ) then
               write (f, "(5A)") '"', &
                  trim(vl%variable_names(source_number)), '" -> "', &
                  trim(vl%variable_names(target_number)), '"'
            end if
         end do
      end do

      write (f, "(A)") "}"
      call write_blank_line( f )
   end subroutine write_graphviz_file

   subroutine write_variable_list_information( f, vl )
      integer, intent(in) :: f
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: i_field, i_var
      character(len=STRING_LENGTH) :: variable_name

      call write_log_heading( f, "List of fields", level=2_SP )
      call write_table_entry( f, "Field no.",    VARIABLE_COLUMN_WIDTH )
      call write_table_entry( f, "Variable no.", VARIABLE_COLUMN_WIDTH )
      call write_table_entry( f, "Name",         NAME_COLUMN_WIDTH, &
         end_row=.true. )
      call write_table_rule_entry( f, VARIABLE_COLUMN_WIDTH, &
         alignment=RIGHT_ALIGNED )
      call write_table_rule_entry( f, VARIABLE_COLUMN_WIDTH, &
         alignment=RIGHT_ALIGNED )
      call write_table_rule_entry( f, NAME_COLUMN_WIDTH, &
         alignment=LEFT_ALIGNED, end_row=.true. )
      do i_field = 1, wb_variable_list_required_number(vl)
         i_var = wb_variable_list_ordered_variable( vl, i_field )
         call wb_variable_list_variable_name( vl, i_var, variable_name )

         call write_table_entry( f, i_field, VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, i_var,   VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, variable_name, NAME_COLUMN_WIDTH, &
            end_row=.true. )
      end do
      call write_blank_line( f )
   end subroutine write_variable_list_information
end module wb_variables
