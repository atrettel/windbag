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
   use wb_text
   implicit none

   private

   public phase_rule,                        &
          wb_variable_list_add,              &
          wb_variable_list_add_dependency,   &
          wb_variable_list_construct,        &
          wb_variable_list_destroy,          &
          wb_variable_list_mark_as_required, &
          wb_variable_list_min_required,     &
          wb_variable_list_require,          &
          wb_variable_list_sequence_index,   &
          wb_variable_list_set_as_minimum,   &
          wb_variable_list_total_required,   &
          wb_variable_list_variable_id,      &
          wb_variable_list_variable_name,    &
          write_graphviz_file,               &
          write_variable_list_information

   integer(SP), public, parameter :: NO_SEQUENCE_INDEX  = -1_SP
   integer(SP), public, parameter :: NUMBER_OF_PHASES   =  1_SP
   integer(SP), public, parameter :: UNUSED_VARIABLE_ID = -1_SP

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

   subroutine wb_variable_list_add_dependency( vl, source_id, target_id )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: source_id, target_id
      
      vl%adjacency_matrix(source_id,target_id) = .true.
   end subroutine wb_variable_list_add_dependency

   subroutine wb_variable_list_add_variable( vl, variable_name, is_required, &
      variable_id )
      type(WB_Variable_List), intent(inout) :: vl
      character(len=*), intent(in) :: variable_name
      logical, intent(in) :: is_required
      integer(SP), intent(out) :: variable_id

      variable_id = UNUSED_VARIABLE_ID

      if ( wb_variable_list_total_variables(vl) .lt. wb_variable_list_max_variables(vl) ) then
         variable_id = wb_variable_list_total_variables(vl) + 1_SP
         vl%variable_names(variable_id) = trim(variable_name)
         if ( is_required ) then
            call wb_variable_list_mark_as_required( vl, variable_id )
         end if
         vl%number_of_variables = wb_variable_list_total_variables(vl) + 1_SP
      end if
   end subroutine wb_variable_list_add_variable

   subroutine wb_variable_list_add_tensor( vl, variable_base_name, m, n, &
      is_required, variable_ids )
      type(WB_Variable_List), intent(inout) :: vl
      character(len=*), intent(in) :: variable_base_name
      integer(SP), intent(in) :: m, n
      logical, intent(in) :: is_required
      integer(SP), dimension(:,:), allocatable, intent(inout) :: &
         variable_ids
      integer(SP) :: i, j, variable_id
      character(len=STRING_LENGTH) :: variable_name

      variable_ids(:,:) = UNUSED_VARIABLE_ID
      do i = 1, m
         do j = 1, n
            write (variable_name,"(A, A, I0, A, I0)") &
               trim(variable_base_name), " ", i, ",", j
            call wb_variable_list_add_variable( vl, &
               variable_name, is_required, variable_id )
            variable_ids(i,j) = variable_id
         end do
      end do
   end subroutine wb_variable_list_add_tensor

   subroutine wb_variable_list_add_vector( vl, variable_base_name, n, &
      is_required, variable_ids )
      type(WB_Variable_List), intent(inout) :: vl
      character(len=*), intent(in) :: variable_base_name
      integer(SP), intent(in) :: n
      logical, intent(in) :: is_required
      integer(SP), dimension(:), allocatable, intent(inout) :: &
         variable_ids
      integer(SP) :: i, variable_id
      character(len=STRING_LENGTH) :: variable_name

      variable_ids(:) = UNUSED_VARIABLE_ID
      do i = 1, n
         write (variable_name,"(A, A, I0)") trim(variable_base_name), " ", i
         call wb_variable_list_add_variable( vl, &
            variable_name, is_required, variable_id )
         variable_ids(i) = variable_id
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
      vl%order_of_evaluation(:)    =  UNUSED_VARIABLE_ID
      vl%variable_names(:)         = DEFAULT_VARIABLE_NAME
   end subroutine wb_variable_list_construct

   subroutine wb_variable_list_destroy( vl )
      type(WB_Variable_List), intent(inout) :: vl

      deallocate( vl%is_a_required_variable, &
                  vl%adjacency_matrix,       &
                  vl%order_of_evaluation,    &
                  vl%variable_names )
   end subroutine wb_variable_list_destroy

   function wb_variable_list_is_dependent( vl, source_id, target_id ) &
   result( is_dependent )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: source_id, target_id
      logical :: is_dependent
      
      is_dependent = vl%adjacency_matrix(source_id,target_id)
   end function wb_variable_list_is_dependent

   function wb_variable_list_is_required( vl, variable_id ) &
   result( is_required )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: variable_id
      logical :: is_required

      is_required = vl%is_a_required_variable(variable_id)
   end function wb_variable_list_is_required

   function wb_variable_list_is_unrequired( vl, variable_id ) &
   result( is_unrequired )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: variable_id
      logical :: is_unrequired

      is_unrequired = wb_variable_list_is_required(vl,variable_id) &
         .eqv. .false.
   end function wb_variable_list_is_unrequired

   subroutine wb_variable_list_mark_as_required( vl, variable_id )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: variable_id

      if ( wb_variable_list_is_unrequired( vl, variable_id ) ) then
         vl%is_a_required_variable(variable_id) = .true.
         vl%order_of_evaluation(wb_variable_list_total_required(vl)) = &
            variable_id
      end if
   end subroutine wb_variable_list_mark_as_required

   function wb_variable_list_max_variables( vl ) &
   result( max_number_of_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: max_number_of_variables

      max_number_of_variables = vl%max_number_of_variables
   end function wb_variable_list_max_variables

   function wb_variable_list_min_required( vl ) &
   result( min_required )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: min_required

      min_required = vl%min_number_of_required_variables
   end function wb_variable_list_min_required

   function wb_variable_list_total_variables( vl ) &
   result( total_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: total_variables

      total_variables = vl%number_of_variables
   end function wb_variable_list_total_variables

   function wb_variable_list_sequence_index( vl, variable_id ) &
   result( sequence_index )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: variable_id
      integer(SP) :: sequence_index, loop_sequence_index

      sequence_index = NO_SEQUENCE_INDEX
      do loop_sequence_index = 1, wb_variable_list_total_required(vl)
         if ( vl%order_of_evaluation(loop_sequence_index) .eq. &
              variable_id ) then
            sequence_index = loop_sequence_index
         end if
      end do
   end function wb_variable_list_sequence_index

   function wb_variable_list_variable_id( vl, sequence_index ) &
   result( variable_id )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: sequence_index
      integer(SP) :: variable_id

      if ( sequence_index .lt. 1_SP .or. &
           sequence_index .gt. wb_variable_list_total_required(vl) ) then
         variable_id = UNUSED_VARIABLE_ID
      else
         variable_id = vl%order_of_evaluation(sequence_index)
      end if
   end function wb_variable_list_variable_id

   recursive subroutine wb_variable_list_require( vl, target_id )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: target_id
      integer(SP) :: source_id

      if ( wb_variable_list_is_unrequired( vl, target_id ) ) then
         do source_id = 1, wb_variable_list_total_variables(vl)
            if ( wb_variable_list_is_dependent( vl, &
                 source_id, target_id ) .and. &
                 wb_variable_list_is_unrequired( vl, source_id ) ) then
               call wb_variable_list_require( vl, source_id )
            end if
         end do

         call wb_variable_list_mark_as_required( vl, target_id )
      end if
   end subroutine wb_variable_list_require

   function wb_variable_list_total_required( vl ) &
   result( number_of_required_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: number_of_required_variables

      number_of_required_variables = int(count(vl%is_a_required_variable),SP)
   end function wb_variable_list_total_required

   subroutine wb_variable_list_set_as_minimum( vl )
      type(WB_Variable_List), intent(inout) :: vl

      vl%min_number_of_required_variables = wb_variable_list_total_required(vl)
   end subroutine wb_variable_list_set_as_minimum

   subroutine wb_variable_list_variable_name( vl, variable_id, variable_name )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: variable_id
      character(len=*), intent(inout) :: variable_name

      variable_name = vl%variable_names(variable_id)
   end subroutine wb_variable_list_variable_name

   subroutine write_graphviz_file( f, vl )
      integer, intent(in) :: f
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: source_id, target_id

      write (f, "(A)") "digraph variables {"
      write (f, "(A)") "concentrate=true"
      write (f, "(A)") "rankdir=LR"
      write (f, "(A)") 'ranksep="1.0"'
      write (f, "(A)") 'node [shape=box, fixedsize=true, width="2.0", height="0.5", margin="0.5"]'

      do source_id = 1, wb_variable_list_total_variables(vl)
         do target_id = 1, wb_variable_list_total_variables(vl)
            if ( wb_variable_list_is_dependent( vl, source_id, &
                                                    target_id ) ) then
               write (f, "(5A)") '"', &
                  trim(vl%variable_names(source_id)), '" -> "', &
                  trim(vl%variable_names(target_id)), '"'
            end if
         end do
      end do

      write (f, "(A)") "}"
      call write_blank_line( f )
   end subroutine write_graphviz_file

   subroutine write_variable_list_information( f, vl )
      integer, intent(in) :: f
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: sequence_index, variable_id
      character(len=STRING_LENGTH) :: variable_name

      call write_log_heading( f, "List of required variables", level=2_SP )
      call write_table_entry( f, "Sequence index", VARIABLE_COLUMN_WIDTH )
      call write_table_entry( f, "Variable ID",    VARIABLE_COLUMN_WIDTH )
      call write_table_entry( f, "Name",         NAME_COLUMN_WIDTH, &
         end_row=.true. )
      call write_table_rule_entry( f, VARIABLE_COLUMN_WIDTH, &
         alignment=RIGHT_ALIGNED )
      call write_table_rule_entry( f, VARIABLE_COLUMN_WIDTH, &
         alignment=RIGHT_ALIGNED )
      call write_table_rule_entry( f, NAME_COLUMN_WIDTH, &
         alignment=LEFT_ALIGNED, end_row=.true. )
      do sequence_index = 1, wb_variable_list_total_required(vl)
         variable_id = wb_variable_list_variable_id( vl, sequence_index )
         call wb_variable_list_variable_name( vl, variable_id, variable_name )

         call write_table_entry( f, sequence_index, VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, variable_id,    VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, variable_name,      NAME_COLUMN_WIDTH, &
            end_row=.true. )
      end do
      call write_blank_line( f )
   end subroutine write_variable_list_information
end module wb_variables
