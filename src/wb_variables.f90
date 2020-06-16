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
      wb_variable_list_destroy, wb_variable_list_required_number, &
      construct_compressible_conservative_variables

   integer(SP), public, parameter :: NUMBER_OF_PHASES        =  1_SP
   integer(SP), public, parameter :: UNUSED_VARIABLE_NUMBER  = -1_SP

   character(len=*), public, parameter :: DEFAULT_VARIABLE_NAME = "Variable"

   type, public :: WB_Variable_List
      private
      integer(SP) :: number_of_variables, max_number_of_variables
      logical, dimension(:), allocatable :: is_a_required_variable
      logical, dimension(:,:), allocatable :: adjacency_matrix
      integer(SP), dimension(:), allocatable :: order_of_evaluation
      character(len=STRING_LENGTH), dimension(:), allocatable :: variable_names
   end type WB_Variable_List

   interface wb_variable_list_add
      module procedure wb_variable_list_add_variable, &
         wb_variable_list_add_vector
   end interface
contains
   subroutine construct_compressible_conservative_variables( vl, nd, nc )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP) :: nc, nd
      integer(SP) :: l_bulk_viscosity,                   &
                     l_dilatational_viscosity,           &
                     l_dynamic_viscosity,                &
                     l_heat_capacity_ratio,              &
                     l_kinematic_viscosity,              &
                     l_mach_number,                      &
                     l_mass_density,                     &
                     l_prandtl_number,                   &
                     l_pressure,                         &
                     l_specific_enthalpy,                &
                     l_specific_entropy,                 &
                     l_specific_internal_energy,         &
                     l_specific_isobaric_heat_capacity,  &
                     l_specific_isochoric_heat_capacity, &
                     l_specific_total_enthalpy,          &
                     l_specific_total_internal_energy,   &
                     l_specific_volume,                  &
                     l_speed,                            &
                     l_speed_of_sound,                   &
                     l_temperature,                      &
                     l_thermal_conductivity,             &
                     l_thermal_diffusivity
      integer(SP), dimension(:), allocatable :: &
         l_amount_fractions,   &
         l_coordinates,        &
         l_mass_fractions,     &
         l_momentum_densities, &
         l_velocities
      integer(SP) :: i_dim, ic, thermodynamic_degrees_of_freedom, &
         number_of_required_mass_fractions
      logical :: density_is_required, energy_is_required

      allocate( l_amount_fractions(nc), &
                     l_coordinates(nd), &
                  l_mass_fractions(nc), &
              l_momentum_densities(nd), &
                      l_velocities(nd)  )

      thermodynamic_degrees_of_freedom = phase_rule( nc, NUMBER_OF_PHASES )

      density_is_required               = .false.
      energy_is_required                = .false.
      number_of_required_mass_fractions = 0_SP
      if ( thermodynamic_degrees_of_freedom .gt. 0_SP ) then
         density_is_required = .true.
      end if
      if ( thermodynamic_degrees_of_freedom .gt. 1_SP ) then
         energy_is_required = .true.
      end if
      if ( thermodynamic_degrees_of_freedom .gt. 2_SP ) then
         number_of_required_mass_fractions = &
            thermodynamic_degrees_of_freedom - 2_SP
      end if

      ! Declarations
      call wb_variable_list_add( vl, "Amount fraction",                       nc, .false., l_amount_fractions                 )
      call wb_variable_list_add( vl, "Bulk viscosity",                            .false., l_bulk_viscosity                   )
      call wb_variable_list_add( vl, "Coordinate",                            nd,  .true., l_coordinates                      )
      call wb_variable_list_add( vl, "Dilatational viscosity",                    .false., l_dilatational_viscosity           )
      call wb_variable_list_add( vl, "Dynamic viscosity",                         .false., l_dynamic_viscosity                )
      call wb_variable_list_add( vl, "Heat capacity ratio",                       .false., l_heat_capacity_ratio              )
      call wb_variable_list_add( vl, "Kinematic viscosity",                       .false., l_kinematic_viscosity              )
      call wb_variable_list_add( vl, "Mach number",                               .false., l_mach_number                      )
      call wb_variable_list_add( vl, "Mass density",                  density_is_required, l_mass_density                     )
      call wb_variable_list_add( vl, "Mass fraction",                         nc, .false., l_mass_fractions                   )
      call wb_variable_list_add( vl, "Momentum density",                      nd,  .true., l_momentum_densities               )
      call wb_variable_list_add( vl, "Prandtl number",                            .false., l_prandtl_number                   )
      call wb_variable_list_add( vl, "Pressure",                                  .false., l_pressure                         )
      call wb_variable_list_add( vl, "Specific enthalpy",                         .false., l_specific_enthalpy                )
      call wb_variable_list_add( vl, "Specific entropy",                          .false., l_specific_entropy                 )
      call wb_variable_list_add( vl, "Specific internal energy",                  .false., l_specific_internal_energy         )
      call wb_variable_list_add( vl, "Specific isobaric heat capacity",           .false., l_specific_isobaric_heat_capacity  )
      call wb_variable_list_add( vl, "Specific isochoric heat capacity",          .false., l_specific_isochoric_heat_capacity )
      call wb_variable_list_add( vl, "Specific total enthalpy",                   .false., l_specific_total_enthalpy          )
      call wb_variable_list_add( vl, "Specific total internal energy", energy_is_required, l_specific_total_internal_energy   )
      call wb_variable_list_add( vl, "Specific volume",                           .false., l_specific_volume                  )
      call wb_variable_list_add( vl, "Speed",                                     .false., l_speed                            )
      call wb_variable_list_add( vl, "Speed of sound",                            .false., l_speed_of_sound                   )
      call wb_variable_list_add( vl, "Temperature",                               .false., l_temperature                      )
      call wb_variable_list_add( vl, "Thermal conductivity",                      .false., l_thermal_conductivity             )
      call wb_variable_list_add( vl, "Thermal diffusivity",                       .false., l_thermal_diffusivity              )
      call wb_variable_list_add( vl, "Velocity",                              nd, .false., l_velocities                       )

      do ic = 1, number_of_required_mass_fractions
         call wb_variable_list_mark_as_required( vl, l_mass_fractions(ic) )
      end do

      ! Dependencies

      ! gamma = c_p / c_v
      call wb_variable_list_add_dependency( vl, l_specific_isobaric_heat_capacity,  l_heat_capacity_ratio )
      call wb_variable_list_add_dependency( vl, l_specific_isochoric_heat_capacity, l_heat_capacity_ratio )
      ! nu = mu / rho
      call wb_variable_list_add_dependency( vl, l_dynamic_viscosity, l_kinematic_viscosity )
      call wb_variable_list_add_dependency( vl, l_mass_density,      l_kinematic_viscosity )
      ! M = V / a
      call wb_variable_list_add_dependency( vl, l_speed,          l_mach_number )
      call wb_variable_list_add_dependency( vl, l_speed_of_sound, l_mach_number )
      ! Pr = nu / alpha
      call wb_variable_list_add_dependency( vl, l_kinematic_viscosity, l_prandtl_number )
      call wb_variable_list_add_dependency( vl, l_thermal_diffusivity, l_prandtl_number )
      ! v = 1 / rho
      call wb_variable_list_add_dependency( vl, l_mass_density, l_specific_volume )
      ! V = sqrt( u^2 + v^2 + w^2 )
      do i_dim = 1, nd
         call wb_variable_list_add_dependency( vl, l_velocities(i_dim), l_speed )
      end do
      ! alpha = k / ( c_p rho )
      call wb_variable_list_add_dependency( vl, l_mass_density,                    l_thermal_diffusivity )
      call wb_variable_list_add_dependency( vl, l_specific_isobaric_heat_capacity, l_thermal_diffusivity )
      call wb_variable_list_add_dependency( vl, l_thermal_conductivity,            l_thermal_diffusivity )
      ! u = (rho u) / rho
      do i_dim = 1, nd
         call wb_variable_list_add_dependency( vl, l_mass_density,              l_velocities(i_dim) )
         call wb_variable_list_add_dependency( vl, l_momentum_densities(i_dim), l_velocities(i_dim) )
      end do

      ! Requirements
      call wb_variable_list_require( vl, l_speed )

      deallocate( l_amount_fractions,   &
                  l_coordinates,        &
                  l_mass_fractions,     &
                  l_momentum_densities, &
                  l_velocities          )
   end subroutine construct_compressible_conservative_variables

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
         write (variable_name,"(A, A, I2.1)") trim(variable_base_name), " ", i
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

   function wb_variable_list_is_unrequired( vl, variable_number ) &
      result( is_unrequired )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP), intent(in) :: variable_number
      logical :: is_unrequired
      
      is_unrequired = vl%is_a_required_variable(variable_number) .eqv. .false.
   end function wb_variable_list_is_unrequired

   subroutine wb_variable_list_mark_as_required( vl, variable_number )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: variable_number

      if ( vl%is_a_required_variable(variable_number) .eqv. .false. ) then
         vl%is_a_required_variable(variable_number) = .true.
         vl%order_of_evaluation(wb_variable_list_required_number(vl)+1_SP) = &
            variable_number
      end if
   end subroutine wb_variable_list_mark_as_required

   function wb_variable_list_max_number( vl ) result( max_number_of_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: max_number_of_variables

      max_number_of_variables = vl%max_number_of_variables
   end function wb_variable_list_max_number

   function wb_variable_list_number( vl ) result( number_of_variables )
      type(WB_Variable_List), intent(in) :: vl
      integer(SP) :: number_of_variables

      number_of_variables = vl%number_of_variables
   end function wb_variable_list_number

   recursive subroutine wb_variable_list_require( vl, target_number )
      type(WB_Variable_List), intent(inout) :: vl
      integer(SP), intent(in) :: target_number
      integer(SP) :: source_number

      if ( vl%is_a_required_variable(target_number) .eqv. .false. ) then
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
end module wb_variables
