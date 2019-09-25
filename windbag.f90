! Copyright (C) 2019 Andrew Trettel.  All rights reserved.
program windbag
   use wbbase

   implicit none
   character(len=STRING_LENGTH) :: input_file_name, casename
   type(WB_Field_Data) :: field_data
   integer(IP) :: number_of_components, number_of_ghost_points
   integer(IP) :: nx_global, ny_global, nz_global
   integer(IP), dimension(:), allocatable :: factors

   call boot_program( input_file_name )

   call read_general_namelist( input_file_name, casename, &
      number_of_components, number_of_ghost_points, &
      nx_global, ny_global, nz_global )

   field_data = WB_Field_Data( number_of_components, number_of_ghost_points, &
      nx_global, ny_global, nz_global )

   if ( field_data%i_proc .eq. ROOT_PROCESS_NUMBER ) then
      write (output_unit, "(A, A, A)") "'", trim(input_file_name), "'"
      write (output_unit, "(A, A, A)") "'", trim(casename),        "'"
      write (output_unit, "(A, I3)") "n_proc                 = ", field_data%n_proc
      write (output_unit, "(A, I3)") "i_proc                 = ", field_data%i_proc
      write (output_unit, "(A, I3)") "number_of_components   = ", field_data%number_of_components
      write (output_unit, "(A, I3)") "number_of_ghost_points = ", field_data%number_of_ghost_points
      write (output_unit, "(A, I3)") "nx_global              = ", field_data%nx_global
      write (output_unit, "(A, I3)") "ny_global              = ", field_data%ny_global
      write (output_unit, "(A, I3)") "nz_global              = ", field_data%nz_global

      call integer_factorization( 6_IP, factors )
      print *, factors
      deallocate( factors )

      call integer_factorization( 30_IP, factors )
      print *, factors
      deallocate( factors )
   end if

   call stop_program( "end of program", EXIT_SUCCESS )
end program windbag
