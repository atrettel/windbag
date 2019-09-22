! Copyright (C) 2019 Andrew Trettel.  All rights reserved.
program windbag
   use wbbase

   implicit none
   character(len=STRING_LENGTH) :: input_file_name
   type(WB_Field_Data) :: field_data

   call boot_program( input_file_name )

   field_data = WB_Field_Data( 128_IP, 128_IP, 128_IP )

   if ( field_data%i_proc .eq. ROOT_PROCESS_NUMBER ) then
      print *, input_file_name
      print *, field_data%n_proc
      print *, field_data%i_proc
      print *, field_data%nx_global
      print *, field_data%ny_global
      print *, field_data%nz_global
   end if

   call stop_program( "end of program", EXIT_SUCCESS )
end program windbag
