! Copyright (C) 2019 Andrew Trettel.  All rights reserved.

program windbag
   implicit none

   integer(4) :: number_of_arguments
   character(len=64) :: input_file_name
   logical :: input_file_exists

   number_of_arguments = command_argument_count()

   if ( number_of_arguments .eq. 0 ) then
      stop "windbag: no argument given"
   else
      call get_command_argument( 1, input_file_name )

      inquire( file=input_file_name, exist=input_file_exists )

      if ( input_file_exists ) then
         stop "windbag: input file exists"
      else
         stop "windbag: input file does not exist"
      end if
   end if
end program windbag
