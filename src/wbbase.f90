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
module wbbase
   use iso_fortran_env, only : compiler_options, compiler_version, real64
   use mpi_f08
   implicit none

   private

   public WB_State, check_input_file, deallocate_state, initialize_state, &
      print_initial_information

   integer, public, parameter ::            FP = real64
   integer, public, parameter ::            ND = 3
   integer, public, parameter ::  BLOCK_MASTER = 0
   integer, public, parameter ::  WORLD_MASTER = 0
   integer, public, parameter :: STRING_LENGTH = 64

   character(len=*), public, parameter :: PROGRAM_NAME = "windbag"
   character(len=*), public, parameter ::      VERSION = "0.0.0"

   type(MPI_Datatype), public, save :: MPI_FP

   type WB_Block
      integer :: block_size
      integer, dimension(ND) :: np, nx
      integer, dimension(ND,2) :: neighbors
      logical, dimension(ND) :: periods
      logical :: reorder = .false.
   end type WB_Block

   type WB_Process
      integer :: block_rank
      integer, dimension(ND) :: block_coords, nx
      integer :: ib
   end type WB_Process

   type WB_State
      character(len=:), allocatable :: case_name
      integer :: block_rank, block_size
      integer :: world_rank, world_size
      integer :: ib, nb
      integer :: nf
      integer :: ng
      integer :: nv = 5
      integer, dimension(ND) :: nx
      integer, dimension(ND) :: block_coords
      integer, dimension(ND,2) :: neighbors
      real(FP) :: t = 0.0_FP
      real(FP), dimension(:,:,:,:), allocatable :: f
      type(MPI_Comm) :: comm_block
      type(WB_Block), dimension(:), allocatable :: blocks
      type(WB_Process), dimension(:), allocatable :: processes
   end type WB_State
contains
   subroutine check_block_neighbors( s )
      type(WB_State), intent(in) :: s
      integer :: ib, id, id2, ierr, neighbor_l, neighbor_u

      if ( s%world_rank .eq. WORLD_MASTER ) then
         ! Ensure that lower and upper pairs exist, and that their number of
         ! points and processes match.
         do ib = 1, s%nb
            do id = 1, ND
               neighbor_l = s%blocks(ib)%neighbors(id,1)
               if ( neighbor_l .ne. 0 ) then
                  neighbor_u = s%blocks(neighbor_l)%neighbors(id,2)
                  if ( ib .ne. neighbor_u ) then
                     write (*,"(A, A, I1, A, I1, A, I1)") &
                        PROGRAM_NAME, ": lower face of block ", ib, &
                        " does not neighbor upper face of block ", &
                        neighbor_l, " in direction ", id
                     call mpi_abort( MPI_COMM_WORLD, MPI_ERR_TOPOLOGY, ierr )
                  else
                     do id2 = 1, ND
                        if ( id2 .ne. id .and. s%blocks(ib)%np(id2) .ne. &
                           s%blocks(neighbor_l)%np(id2) ) then
                           write (*,"(A, A, I1, A, I1, A, I1, A, I1)") &
                              PROGRAM_NAME, ": face in direction ", id, &
                              " shared by blocks ", ib, " and ", neighbor_l, &
                              " does not match processes in direction ", id2
                           call mpi_abort( MPI_COMM_WORLD, MPI_ERR_TOPOLOGY, &
                              ierr )
                        end if
                        if ( id2 .ne. id .and. s%blocks(ib)%nx(id2) .ne. &
                           s%blocks(neighbor_l)%nx(id2) ) then
                           write (*,"(A, A, I1, A, I1, A, I1, A, I1)") &
                              PROGRAM_NAME, ": face in direction ", id, &
                              " shared by blocks ", ib, " and ", neighbor_l, &
                              " does not match points in direction ", id2
                           call mpi_abort( MPI_COMM_WORLD, MPI_ERR_TOPOLOGY, &
                              ierr )
                        end if
                     end do
                  end if
               end if
            end do
         end do
      end if
   end subroutine check_block_neighbors

   subroutine check_input_file( filename )
      character(len=STRING_LENGTH), intent(out)  :: filename
      integer :: argc, filename_length, ierr, world_rank
      logical :: file_exists
      call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr )
      call find_mpi_fp
      if ( world_rank .eq. WORLD_MASTER ) then
         argc = command_argument_count()
         if ( argc .eq. 0 ) then
            write (*,"(A, A, A)") "Usage: ", PROGRAM_NAME, " [INPUT_FILE]"
            call mpi_abort( MPI_COMM_WORLD, MPI_SUCCESS, ierr )
         end if
         call get_command_argument( 1, filename, filename_length, ierr )
         inquire( file=filename, exist=file_exists, iostat=ierr )
         if ( file_exists .eqv. .false. ) then
            write (*,"(A, A)") PROGRAM_NAME, ": input file does not exist"
            call mpi_abort( MPI_COMM_WORLD, MPI_ERR_NO_SUCH_FILE, ierr )
         end if
      end if
      call mpi_barrier( MPI_COMM_WORLD, ierr )
   end subroutine check_input_file

   subroutine deallocate_state( s )
      integer :: ierr
      type(WB_State), intent(inout) :: s

      deallocate( s%case_name )
      deallocate( s%blocks )
      deallocate( s%processes )
      call mpi_comm_free( s%comm_block, ierr )
   end subroutine deallocate_state

   subroutine find_mpi_fp
      integer :: mpi_float_size, ierr
      call mpi_sizeof( 1.0_FP, mpi_float_size, ierr )
      call mpi_type_match_size( MPI_TYPECLASS_REAL, mpi_float_size, MPI_FP, &
         ierr )
   end subroutine find_mpi_fp

   subroutine identify_process_neighbors( s )
      integer :: id, idir, ierr, block_neighbor, world_rank
      integer, dimension(2) :: block_ranks
      integer, dimension(ND) :: block_coords
      type(WB_State), intent(inout) :: s

      do id = 1, ND
         call mpi_cart_shift( s%comm_block, id-1, 1, &
            block_ranks(1), block_ranks(2), ierr )
         do idir = 1, 2
            if ( block_ranks(idir) .ne. MPI_PROC_NULL ) then
               do world_rank = 0, s%world_size-1
                  if ( s%processes(world_rank)%ib .eq. s%ib .and. &
                     s%processes(world_rank)%block_rank .eq. &
                     block_ranks(idir) ) then
                     exit
                  end if
               end do
               s%neighbors(id,idir) = world_rank
            else
               block_neighbor = s%blocks(s%ib)%neighbors(id,idir)
               if ( block_neighbor .eq. 0 ) then
                  s%neighbors(id,idir) = MPI_PROC_NULL
               else
                  ! This block neighbors another block.  This neighboring block
                  ! sits on the opposite side of the current block.  Both
                  ! processes share the same block coordinates except for the
                  ! current direction id.  The block coordinates for the id
                  ! direction would be opposites for each.
                  block_coords = s%block_coords
                  if ( idir .eq. 1 ) then
                     block_coords(id) = s%blocks(block_neighbor)%np(id)-1
                  else
                     block_coords(id) = 0
                  end if
                  do world_rank = 0, s%world_size-1
                     if ( s%processes(world_rank)%ib .eq. block_neighbor &
                        .and. all( s%processes(world_rank)%block_coords &
                        .eq. block_coords ) ) then
                        exit
                     end if
                  end do
                  s%neighbors(id,idir) = world_rank
               end if
            end if
         end do
      end do
   end subroutine identify_process_neighbors

   subroutine initialize_state( s, filename )
      character(len=STRING_LENGTH) :: filename
      integer :: ierr
      type(WB_State), intent(out) :: s

      call mpi_comm_rank( MPI_COMM_WORLD, s%world_rank, ierr )
      call mpi_comm_size( MPI_COMM_WORLD, s%world_size, ierr )
      call read_general_namelist( s, filename )
      call read_block_namelists( s, filename )
      call check_block_neighbors( s )
      call setup_processes( s )
   end subroutine initialize_state

   subroutine print_block_information( s )
      integer :: ib, id
      type(WB_State), intent(in) :: s

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") "## Block information"
         write (*,"(A)") ""

         write (*,"(A)", advance="no") "| `ib` |    size | `nd`           "
         write (*,"(A)", advance="yes") "|       points | `nx`              |"

         write (*,"(A)", advance="no") "| ---: | ------: | :------------- "
         write (*,"(A)", advance="yes") "| -----------: | :---------------- |"

         do ib = 1, s%nb
            write (*,"(A, I4, A)", advance="no") "| ", ib, " "
            write (*,"(A, I7, A)", advance="no") "| ", &
               s%blocks(ib)%block_size, " | ("
            do id = 1, ND
               write (*,"(I3, A)", advance="no") s%blocks(ib)%np(id), ","
            end do
            write (*,"(A, I12, A)", advance="no") ") | ", &
               product(s%blocks(ib)%nx), " | ("
            do id = 1, ND
               write (*,"(I4, A)", advance="no") s%blocks(ib)%nx(id), ","
            end do
            write (*,"(A)") ") |"
         end do
         write (*,"(A)") ""
      end if
   end subroutine

   subroutine print_initial_information( s )
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: lib_version
      integer :: ierr, mpi_major_version_number, mpi_minor_version_number, &
         string_length
      type(WB_State), intent(in) :: s

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A, A, A, A, A, A, A)") "# ", PROGRAM_NAME, " ", VERSION, &
            ", case `", s%case_name, "`"
         write (*,"(A)") ""

         call mpi_get_version( mpi_major_version_number, &
            mpi_minor_version_number, ierr )
         write (*,"(A, I1, A, I1)") "- MPI version: ", &
            mpi_major_version_number, ".", mpi_minor_version_number
         write (*,"(A)") ""

         call mpi_get_library_version( lib_version, string_length, ierr )
         write (*,"(A, A)") "- MPI library version: ", trim(lib_version)
         write (*,"(A)") ""

         write (*,"(A, A, A, A, A)") "- Compiled using ", compiler_version(), &
            " using the following options: `", compiler_options(), "`"
         write (*,"(A)") ""

         if ( MPI_SUBARRAYS_SUPPORTED ) then
            write (*,"(A)") "- MPI subarrays are supported."
         else
            write (*,"(A)") "- MPI subarrays are not supported."
         end if
         write (*,"(A)") ""
      end if
      call print_block_information( s )
      call print_process_information( s )
      call print_process_neighbors( s )
   end subroutine print_initial_information

   subroutine print_process_information( s )
      integer :: id, ierr, string_length, world_rank
      character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name
      type(WB_State), intent(in) :: s

      call mpi_get_processor_name( processor_name, string_length, ierr )

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") "## Process information"
         write (*,"(A)") ""

         write (*,"(A)", advance="no") &
            "| `world_rank` | hostname     | `ib` | `block_rank` "
         write (*,"(A)", advance="yes") &
            "| `block_coords` |    points |           `nx` |"

         write (*,"(A)", advance="no") &
            "| -----------: | :----------- | ---: | -----------: "
         write (*,"(A)", advance="yes") &
            "| :------------- | --------: | :------------- |"
      end if

      do world_rank = 0, s%world_size-1
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            write (*,"(A, I12, A)", advance="no") "| ", s%world_rank, " "
            write (*,"(A, A12, A)", advance="no") "| ", &
               trim(processor_name), " "
            write (*,"(A, I4, A)", advance="no") "| ", &
               s%ib, " "
            write (*,"(A, I12, A)", advance="no") "| ", &
               s%block_rank, " "
            write (*,"(A)", advance="no") "| ("
            do id = 1, ND
               write (*,"(I3, A)", advance="no") &
                  s%block_coords(id), ","
            end do
            write (*,"(A, I9, A)", advance="no") ") | ", &
               product(s%nx), " | ("
            do id = 1, ND
               write (*,"(I3, A)", advance="no") &
                  s%nx(id), ","
            end do
            write (*,"(A)") ") |"
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") ""
      end if
   end subroutine print_process_information

   subroutine print_process_neighbors( s )
      integer :: id, idir, ierr, world_rank
      type(WB_State), intent(in) :: s

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") "## Process neighbors"
         write (*,"(A)") ""

         write (*,"(A)", advance="no") "| `world_rank` "
         do id = 1, ND
            write (*,"(A, I1, A)", advance="no") "|       ", id, "L "
            write (*,"(A, I1, A)", advance="no") "|       ", id, "U "
         end do
         write (*,"(A)") "|"

         write (*,"(A)", advance="no") "| -----------: "
         do id = 1, ND
            write (*,"(A, I1, A)", advance="no") "| -------: "
            write (*,"(A, I1, A)", advance="no") "| -------: "
         end do
         write (*,"(A)") "|"
      end if

      do world_rank = 0, s%world_size-1
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            write (*,"(A, I12, A)", advance="no") "| ", s%world_rank, " "
            do id = 1, ND
               do idir = 1, 2
                  if ( s%neighbors(id,idir) .eq. MPI_PROC_NULL ) then
                     write (*,"(A)", advance="no") "|          "
                  else
                     write (*,"(A, I8, A)", advance="no") "| ", &
                        s%neighbors(id,idir), " "
                  end if
               end do
            end do
            write (*,"(A)") "|"
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") ""
      end if
   end subroutine print_process_neighbors

   subroutine read_general_namelist( s, filename )
      character(len=STRING_LENGTH), intent(in) :: filename
      character(len=STRING_LENGTH) :: case_name
      integer :: ierr, file_unit, nb, ng
      type(WB_State), intent(inout) :: s
      namelist /general/ case_name, nb, ng

      if ( s%world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         read( unit=file_unit, nml=general )
         close( unit=file_unit )

         s%nb = nb
         s%ng = ng

         if ( s%nb .gt. s%world_size ) then
            write (*,"(A, A)") PROGRAM_NAME, &
               ": number of blocks is greater than world size"
            call mpi_abort( MPI_COMM_WORLD, MPI_ERR_RANK, ierr )
         end if
      end if

      call mpi_bcast( case_name, STRING_LENGTH, MPI_CHARACTER, &
         WORLD_MASTER, MPI_COMM_WORLD, ierr )

      s%case_name = trim(case_name)

      call mpi_bcast( s%nb, 1, MPI_INTEGER, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
      call mpi_bcast( s%ng, 1, MPI_INTEGER, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
   end subroutine read_general_namelist

   subroutine read_block_namelists( s, filename )
      character(len=STRING_LENGTH), intent(in) :: filename
      integer :: ierr, file_unit, ib, ib_loop, id
      type(WB_State), intent(inout) :: s
      integer, dimension(ND) :: np, nx, neighbors_l, neighbors_u
      namelist /block/ ib, np, nx, neighbors_l, neighbors_u

      allocate( s%blocks(s%nb) )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         do ib_loop = 1, s%nb
            read( unit=file_unit, nml=block )

            if ( ib .gt. s%nb .or. ib .lt. 1 ) then
               write (*,"(A, A, I1, A, I1, A)") &
                  PROGRAM_NAME, ": block ", ib, &
                  " is out of range (min = 1 and max = ", s%nb, ")"
               call mpi_abort( MPI_COMM_WORLD, MPI_ERR_RANK, ierr )
            end if

            s%blocks(ib)%block_size = product(np)
            s%blocks(ib)%np = np
            s%blocks(ib)%neighbors(:,1) = neighbors_l
            s%blocks(ib)%neighbors(:,2) = neighbors_u
            s%blocks(ib)%nx = nx

            do id = 1, ND
               if ( neighbors_l(id) .eq. ib .and. neighbors_u(id) .eq. ib ) then
                  s%blocks(ib)%periods(id) = .true.
               else
                  s%blocks(ib)%periods(id) = .false.
               end if
            end do
         end do
         close( unit=file_unit )

         if ( sum( s%blocks(:)%block_size ) .ne. s%world_size ) then
            write (*,"(A, A)") PROGRAM_NAME, &
               ": block domain decomposition does not match world size"
            call mpi_abort( MPI_COMM_WORLD, MPI_ERR_RANK, ierr )
         end if
      end if

      do ib = 1, s%nb
         call mpi_bcast( s%blocks(ib)%block_size, 1, MPI_INTEGER, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%np, ND, MPI_INTEGER, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%neighbors, ND*2, MPI_INTEGER, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%nx, ND, MPI_INTEGER, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%periods, ND, MPI_LOGICAL, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
      end do
   end subroutine read_block_namelists

   subroutine setup_processes( s )
      integer :: assigned_processes, ib, id, ierr, total_points, world_rank
      type(MPI_Comm) :: comm_split
      type(WB_State), intent(inout) :: s

      allocate( s%processes(0:s%world_size-1) )
      ib = 1
      assigned_processes = 0
      do world_rank = 0, s%world_size-1
         s%processes(world_rank)%ib = ib
         assigned_processes = assigned_processes + 1
         if ( assigned_processes .eq. s%blocks(ib)%block_size ) then
            assigned_processes = 0
            ib = ib + 1
         end if
      end do

      s%ib = s%processes(s%world_rank)%ib
      call mpi_comm_split( MPI_COMM_WORLD, s%ib, 0, comm_split, ierr )
      call mpi_cart_create( comm_split, ND, s%blocks(s%ib)%np, &
         s%blocks(s%ib)%periods, s%blocks(s%ib)%reorder, s%comm_block, ierr )
      call mpi_comm_free( comm_split, ierr )
      call mpi_comm_rank( s%comm_block, s%block_rank, ierr )
      call mpi_comm_size( s%comm_block, s%block_size, ierr )
      call mpi_cart_coords( s%comm_block, s%block_rank, ND, s%block_coords, &
         ierr )

      do id = 1, ND
         s%nx(id) = s%blocks(s%ib)%nx(id) / s%blocks(s%ib)%np(id)
         if ( s%block_coords(id) .eq. s%blocks(s%ib)%np(id)-1 ) then
            s%nx(id) = s%nx(id) + modulo( s%blocks(s%ib)%nx(id), &
               s%blocks(s%ib)%np(id) )
         end if
      end do

      s%processes(s%world_rank)%block_rank = s%block_rank
      s%processes(s%world_rank)%block_coords = s%block_coords
      s%processes(s%world_rank)%nx = s%nx
      do world_rank = 0, s%world_size-1
         call mpi_bcast( s%processes(world_rank)%block_rank, 1, &
            MPI_INTEGER, world_rank, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%processes(world_rank)%block_coords, 3, &
            MPI_INTEGER, world_rank, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%processes(world_rank)%nx, 3, &
            MPI_INTEGER, world_rank, MPI_COMM_WORLD, ierr )
      end do

      ! Check if the sum of the points in a block's processes equals the total
      ! number of points.
      call mpi_reduce( product(s%nx), total_points, ND, MPI_INTEGER, MPI_SUM, &
         BLOCK_MASTER, s%comm_block, ierr )
      if ( s%block_rank .eq. BLOCK_MASTER .and. &
         product(s%blocks(s%ib)%nx) .ne. total_points ) then
         write (*,"(A, A, I2, A, I8, A, I8, A)") &
            PROGRAM_NAME, ": total points in block ", s%ib, &
            " (", product(s%blocks(s%ib)%nx), &
            ") does not match sum of points in processes (", total_points, ")"
         call mpi_abort( MPI_COMM_WORLD, MPI_ERR_SIZE, ierr )
      end if

      call identify_process_neighbors( s )
   end subroutine setup_processes
end module wbbase
