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
   use iso_fortran_env, only : compiler_options, compiler_version, &
      error_unit, real64
   use mpi_f08
   implicit none

   private

   public WB_State, check_input_file, deallocate_state, initialize_state, &
      print_initial_information, wb_abort

   integer, public, parameter ::            FP = real64
   integer, public, parameter ::         N_DIM = 3
   integer, public, parameter ::  BLOCK_MASTER = 0
   integer, public, parameter ::  WORLD_MASTER = 0
   integer, public, parameter :: STRING_LENGTH = 64

   integer, public, parameter :: NO_BLOCK_NEIGHBOR = 0

   integer, public, parameter ::     N_DIR = 2
   integer, public, parameter :: LOWER_DIR = 1
   integer, public, parameter :: UPPER_DIR = 2

   character(len=*), public, parameter :: PROGRAM_NAME = "windbag"
   character(len=*), public, parameter ::      VERSION = "0.0.0"

   type(MPI_Datatype), public, save :: MPI_FP

   type WB_Block
      private
      integer :: block_size
      integer, dimension(:), allocatable :: np, nx
      integer, dimension(:,:), allocatable :: neighbors
      logical, dimension(:), allocatable :: periods
      logical :: reorder = .false.
   end type WB_Block

   type WB_Process
      private
      integer :: block_rank
      integer, dimension(:), allocatable :: block_coords, nx
      integer :: ib
   end type WB_Process

   type WB_State
      character(len=:), allocatable, public :: case_name
      integer, private :: block_rank, block_size
      integer, public :: world_rank, world_size
      integer, public :: ib, nb
      integer, public :: nf
      integer, public :: ng
      integer, public :: nv = 5
      integer, dimension(:), allocatable, public :: nx
      integer, dimension(:), allocatable, private :: block_coords
      integer, dimension(:,:), allocatable, public :: neighbors
      real(FP), public :: t = 0.0_FP
      real(FP), dimension(:,:,:,:), allocatable, public :: f
      type(MPI_Comm), private :: comm_block
      type(WB_Block), dimension(:), allocatable, private :: blocks
      type(WB_Process), dimension(:), allocatable, private :: processes
   end type WB_State
contains
   subroutine check_block_neighbors( s )
      type(WB_State), intent(in) :: s
      integer :: ib, i_dim, i_dim_d, neighbor_l, neighbor_u

      if ( s%world_rank .eq. WORLD_MASTER ) then
         ! Ensure that lower and upper pairs exist, and that their number of
         ! points and processes match.
         do ib = 1, s%nb
            do i_dim = 1, N_DIM
               neighbor_l = s%blocks(ib)%neighbors(i_dim,LOWER_DIR)
               if ( neighbor_l .ne. NO_BLOCK_NEIGHBOR ) then
                  neighbor_u = s%blocks(neighbor_l)%neighbors(i_dim,UPPER_DIR)
                  if ( ib .ne. neighbor_u ) then
                     call wb_abort( "lower face of block N1 does not neighbor &
                                    &upper face of block N2 in direction N3", &
                                     MPI_ERR_TOPOLOGY, &
                                     (/ ib, neighbor_l, i_dim /) )
                  else
                     do i_dim_d = 1, N_DIM
                        if ( i_dim_d .ne. i_dim .and. &
                           s%blocks(ib)%np(i_dim_d) .ne. &
                           s%blocks(neighbor_l)%np(i_dim_d) ) then
                           call wb_abort( "face in direction N1 shared by &
                                          &blocks N2 and N3 does not match &
                                          &processes in direction N4", &
                              MPI_ERR_TOPOLOGY, &
                              (/ i_dim, ib, neighbor_l, i_dim_d /) )
                        end if
                        if ( i_dim_d .ne. i_dim .and. &
                           s%blocks(ib)%nx(i_dim_d) .ne. &
                           s%blocks(neighbor_l)%nx(i_dim_d) ) then
                           call wb_abort( "face in direction N1 shared by &
                                          &blocks N2 and N3 does not match &
                                          &points in direction N4", &
                              MPI_ERR_TOPOLOGY, &
                              (/ i_dim, ib, neighbor_l, i_dim_d /) )
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
            call wb_abort( "no input file given", MPI_ERR_BAD_FILE )
         end if
         call get_command_argument( 1, filename, filename_length, ierr )
         inquire( file=filename, exist=file_exists, iostat=ierr )
         if ( file_exists .eqv. .false. ) then
            call wb_abort( "input file does not exist", MPI_ERR_NO_SUCH_FILE )
         end if
      end if
      call mpi_barrier( MPI_COMM_WORLD, ierr )
   end subroutine check_input_file

   subroutine deallocate_state( s )
      integer :: ib, ierr, world_rank
      type(WB_State), intent(inout) :: s

      deallocate( s%case_name )
      deallocate( s%nx )
      deallocate( s%block_coords )
      deallocate( s%neighbors )
      call mpi_comm_free( s%comm_block, ierr )

      do ib = 1, s%nb
         deallocate( s%blocks(ib)%np )
         deallocate( s%blocks(ib)%nx )
         deallocate( s%blocks(ib)%neighbors )
         deallocate( s%blocks(ib)%periods )
      end do
      deallocate( s%blocks )

      do world_rank = 0, s%world_size-1
         deallocate( s%processes(world_rank)%block_coords )
         deallocate( s%processes(world_rank)%nx )
      end do
      deallocate( s%processes )
   end subroutine deallocate_state

   subroutine find_mpi_fp
      integer :: mpi_float_size, ierr
      call mpi_sizeof( 1.0_FP, mpi_float_size, ierr )
      call mpi_type_match_size( MPI_TYPECLASS_REAL, mpi_float_size, MPI_FP, &
         ierr )
   end subroutine find_mpi_fp

   subroutine identify_process_neighbors( s )
      integer :: i_dim, i_dir, ierr, block_neighbor, world_rank
      integer, dimension(N_DIR) :: block_ranks
      integer, dimension(:), allocatable :: block_coords
      type(WB_State), intent(inout) :: s

      allocate( block_coords(N_DIM) )
      do i_dim = 1, N_DIM
         call mpi_cart_shift( s%comm_block, i_dim-1, 1, &
            block_ranks(LOWER_DIR), block_ranks(UPPER_DIR), ierr )
         do i_dir = 1, N_DIR
            if ( block_ranks(i_dir) .ne. MPI_PROC_NULL ) then
               do world_rank = 0, s%world_size-1
                  if ( s%processes(world_rank)%ib .eq. s%ib .and. &
                     s%processes(world_rank)%block_rank .eq. &
                     block_ranks(i_dir) ) then
                     exit
                  end if
               end do
               s%neighbors(i_dim,i_dir) = world_rank
            else
               block_neighbor = s%blocks(s%ib)%neighbors(i_dim,i_dir)
               if ( block_neighbor .eq. NO_BLOCK_NEIGHBOR ) then
                  s%neighbors(i_dim,i_dir) = MPI_PROC_NULL
               else
                  ! This block neighbors another block.  This neighboring block
                  ! sits on the opposite side of the current block.  Both
                  ! processes share the same block coordinates except for the
                  ! current spatial dimension i_dim.  The block coordinates for
                  ! dimension i_dim would be opposites for each.
                  block_coords = s%block_coords
                  if ( i_dir .eq. LOWER_DIR ) then
                     block_coords(i_dim) = s%blocks(block_neighbor)%np(i_dim)-1
                  else
                     block_coords(i_dim) = 0
                  end if
                  do world_rank = 0, s%world_size-1
                     if ( s%processes(world_rank)%ib .eq. block_neighbor &
                        .and. all( s%processes(world_rank)%block_coords &
                        .eq. block_coords ) ) then
                        exit
                     end if
                  end do
                  s%neighbors(i_dim,i_dir) = world_rank
               end if
            end if
         end do
      end do
      deallocate( block_coords )
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
      integer :: ib, i_dim
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
            do i_dim = 1, N_DIM
               write (*,"(I3, A)", advance="no") s%blocks(ib)%np(i_dim), ","
            end do
            write (*,"(A, I12, A)", advance="no") ") | ", &
               product(s%blocks(ib)%nx), " | ("
            do i_dim = 1, N_DIM
               write (*,"(I4, A)", advance="no") s%blocks(ib)%nx(i_dim), ","
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
      integer :: i_dim, ierr, string_length, world_rank
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
            do i_dim = 1, N_DIM
               write (*,"(I3, A)", advance="no") &
                  s%block_coords(i_dim), ","
            end do
            write (*,"(A, I9, A)", advance="no") ") | ", &
               product(s%nx), " | ("
            do i_dim = 1, N_DIM
               write (*,"(I3, A)", advance="no") &
                  s%nx(i_dim), ","
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
      integer :: i_dim, i_dir, ierr, world_rank
      type(WB_State), intent(in) :: s

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") "## Process neighbors"
         write (*,"(A)") ""

         write (*,"(A)", advance="no") "| `world_rank` "
         do i_dim = 1, N_DIM
            write (*,"(A, I1, A)", advance="no") "|       ", i_dim, "L "
            write (*,"(A, I1, A)", advance="no") "|       ", i_dim, "U "
         end do
         write (*,"(A)") "|"

         write (*,"(A)", advance="no") "| -----------: "
         do i_dim = 1, N_DIM
            write (*,"(A, I1, A)", advance="no") "| -------: "
            write (*,"(A, I1, A)", advance="no") "| -------: "
         end do
         write (*,"(A)") "|"
      end if

      do world_rank = 0, s%world_size-1
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            write (*,"(A, I12, A)", advance="no") "| ", s%world_rank, " "
            do i_dim = 1, N_DIM
               do i_dir = 1, N_DIR
                  if ( s%neighbors(i_dim,i_dir) .eq. MPI_PROC_NULL ) then
                     write (*,"(A)", advance="no") "|          "
                  else
                     write (*,"(A, I8, A)", advance="no") "| ", &
                        s%neighbors(i_dim,i_dir), " "
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
            call wb_abort( "number of blocks is greater than world size", &
               MPI_ERR_RANK )
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
      integer :: ierr, file_unit, ib, ib_loop, i_dim
      type(WB_State), intent(inout) :: s
      integer, dimension(N_DIM) :: np, nx, neighbors_l, neighbors_u
      namelist /block/ ib, np, nx, neighbors_l, neighbors_u

      allocate( s%blocks(s%nb) )
      do ib = 1, s%nb
         allocate( s%blocks(ib)%np(N_DIM) )
         allocate( s%blocks(ib)%nx(N_DIM) )
         allocate( s%blocks(ib)%neighbors(N_DIM,N_DIR) )
         allocate( s%blocks(ib)%periods(N_DIM) )
      end do

      if ( s%world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         do ib_loop = 1, s%nb
            read( unit=file_unit, nml=block )

            if ( ib .gt. s%nb .or. ib .lt. 1 ) then
               call wb_abort( "block N1 is out of acceptable range", &
                  MPI_ERR_RANK, (/ ib /) )
            end if

            s%blocks(ib)%block_size = product(np)
            s%blocks(ib)%np = np
            s%blocks(ib)%neighbors(:,LOWER_DIR) = neighbors_l
            s%blocks(ib)%neighbors(:,UPPER_DIR) = neighbors_u
            s%blocks(ib)%nx = nx

            do i_dim = 1, N_DIM
               if ( neighbors_l(i_dim) .eq. ib .and. &
                  neighbors_u(i_dim) .eq. ib ) then
                  s%blocks(ib)%periods(i_dim) = .true.
               else
                  s%blocks(ib)%periods(i_dim) = .false.
               end if
            end do
         end do
         close( unit=file_unit )

         if ( sum( s%blocks(:)%block_size ) .ne. s%world_size ) then
            call wb_abort( &
               "block domain decomposition does not match world size", &
               MPI_ERR_RANK )
         end if
      end if

      do ib = 1, s%nb
         call mpi_bcast( s%blocks(ib)%block_size, 1, MPI_INTEGER, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%np, N_DIM, MPI_INTEGER, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%neighbors, N_DIM*N_DIR, MPI_INTEGER, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%nx, N_DIM, MPI_INTEGER, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%periods, N_DIM, MPI_LOGICAL, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
      end do
   end subroutine read_block_namelists

   subroutine setup_processes( s )
      integer :: assigned_processes, ib, i_dim, ierr, total_points, world_rank
      type(MPI_Comm) :: comm_split
      type(WB_State), intent(inout) :: s

      allocate( s%nx(N_DIM) )
      allocate( s%block_coords(N_DIM) )
      allocate( s%neighbors(N_DIM,N_DIR) )

      allocate( s%processes(0:s%world_size-1) )
      do world_rank = 0, s%world_size-1
         allocate( s%processes(world_rank)%block_coords(N_DIM) )
         allocate( s%processes(world_rank)%nx(N_DIM) )
      end do

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
      call mpi_cart_create( comm_split, N_DIM, s%blocks(s%ib)%np, &
         s%blocks(s%ib)%periods, s%blocks(s%ib)%reorder, s%comm_block, ierr )
      call mpi_comm_free( comm_split, ierr )
      call mpi_comm_rank( s%comm_block, s%block_rank, ierr )
      call mpi_comm_size( s%comm_block, s%block_size, ierr )
      call mpi_cart_coords( s%comm_block, s%block_rank, N_DIM, &
         s%block_coords, ierr )

      do i_dim = 1, N_DIM
         s%nx(i_dim) = s%blocks(s%ib)%nx(i_dim) / s%blocks(s%ib)%np(i_dim)
         if ( s%block_coords(i_dim) .eq. s%blocks(s%ib)%np(i_dim)-1 ) then
            s%nx(i_dim) = s%nx(i_dim) + modulo( s%blocks(s%ib)%nx(i_dim), &
               s%blocks(s%ib)%np(i_dim) )
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
      call mpi_reduce( product(s%nx), total_points, N_DIM, MPI_INTEGER, &
         MPI_SUM, BLOCK_MASTER, s%comm_block, ierr )
      if ( s%block_rank .eq. BLOCK_MASTER .and. &
         product(s%blocks(s%ib)%nx) .ne. total_points ) then
         call wb_abort( "total points in block N1 does not match sum of &
                        &points in individual processes", &
            MPI_ERR_SIZE, &
            (/ s%ib /) )
      end if

      call identify_process_neighbors( s )
   end subroutine setup_processes

   subroutine wb_abort( message, errno, ints, floats )
      character(len=*), intent(in) :: message
      integer, intent(in) :: errno
      integer :: i, ierr
      integer, dimension(:), optional, intent(in) :: ints
      real(FP), dimension(:), optional, intent(in) :: floats

      write (error_unit, "(A, A, A)") PROGRAM_NAME, ": ", message
      if ( present(ints) ) then
         do i = 1, size(ints)
            write (error_unit, "(A, I1, A, I8)") "N", i, " = ", ints(i)
         end do
      end if
      if ( present(floats) ) then
         do i = 1, size(floats)
            write (error_unit, "(A, I1, A, ES9.2)") "F", i, " = ", floats(i)
         end do
      end if
      call mpi_abort( MPI_COMM_WORLD, errno, ierr )
   end subroutine wb_abort
end module wbbase
