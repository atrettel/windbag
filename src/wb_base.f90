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
module wb_base
   use iso_fortran_env
   use mpi_f08
   use wb_text
   implicit none

   private

   public WB_State, check_input_file, deallocate_state, initialize_state, &
      print_initial_information, wb_abort

   integer, public, parameter       ::     FP = real64
   type(MPI_Datatype), public, save :: MPI_FP

   integer, public, parameter       ::     SP = int64
   type(MPI_Datatype), public, save :: MPI_SP

   integer, public, parameter :: BLOCK_MASTER = 0
   integer, public, parameter :: WORLD_MASTER = 0

   integer, public, parameter :: NO_BLOCK_NEIGHBOR = 0

   integer, public, parameter ::     N_DIR = 2
   integer, public, parameter :: LOWER_DIR = 1
   integer, public, parameter :: UPPER_DIR = 2

   integer, public, parameter :: DEFAULT_BLOCK_NEIGHBOR = 1
   integer, public, parameter :: DEFAULT_IB             = 1
   integer, public, parameter :: DEFAULT_N_DIM          = 3
   integer, public, parameter :: DEFAULT_NB             = 1
   integer, public, parameter :: DEFAULT_NG             = 3
   integer, public, parameter :: DEFAULT_NP             = 0
   integer, public, parameter :: DEFAULT_NX             = 0

   integer, public, parameter :: STRING_LENGTH = 64

   character(len=*), public, parameter ::      PROGRAM_NAME = "windbag"
   character(len=*), public, parameter ::           VERSION = "0.0.0"
   character(len=*), public, parameter :: DEFAULT_CASE_NAME = "casename"

   ! Big-endian architectures put the most significant byte first, and
   ! little-endian architectures put the least significant byte first.  Binary
   ! output requires knowing the architecture's endianness.  This statement
   ! casts an integer for 1 as a single character (the X is arbitrary) and then
   ! gets the ASCII code for the first character (the result only contains the
   ! leading bits of the integer).  If the result is 1, then the architecture
   ! is little-endian, since the least significant byte came first.  If the
   ! result is 0, then the architecture is big-endian.  Little-endian
   ! architectures are more common nowadays.
   logical, public, parameter :: ARCH_IS_BIG_ENDIAN = &
      ichar( transfer( 1_SP, "X" ) ) .eq. 0

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
      integer, public :: n_dim
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
      integer :: ib, i_dim, j_dim, neighbor_l, neighbor_u

      if ( s%world_rank .eq. WORLD_MASTER ) then
         ! Ensure that lower and upper pairs exist, and that their number of
         ! points and processes match.
         do ib = 1, s%nb
            do i_dim = 1, s%n_dim
               neighbor_l = s%blocks(ib)%neighbors(i_dim,LOWER_DIR)
               if ( neighbor_l .ne. NO_BLOCK_NEIGHBOR ) then
                  neighbor_u = s%blocks(neighbor_l)%neighbors(i_dim,UPPER_DIR)
                  if ( ib .ne. neighbor_u ) then
                     call wb_abort( "lower face of block N1 does not neighbor &
                                    &upper face of block N2 in direction N3", &
                                     MPI_ERR_TOPOLOGY, &
                                     int( (/ ib, neighbor_l, i_dim /), SP ) )
                  else
                     do j_dim = 1, s%n_dim
                        if ( j_dim .ne. i_dim .and. &
                           s%blocks(ib)%np(j_dim) .ne. &
                           s%blocks(neighbor_l)%np(j_dim) ) then
                           call wb_abort( "face in direction N1 shared by &
                                          &blocks N2 and N3 does not match &
                                          &processes in direction N4", &
                              MPI_ERR_TOPOLOGY, &
                              int( (/ i_dim, ib, neighbor_l, j_dim /), SP ) )
                        end if
                        if ( j_dim .ne. i_dim .and. &
                           s%blocks(ib)%nx(j_dim) .ne. &
                           s%blocks(neighbor_l)%nx(j_dim) ) then
                           call wb_abort( "face in direction N1 shared by &
                                          &blocks N2 and N3 does not match &
                                          &points in direction N4", &
                              MPI_ERR_TOPOLOGY, &
                              int( (/ i_dim, ib, neighbor_l, j_dim /), SP ) )
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
      call find_mpi_precisions
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

   subroutine find_mpi_precisions
      integer :: mpi_float_size, mpi_int_size, ierr

      call mpi_sizeof( 1.0_FP, mpi_float_size, ierr )
      call mpi_type_match_size( MPI_TYPECLASS_REAL, mpi_float_size, MPI_FP, &
         ierr )

      call mpi_sizeof( 1_SP, mpi_int_size, ierr )
      call mpi_type_match_size( MPI_TYPECLASS_INTEGER, mpi_int_size, MPI_SP, &
         ierr )
   end subroutine find_mpi_precisions

   subroutine identify_process_neighbors( s )
      integer :: i_dim, i_dir, ierr, block_neighbor, world_rank
      integer, dimension(N_DIR) :: block_ranks
      integer, dimension(:), allocatable :: block_coords
      type(WB_State), intent(inout) :: s

      allocate( block_coords(s%n_dim) )
      do i_dim = 1, s%n_dim
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

   subroutine print_initial_information( s )
      type(WB_State), intent(in) :: s

      call write_global_information( output_unit, s )
      call write_block_information( output_unit, s )
      call write_process_information( output_unit, s )
      call write_process_neighbors( output_unit, s )
   end subroutine print_initial_information

   subroutine read_general_namelist( s, filename )
      character(len=STRING_LENGTH), intent(in) :: filename
      character(len=STRING_LENGTH) :: case_name=DEFAULT_CASE_NAME
      integer :: ierr, file_unit, n_dim=DEFAULT_N_DIM, nb=DEFAULT_NB, &
         ng=DEFAULT_NG
      type(WB_State), intent(inout) :: s
      namelist /general/ case_name, nb, ng, n_dim

      if ( s%world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         read( unit=file_unit, nml=general )
         close( unit=file_unit )

         s%nb    = nb
         s%ng    = ng
         s%n_dim = n_dim

         if ( s%nb .lt. 1 ) then
            call wb_abort( "number of blocks is less than 1", &
               MPI_ERR_COUNT )
         end if
         if ( s%nb .gt. s%world_size ) then
            call wb_abort( "number of blocks is greater than world size", &
               MPI_ERR_COUNT )
         end if
         if ( s%ng .lt. 1 ) then
            call wb_abort( "number of ghost points is less than 1", &
               MPI_ERR_COUNT )
         end if
         if ( s%n_dim .lt. 1 .or. s%n_dim .gt. 3 ) then
            call wb_abort( "number of dimensions must be 1, 2, or 3", &
               MPI_ERR_COUNT )
         end if
      end if

      call mpi_bcast( case_name, STRING_LENGTH, MPI_CHARACTER, &
         WORLD_MASTER, MPI_COMM_WORLD, ierr )

      s%case_name = trim(case_name)

      call mpi_bcast( s%nb, 1, MPI_INTEGER, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
      call mpi_bcast( s%ng, 1, MPI_INTEGER, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
      call mpi_bcast( s%n_dim, 1, MPI_INTEGER, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
   end subroutine read_general_namelist

   subroutine read_block_namelists( s, filename )
      character(len=STRING_LENGTH), intent(in) :: filename
      integer :: ierr, file_unit, ib=DEFAULT_IB, ib_loop, i_dim, i_dir
      type(WB_State), intent(inout) :: s
      integer, dimension(:), allocatable :: np, nx, neighbors_l, neighbors_u
      namelist /block/ ib, np, nx, neighbors_l, neighbors_u

      allocate( np(s%n_dim) )
      allocate( nx(s%n_dim) )
      allocate( neighbors_l(s%n_dim) )
      allocate( neighbors_u(s%n_dim) )

      do i_dim = 1, s%n_dim
         np(i_dim)          = DEFAULT_NP
         nx(i_dim)          = DEFAULT_NX
         neighbors_l(i_dim) = DEFAULT_BLOCK_NEIGHBOR
         neighbors_u(i_dim) = DEFAULT_BLOCK_NEIGHBOR
      end do

      allocate( s%blocks(s%nb) )
      do ib_loop = 1, s%nb
         allocate( s%blocks(ib_loop)%np(s%n_dim) )
         allocate( s%blocks(ib_loop)%nx(s%n_dim) )
         allocate( s%blocks(ib_loop)%neighbors(s%n_dim,N_DIR) )
         allocate( s%blocks(ib_loop)%periods(s%n_dim) )
      end do

      if ( s%world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         do ib_loop = 1, s%nb
            read( unit=file_unit, nml=block )

            if ( ib .gt. s%nb .or. ib .lt. 1 ) then
               call wb_abort( "block N1 is out of acceptable range", &
                  MPI_ERR_RANK, int( (/ ib /), SP ) )
            end if

            s%blocks(ib)%block_size = product(np)
            s%blocks(ib)%np = np
            s%blocks(ib)%neighbors(:,LOWER_DIR) = neighbors_l
            s%blocks(ib)%neighbors(:,UPPER_DIR) = neighbors_u
            s%blocks(ib)%nx = nx

            do i_dim = 1, s%n_dim
               if ( s%blocks(ib)%np(i_dim) .lt. 1 ) then
                  call wb_abort( "number of processes in direction N1 of &
                                 &block N2 is less than 1", &
                                 MPI_ERR_COUNT, int( (/ i_dim, ib /), SP ) )
               end if
               if ( s%blocks(ib)%nx(i_dim) .lt. s%ng ) then
                  call wb_abort( "number of points in direction N1 of block &
                                 &N2 is less than number of ghost points N3", &
                                 MPI_ERR_COUNT, &
                                 int( (/ i_dim, ib, s%ng /), SP ) )
               end if
               do i_dir = 1, N_DIR
                  if ( s%blocks(ib)%neighbors(i_dim,i_dir) .lt. &
                     NO_BLOCK_NEIGHBOR ) then
                     call wb_abort( "neighbor to block N1 in direction N2 and &
                                    &dimension N3 is negative", &
                                    MPI_ERR_COUNT, &
                                    int( (/ ib, i_dir, i_dim /), SP ) )
                  end if
               end do
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
         call mpi_bcast( s%blocks(ib)%np, s%n_dim, MPI_INTEGER, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%neighbors, s%n_dim*N_DIR, MPI_INTEGER, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%nx, s%n_dim, MPI_INTEGER, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%periods, s%n_dim, MPI_LOGICAL, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
      end do

      deallocate( np )
      deallocate( nx )
      deallocate( neighbors_l )
      deallocate( neighbors_u )
   end subroutine read_block_namelists

   subroutine setup_processes( s )
      integer :: assigned_processes, ib, i_dim, ierr, total_points, world_rank
      type(MPI_Comm) :: comm_split
      type(WB_State), intent(inout) :: s

      allocate( s%nx(s%n_dim) )
      allocate( s%block_coords(s%n_dim) )
      allocate( s%neighbors(s%n_dim,N_DIR) )

      allocate( s%processes(0:s%world_size-1) )
      do world_rank = 0, s%world_size-1
         allocate( s%processes(world_rank)%block_coords(s%n_dim) )
         allocate( s%processes(world_rank)%nx(s%n_dim) )
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
      call mpi_cart_create( comm_split, s%n_dim, s%blocks(s%ib)%np, &
         s%blocks(s%ib)%periods, s%blocks(s%ib)%reorder, s%comm_block, ierr )
      call mpi_comm_free( comm_split, ierr )
      call mpi_comm_rank( s%comm_block, s%block_rank, ierr )
      call mpi_comm_size( s%comm_block, s%block_size, ierr )
      call mpi_cart_coords( s%comm_block, s%block_rank, s%n_dim, &
         s%block_coords, ierr )

      do i_dim = 1, s%n_dim
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
      call mpi_reduce( product(s%nx), total_points, s%n_dim, MPI_INTEGER, &
         MPI_SUM, BLOCK_MASTER, s%comm_block, ierr )
      if ( s%block_rank .eq. BLOCK_MASTER .and. &
         product(s%blocks(s%ib)%nx) .ne. total_points ) then
         call wb_abort( "total points in block N1 does not match sum of &
                        &points in individual processes", &
            MPI_ERR_SIZE, int( (/ s%ib /), SP ) )
      end if

      call identify_process_neighbors( s )
   end subroutine setup_processes

   subroutine wb_abort( message, errno, ints, floats )
      character(len=*), intent(in) :: message
      integer, intent(in) :: errno
      integer :: i, ierr
      integer(SP), dimension(:), optional, intent(in) :: ints
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

   subroutine write_block_information( f, s )
      integer, intent(in) :: f
      type(WB_State), intent(in) :: s
      integer :: ib, i_dim
      character(len=STRING_LENGTH) :: label

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (f,"(A)") "## Block information"
         write (f,"(A)") ""

         call write_string_table_entry( f, "`ib`", IB_COLUMN_WIDTH )
         call write_string_table_entry( f, "`block_size`", SIZE_COLUMN_WIDTH )
         do i_dim = 1, s%n_dim
            write (label,"(A, I1, A)") "`np(", i_dim, ")`"
            call write_string_table_entry( f, trim(label), NP_COLUMN_WIDTH )
         end do
         call write_string_table_entry( f, "points", POINTS_COLUMN_WIDTH )
         do i_dim = 1, s%n_dim
            write (label,"(A, I1, A)") "`nx(", i_dim, ")`"
            call write_string_table_entry( f, trim(label), NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. s%n_dim) )
         end do

         call write_rule_table_entry( f, IB_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_rule_table_entry( f, SIZE_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1, s%n_dim
            call write_rule_table_entry( f, NP_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED )
         end do
         call write_rule_table_entry( f, POINTS_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1, s%n_dim
            call write_rule_table_entry( f, NX_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED, end_row=(i_dim .eq. s%n_dim) )
         end do

         do ib = 1, s%nb
            call write_integer_table_entry( f, int(ib,int64), IB_COLUMN_WIDTH )
            call write_integer_table_entry( f, &
               int(s%blocks(ib)%block_size,int64), SIZE_COLUMN_WIDTH )
            do i_dim = 1, s%n_dim
               call write_integer_table_entry( f, &
               int(s%blocks(ib)%np(i_dim),int64), NP_COLUMN_WIDTH )
            end do
            call write_integer_table_entry( f, &
               int(product(s%blocks(ib)%nx),int64), POINTS_COLUMN_WIDTH )
            do i_dim = 1, s%n_dim
               call write_integer_table_entry( f, &
               int(s%blocks(ib)%nx(i_dim),int64), NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. s%n_dim) )
            end do
         end do

         write (f,"(A)") ""
      end if
   end subroutine write_block_information

   subroutine write_global_information( f, s )
      integer, intent(in) :: f
      type(WB_State), intent(in) :: s
      integer :: ierr, mpi_major_version_number, mpi_minor_version_number, &
         version_length
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: lib_version

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (f,"(A, A, A, A, A, A, A)") "# ", PROGRAM_NAME, " ", VERSION, &
            ", case `", s%case_name, "`"
         write (f,"(A)") ""

         if ( ARCH_IS_BIG_ENDIAN ) then
            write (f,"(A)") "- Architecture is big-endian."
         else
            write (f,"(A)") "- Architecture is little-endian."
         end if
         write (f,"(A)") ""

         if ( FP .eq. real64 ) then
            write (f,"(A)") "- Floating point numbers are double precision."
         else if (FP .eq. real32 ) then
            write (f,"(A)") "- Floating point numbers are single precision."
         else
            write (f,"(A)") "- Floating point numbers are an unknown &
                             &precision."
         end if
         write (f,"(A)") ""

         if ( MPI_FP .eq. MPI_REAL8 ) then
            write (f,"(A)") "- MPI floating point precision is `MPI_REAL8`."
         else if ( MPI_FP .eq. MPI_REAL4 ) then
            write (f,"(A)") "- MPI floating point precision is `MPI_REAL4`."
         else if ( MPI_FP .eq. MPI_REAL ) then
            write (f,"(A)") "- MPI floating point precision is `MPI_REAL`."
         else
            write (f,"(A)") "- MPI floating point precision is unknown."
         end if
         write (f,"(A)") ""

         if ( SP .eq. int64 ) then
            write (f,"(A)") "- Signed integers are 64-bit precision."
         else if (SP .eq. int32 ) then
            write (f,"(A)") "- Signed integers are 32-bit precision."
         else
            write (f,"(A)") "- Signed integers are of unknown precision."
         end if
         write (f,"(A)") ""

         if ( MPI_SP .eq. MPI_INTEGER8 ) then
            write (f,"(A)") "- MPI integer precision is `MPI_INTEGER8`."
         else if ( MPI_SP .eq. MPI_INTEGER4 ) then
            write (f,"(A)") "- MPI integer precision is `MPI_INTEGER4`."
         else if ( MPI_SP .eq. MPI_INTEGER ) then
            write (f,"(A)") "- MPI integer precision is `MPI_INTEGER`."
         else
            write (f,"(A)") "- MPI integer precision is unknown."
         end if
         write (f,"(A)") ""

         call mpi_get_version( mpi_major_version_number, &
            mpi_minor_version_number, ierr )
         write (f,"(A, I1, A, I1)") "- MPI version: ", &
            mpi_major_version_number, ".", mpi_minor_version_number
         write (f,"(A)") ""

         call mpi_get_library_version( lib_version, version_length, ierr )
         write (f,"(A, A)") "- MPI library version: ", trim(lib_version)
         write (f,"(A)") ""

         write (f,"(A, A, A, A, A)") "- Compiled using ", compiler_version(), &
            " using the following options: `", compiler_options(), "`"
         write (f,"(A)") ""

         if ( MPI_SUBARRAYS_SUPPORTED ) then
            write (f,"(A)") "- MPI subarrays are supported."
         else
            write (f,"(A)") "- MPI subarrays are not supported."
         end if
         write (f,"(A)") ""
      end if
   end subroutine write_global_information

   subroutine write_process_information( f, s )
      integer, intent(in) :: f
      type(WB_State), intent(in) :: s
      integer :: i_dim, ierr, processor_length, world_rank
      character(len=STRING_LENGTH) :: label
      character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

      call mpi_get_processor_name( processor_name, processor_length, ierr )

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") "## Process information"
         write (*,"(A)") ""

         call write_string_table_entry( f, "`world_rank`", RANK_COLUMN_WIDTH )
         call write_string_table_entry( f, "hostname", HOSTNAME_COLUMN_WIDTH )
         call write_string_table_entry( f, "`ib`",           IB_COLUMN_WIDTH )
         call write_string_table_entry( f, "`block_rank`", RANK_COLUMN_WIDTH )
         do i_dim = 1, s%n_dim
            write (label,"(A, I1, A)") "(", i_dim, ")"
            call write_string_table_entry( f, trim(label), &
               COORDS_COLUMN_WIDTH )
         end do
         call write_string_table_entry( f, "points", POINTS_COLUMN_WIDTH )
         do i_dim = 1, s%n_dim
            write (label,"(A, I1, A)") "`nx(", i_dim, ")`"
            call write_string_table_entry( f, trim(label), NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. s%n_dim) )
         end do

         call write_rule_table_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_rule_table_entry( f, HOSTNAME_COLUMN_WIDTH, &
            alignment=LEFT_ALIGNED  )
         call write_rule_table_entry( f, IB_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_rule_table_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1, s%n_dim
            call write_rule_table_entry( f, COORDS_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED )
         end do
         call write_rule_table_entry( f, POINTS_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1, s%n_dim
            call write_rule_table_entry( f, NX_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED, end_row=(i_dim .eq. s%n_dim) )
         end do
      end if

      do world_rank = 0, s%world_size-1
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            call write_integer_table_entry( f, int(s%world_rank,int64), &
               RANK_COLUMN_WIDTH )
            call write_string_table_entry(  f, trim(processor_name), &
               HOSTNAME_COLUMN_WIDTH )
            call write_integer_table_entry( f, int(s%ib,int64), &
               IB_COLUMN_WIDTH )
            call write_integer_table_entry( f, int(s%block_rank,int64), &
               RANK_COLUMN_WIDTH )
            do i_dim = 1, s%n_dim
               call write_integer_table_entry( f, &
               int(s%block_coords(i_dim),int64), COORDS_COLUMN_WIDTH )
            end do
            call write_integer_table_entry( f, int(product(s%nx),int64), &
               POINTS_COLUMN_WIDTH )
            do i_dim = 1, s%n_dim
               call write_integer_table_entry( f, int(s%nx(i_dim),int64), &
                  NX_COLUMN_WIDTH, end_row=(i_dim .eq. s%n_dim) )
            end do
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") ""
      end if
   end subroutine write_process_information

   subroutine write_process_neighbors( f, s )
      integer, intent(in) :: f
      type(WB_State), intent(in) :: s
      integer :: i_dim, i_dir, ierr, j_dir, world_rank, face_count=0, neighbor
      integer, dimension(N_DIR) :: dirs = (/ LOWER_DIR, UPPER_DIR /)
      character(len=STRING_LENGTH) :: label

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") "## Process neighbors"
         write (*,"(A)") ""

         call write_string_table_entry( f, "`world_rank`", RANK_COLUMN_WIDTH )
         do i_dim = 1, s%n_dim
            do i_dir = 1, N_DIR
               j_dir = dirs(i_dir)
               face_count = face_count + 1
               if ( j_dir .eq. LOWER_DIR ) then
                  write (label,"(I1, A)") i_dim, "L"
               else
                  write (label,"(I1, A)") i_dim, "U"
               end if
               call write_string_table_entry( f, trim(label), &
                  RANK_COLUMN_WIDTH, &
                  end_row=( face_count .eq. s%n_dim*N_DIR ) )
            end do
         end do

         face_count = 0
         call write_rule_table_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1, s%n_dim
            do i_dir = 1, N_DIR
               face_count = face_count + 1
               call write_rule_table_entry( f, RANK_COLUMN_WIDTH, &
                  alignment=RIGHT_ALIGNED, &
                  end_row=( face_count .eq. s%n_dim*N_DIR ) )
            end do
         end do
      end if

      do world_rank = 0, s%world_size-1
         face_count = 0
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            call write_integer_table_entry( f, int(s%world_rank,int64), &
               RANK_COLUMN_WIDTH )
            do i_dim = 1, s%n_dim
               do i_dir = 1, N_DIR
                  j_dir = dirs(i_dir)
                  face_count = face_count + 1
                  neighbor = s%neighbors(i_dim,j_dir)
                  if ( neighbor .eq. MPI_PROC_NULL ) then
                     call write_string_table_entry( f, "-", &
                        RANK_COLUMN_WIDTH, &
                        end_row=( face_count .eq. s%n_dim*N_DIR ) )
                  else
                     call write_integer_table_entry( f, int(neighbor,int64), &
                        RANK_COLUMN_WIDTH, &
                        end_row=( face_count .eq. s%n_dim*N_DIR ) )
                  end if
               end do
            end do
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") ""
      end if
   end subroutine write_process_neighbors
end module wb_base
