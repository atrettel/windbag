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
   use wb_representation
   use wb_text
   implicit none

   private

   public WB_State, check_input_file, deallocate_state, initialize_state, &
      print_initial_information, wb_abort

   integer(MP), public, parameter :: BLOCK_MASTER = 0_MP
   integer(MP), public, parameter :: WORLD_MASTER = 0_MP

   integer(SP), public, parameter :: NO_BLOCK_NEIGHBOR = 0_SP

   integer(SP), public, parameter ::     N_DIR = 2_SP
   integer(SP), public, parameter :: LOWER_DIR = 1_SP
   integer(SP), public, parameter :: UPPER_DIR = 2_SP

   integer(SP), public, parameter :: DEFAULT_BLOCK_NEIGHBOR = 1_SP
   integer(SP), public, parameter :: DEFAULT_IB             = 1_SP
   integer(SP), public, parameter :: DEFAULT_N_DIM          = 3_SP
   integer(SP), public, parameter :: DEFAULT_NB             = 1_SP
   integer(SP), public, parameter :: DEFAULT_NG             = 3_SP
   integer(SP), public, parameter :: DEFAULT_NX             = 0_SP

   integer(MP), public, parameter :: DEFAULT_NP             = 0_MP

   character(len=*), public, parameter ::      PROGRAM_NAME = "windbag"
   character(len=*), public, parameter ::           VERSION = "0.0.0"
   character(len=*), public, parameter :: DEFAULT_CASE_NAME = "casename"

   type WB_Block
      private
      integer(MP) :: block_size
      integer(MP), dimension(:), allocatable :: np
      integer(SP), dimension(:), allocatable :: nx
      integer(SP), dimension(:,:), allocatable :: neighbors
      logical, dimension(:), allocatable :: periods
      logical :: reorder = .false.
   end type WB_Block

   type WB_Process
      private
      integer(MP) :: block_rank
      integer(MP), dimension(:), allocatable :: block_coords
      integer(SP), dimension(:), allocatable :: nx
      integer(SP) :: ib
   end type WB_Process

   type WB_State
      character(len=:), allocatable, public :: case_name
      integer(MP), private :: block_rank, block_size
      integer(MP), public :: world_rank, world_size
      integer(SP), public :: ib, nb
      integer(SP), public :: n_dim
      integer(SP), public :: nf
      integer(SP), public :: ng
      integer(SP), public :: nv = 5_SP
      integer(SP), dimension(:), allocatable, public :: nx
      integer(MP), dimension(:), allocatable, private :: block_coords
      integer(MP), dimension(:,:), allocatable, public :: neighbors
      real(FP), public :: t = 0.0_FP
      real(FP), dimension(:,:,:,:), allocatable, public :: f
      type(MPI_Comm), private :: comm_block
      type(WB_Block), dimension(:), allocatable, private :: blocks
      type(WB_Process), dimension(:), allocatable, private :: processes
   end type WB_State
contains
   subroutine check_block_neighbors( s )
      type(WB_State), intent(in) :: s
      integer(SP) :: ib, i_dim, j_dim, neighbor_l, neighbor_u

      ! Ensure that lower and upper pairs exist, and that their number of
      ! points and processes match.  Since this is just checking and not
      ! calculation, it must occur on the world master.  It is impossible to
      ! have each block master check this for their own blocks since the block
      ! communicators do not exist yet.  If that were possible, it would
      ! eliminate the loop over all blocks.
      if ( s%world_rank .eq. WORLD_MASTER ) then
         do ib = 1_SP, s%nb
            do i_dim = 1_SP, s%n_dim
               neighbor_l = s%blocks(ib)%neighbors(i_dim,LOWER_DIR)
               if ( neighbor_l .ne. NO_BLOCK_NEIGHBOR ) then
                  neighbor_u = s%blocks(neighbor_l)%neighbors(i_dim,UPPER_DIR)
                  if ( ib .ne. neighbor_u ) then
                     call wb_abort( "lower face of block N1 does not neighbor &
                                    &upper face of block N2 in direction N3", &
                        MPI_ERR_TOPOLOGY, &
                        (/ ib, neighbor_l, i_dim /) )
                  else
                     do j_dim = 1_SP, s%n_dim
                        if ( j_dim .ne. i_dim .and. &
                           s%blocks(ib)%np(j_dim) .ne. &
                           s%blocks(neighbor_l)%np(j_dim) ) then
                           call wb_abort( "face in direction N1 shared by &
                                          &blocks N2 and N3 does not match &
                                          &processes in direction N4", &
                              MPI_ERR_TOPOLOGY, &
                              (/ i_dim, ib, neighbor_l, j_dim /) )
                        end if
                        if ( j_dim .ne. i_dim .and. &
                           s%blocks(ib)%nx(j_dim) .ne. &
                           s%blocks(neighbor_l)%nx(j_dim) ) then
                           call wb_abort( "face in direction N1 shared by &
                                          &blocks N2 and N3 does not match &
                                          &points in direction N4", &
                              MPI_ERR_TOPOLOGY, &
                              (/ i_dim, ib, neighbor_l, j_dim /) )
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
      integer :: argc, filename_length, ierr
      integer(MP) :: ierr_mpi, world_rank
      logical :: file_exists

      call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr_mpi )
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
      integer(SP) :: ib
      integer(MP) :: ierr, world_rank
      type(WB_State), intent(inout) :: s

      deallocate( s%case_name )
      deallocate( s%nx )
      deallocate( s%block_coords )
      deallocate( s%neighbors )
      call mpi_comm_free( s%comm_block, ierr )

      do ib = 1_SP, s%nb
         deallocate( s%blocks(ib)%np )
         deallocate( s%blocks(ib)%nx )
         deallocate( s%blocks(ib)%neighbors )
         deallocate( s%blocks(ib)%periods )
      end do
      deallocate( s%blocks )

      do world_rank = 0_MP, s%world_size-1_MP
         deallocate( s%processes(world_rank)%block_coords )
         deallocate( s%processes(world_rank)%nx )
      end do
      deallocate( s%processes )
   end subroutine deallocate_state

   subroutine identify_process_neighbors( s )
      integer(MP) :: ierr, world_rank
      integer(SP) :: block_neighbor, i_dir, i_dim
      integer(MP), dimension(N_DIR) :: block_ranks
      integer(MP), dimension(:), allocatable :: block_coords
      type(WB_State), intent(inout) :: s

      allocate( block_coords(s%n_dim) )
      do i_dim = 1_SP, s%n_dim
         call mpi_cart_shift( s%comm_block, int(i_dim,MP)-1_MP, 1_MP, &
            block_ranks(LOWER_DIR), block_ranks(UPPER_DIR), ierr )
         do i_dir = 1_SP, N_DIR
            if ( block_ranks(i_dir) .ne. MPI_PROC_NULL ) then
               do world_rank = 0_MP, s%world_size-1_MP
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
                     block_coords(i_dim) = &
                        s%blocks(block_neighbor)%np(i_dim)-1_MP
                  else
                     block_coords(i_dim) = 0_MP
                  end if
                  do world_rank = 0_MP, s%world_size-1_MP
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
      integer(MP) :: ierr
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

      call write_global_information(  output_unit, s )
      call write_block_information(   output_unit, s )
      call write_process_information( output_unit, s )
      call write_process_neighbors(   output_unit, s )
   end subroutine print_initial_information

   subroutine read_general_namelist( s, filename )
      character(len=STRING_LENGTH), intent(in) :: filename
      character(len=STRING_LENGTH) :: case_name=DEFAULT_CASE_NAME
      integer :: file_unit
      integer(MP) :: ierr
      integer(SP) :: nb=DEFAULT_NB, n_dim=DEFAULT_N_DIM, ng=DEFAULT_NG
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

         if ( s%nb .lt. 1_SP ) then
            call wb_abort( "number of blocks is less than 1", &
               MPI_ERR_COUNT )
         end if
         if ( s%nb .gt. int(s%world_size,SP) ) then
            ! This is a feature of the code and not a bug.  It allows the code
            ! to treat communication between blocks as the same as
            ! communication between processes, but that plan only works if each
            ! block uses at least one process.
            call wb_abort( "number of blocks is greater than world size", &
               MPI_ERR_COUNT )
         end if
         if ( s%ng .lt. 1_SP ) then
            call wb_abort( "number of ghost points is less than 1", &
               MPI_ERR_COUNT )
         end if
         if ( s%n_dim .lt. 1_SP .or. s%n_dim .gt. 3_SP ) then
            call wb_abort( "number of dimensions must be 1, 2, or 3", &
               MPI_ERR_COUNT )
         end if
      end if

      call mpi_bcast( case_name, int(STRING_LENGTH,MP), MPI_CHARACTER, &
         WORLD_MASTER, MPI_COMM_WORLD, ierr )

      s%case_name = trim(case_name)

      call mpi_bcast( s%nb, 1_MP, MPI_SP, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
      call mpi_bcast( s%ng, 1_MP, MPI_SP, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
      call mpi_bcast( s%n_dim, 1_MP, MPI_SP, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
   end subroutine read_general_namelist

   subroutine read_block_namelists( s, filename )
      character(len=STRING_LENGTH), intent(in) :: filename
      integer :: file_unit
      integer(MP) :: ierr
      integer(SP) :: ib=DEFAULT_IB, jb, i_dim, i_dir
      type(WB_State), intent(inout) :: s
      integer(MP), dimension(:), allocatable :: np
      integer(SP), dimension(:), allocatable :: nx, neighbors_l, neighbors_u
      namelist /block/ ib, np, nx, neighbors_l, neighbors_u

      allocate( np(s%n_dim) )
      allocate( nx(s%n_dim) )
      allocate( neighbors_l(s%n_dim) )
      allocate( neighbors_u(s%n_dim) )

      np(:)          = DEFAULT_NP
      nx(:)          = DEFAULT_NX
      neighbors_l(:) = DEFAULT_BLOCK_NEIGHBOR
      neighbors_u(:) = DEFAULT_BLOCK_NEIGHBOR

      allocate( s%blocks(s%nb) )
      do jb = 1_SP, s%nb
         allocate( s%blocks(jb)%np(s%n_dim) )
         allocate( s%blocks(jb)%nx(s%n_dim) )
         allocate( s%blocks(jb)%neighbors(s%n_dim,N_DIR) )
         allocate( s%blocks(jb)%periods(s%n_dim) )
      end do

      if ( s%world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         do jb = 1_SP, s%nb
            read( unit=file_unit, nml=block )

            if ( ib .gt. s%nb .or. ib .lt. 1_SP ) then
               call wb_abort( "block N1 is out of acceptable range", &
                  MPI_ERR_RANK, (/ ib /) )
            end if

            s%blocks(ib)%block_size = product(np)
            s%blocks(ib)%np = np
            s%blocks(ib)%neighbors(:,LOWER_DIR) = neighbors_l
            s%blocks(ib)%neighbors(:,UPPER_DIR) = neighbors_u
            s%blocks(ib)%nx = nx

            do i_dim = 1_SP, s%n_dim
               if ( s%blocks(ib)%np(i_dim) .lt. 1_MP ) then
                  call wb_abort( "number of processes in direction N1 of &
                                 &block N2 is less than 1", &
                                 MPI_ERR_COUNT, (/ i_dim, ib /) )
               end if
               if ( s%blocks(ib)%nx(i_dim) .lt. s%ng ) then
                  call wb_abort( "number of points in direction N1 of block &
                                 &N2 is less than number of ghost points N3", &
                                 MPI_ERR_COUNT, &
                                 (/ i_dim, ib, s%ng /) )
               end if
               do i_dir = 1_SP, N_DIR
                  if ( s%blocks(ib)%neighbors(i_dim,i_dir) .lt. &
                     NO_BLOCK_NEIGHBOR ) then
                     call wb_abort( "neighbor to block N1 in direction N2 and &
                                    &dimension N3 is negative", &
                                    MPI_ERR_COUNT, &
                                    (/ ib, i_dir, i_dim /) )
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

      do ib = 1_SP, s%nb
         call mpi_bcast( s%blocks(ib)%block_size, 1_MP, MPI_INTEGER, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%np, int(s%n_dim,MP), MPI_INTEGER, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%neighbors, int(s%n_dim*N_DIR,MP), &
            MPI_SP, WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%nx, int(s%n_dim,MP), MPI_SP, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%blocks(ib)%periods, int(s%n_dim,MP), MPI_LOGICAL, &
            WORLD_MASTER, MPI_COMM_WORLD, ierr )
      end do

      deallocate( np )
      deallocate( nx )
      deallocate( neighbors_l )
      deallocate( neighbors_u )
   end subroutine read_block_namelists

   subroutine setup_processes( s )
      integer(MP) :: assigned_processes, ierr, world_rank
      integer(SP) :: ib, i_dim, total_points = 0_SP
      type(MPI_Comm) :: comm_split
      type(WB_State), intent(inout) :: s

      allocate( s%nx(s%n_dim) )
      allocate( s%block_coords(s%n_dim) )
      allocate( s%neighbors(s%n_dim,N_DIR) )

      allocate( s%processes(0_MP:s%world_size-1_MP) )
      do world_rank = 0_MP, s%world_size-1_MP
         allocate( s%processes(world_rank)%block_coords(s%n_dim) )
         allocate( s%processes(world_rank)%nx(s%n_dim) )
      end do

      ib = 1_SP
      assigned_processes = 0_MP
      do world_rank = 0_MP, s%world_size-1_MP
         s%processes(world_rank)%ib = ib
         assigned_processes = assigned_processes + 1_MP
         if ( assigned_processes .eq. s%blocks(ib)%block_size ) then
            assigned_processes = 0_MP
            ib = ib + 1_SP
         end if
      end do

      s%ib = s%processes(s%world_rank)%ib
      call mpi_comm_split( MPI_COMM_WORLD, int(s%ib,MP), 0_MP, comm_split, &
         ierr )
      call mpi_cart_create( comm_split, int(s%n_dim,MP), s%blocks(s%ib)%np, &
         s%blocks(s%ib)%periods, s%blocks(s%ib)%reorder, s%comm_block, ierr )
      call mpi_comm_free( comm_split, ierr )
      call mpi_comm_rank( s%comm_block, s%block_rank, ierr )
      call mpi_comm_size( s%comm_block, s%block_size, ierr )
      call mpi_cart_coords( s%comm_block, s%block_rank, int(s%n_dim,MP), &
         s%block_coords, ierr )

      s%nx = s%blocks(s%ib)%nx / int(s%blocks(s%ib)%np,SP)
      do i_dim = 1_SP, s%n_dim
         if ( s%block_coords(i_dim) .eq. s%blocks(s%ib)%np(i_dim)-1_MP ) then
            s%nx(i_dim) = s%nx(i_dim) + modulo( s%blocks(s%ib)%nx(i_dim), &
               int(s%blocks(s%ib)%np(i_dim),SP) )
         end if
      end do

      s%processes(s%world_rank)%block_rank = s%block_rank
      s%processes(s%world_rank)%block_coords = s%block_coords
      s%processes(s%world_rank)%nx = s%nx
      do world_rank = 0_MP, s%world_size-1_MP
         call mpi_bcast( s%processes(world_rank)%block_rank, 1_MP, &
            MPI_INTEGER, world_rank, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%processes(world_rank)%block_coords, int(s%n_dim,MP), &
            MPI_INTEGER, world_rank, MPI_COMM_WORLD, ierr )
         call mpi_bcast( s%processes(world_rank)%nx, int(s%n_dim,MP), &
            MPI_SP, world_rank, MPI_COMM_WORLD, ierr )
      end do

      ! Check if the sum of the points in a block's processes equals the total
      ! number of points.
      call mpi_reduce( product(s%nx), total_points, int(s%n_dim,MP), MPI_SP, &
         MPI_SUM, BLOCK_MASTER, s%comm_block, ierr )
      if ( s%block_rank .eq. BLOCK_MASTER .and. &
         product(s%blocks(s%ib)%nx) .ne. total_points ) then
         call wb_abort( "total points in block N1 does not match sum of &
                        &points in individual processes", &
            MPI_ERR_SIZE, (/ s%ib /) )
      end if

      call identify_process_neighbors( s )
   end subroutine setup_processes

   subroutine wb_abort( message, errno, ints, floats )
      character(len=*), intent(in) :: message
      integer(MP), intent(in) :: errno
      integer(SP) :: i
      integer(MP) :: ierr
      integer(SP), dimension(:), optional, intent(in) :: ints
      real(FP), dimension(:), optional, intent(in) :: floats

      write (error_unit, "(A, A, A)") PROGRAM_NAME, ": ", message
      if ( present(ints) ) then
         do i = 1_SP, size(ints)
            write (error_unit, "(A, I1, A, I8)") "N", i, " = ", ints(i)
         end do
      end if
      if ( present(floats) ) then
         do i = 1_SP, size(floats)
            write (error_unit, "(A, I1, A, ES9.2)") "F", i, " = ", floats(i)
         end do
      end if
      call mpi_abort( MPI_COMM_WORLD, errno, ierr )
   end subroutine wb_abort

   subroutine write_block_information( f, s )
      integer, intent(in) :: f
      type(WB_State), intent(in) :: s
      integer(SP) :: ib, i_dim
      character(len=STRING_LENGTH) :: label

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (f,"(A)") "## Block information"
         write (f,"(A)") ""

         call write_string_table_entry( f, "`ib`", IB_COLUMN_WIDTH )
         call write_string_table_entry( f, "`block_size`", SIZE_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
            write (label,"(A, I1, A)") "`np(", i_dim, ")`"
            call write_string_table_entry( f, trim(label), NP_COLUMN_WIDTH )
         end do
         call write_string_table_entry( f, "points", POINTS_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
            write (label,"(A, I1, A)") "`nx(", i_dim, ")`"
            call write_string_table_entry( f, trim(label), NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. s%n_dim) )
         end do

         call write_rule_table_entry( f, IB_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_rule_table_entry( f, SIZE_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            call write_rule_table_entry( f, NP_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED )
         end do
         call write_rule_table_entry( f, POINTS_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            call write_rule_table_entry( f, NX_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED, end_row=(i_dim .eq. s%n_dim) )
         end do

         do ib = 1_SP, s%nb
            call write_integer_table_entry( f, ib, IB_COLUMN_WIDTH )
            call write_integer_table_entry( f, &
               int(s%blocks(ib)%block_size,SP), SIZE_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               call write_integer_table_entry( f, &
               int(s%blocks(ib)%np(i_dim),SP), NP_COLUMN_WIDTH )
            end do
            call write_integer_table_entry( f, &
               product(s%blocks(ib)%nx), POINTS_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               call write_integer_table_entry( f, &
               s%blocks(ib)%nx(i_dim), NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. s%n_dim) )
            end do
         end do

         write (f,"(A)") ""
      end if
   end subroutine write_block_information

   subroutine write_global_information( f, s )
      integer, intent(in) :: f
      type(WB_State), intent(in) :: s
      integer(MP) :: ierr, mpi_major_version_number, &
         mpi_minor_version_number, version_length
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

         write (f,"(A)", advance="no") "- Floating point precision `FP = "
         if ( FP .eq. real64 ) then
            write (f,"(A)") "real64`."
         else if ( FP .eq. real32 ) then
            write (f,"(A)") "real32`."
         else
            write (f,"(A)") "?`."
         end if
         write (f,"(A)") ""

         write (f,"(A)", advance="no") &
            "    - MPI floating point precision `MPI_FP = "
         if ( MPI_FP .eq. MPI_REAL8 ) then
            write (f,"(A)") "MPI_REAL8`."
         else if ( MPI_FP .eq. MPI_REAL4 ) then
            write (f,"(A)") "MPI_REAL4`."
         else if ( MPI_FP .eq. MPI_REAL ) then
            write (f,"(A)") "MPI_REAL`."
         else
            write (f,"(A)") "?`."
         end if
         write (f,"(A)") ""

         write (f,"(A)", advance="no") "- Signed integer precision `SP = "
         if ( SP .eq. int64 ) then
            write (f,"(A)") "int64`."
         else if ( SP .eq. int32 ) then
            write (f,"(A)") "int32`."
         else
            write (f,"(A)") "?`."
         end if
         write (f,"(A)") ""

         write (f,"(A)", advance="no") &
            "    - MPI signed integer precision `MPI_SP = "
         if ( MPI_SP .eq. MPI_INTEGER8 ) then
            write (f,"(A)") "MPI_INTEGER8`."
         else if ( MPI_SP .eq. MPI_INTEGER4 ) then
            write (f,"(A)") "MPI_INTEGER4`."
         else if ( MPI_SP .eq. MPI_INTEGER ) then
            write (f,"(A)") "MPI_INTEGER`."
         else
            write (f,"(A)") "?`."
         end if
         write (f,"(A)") ""

         write (f,"(A)", advance="no") "- Signed integer precision `MP = "
         if ( MP .eq. int64 ) then
            write (f,"(A)") "int64`."
         else if ( MP .eq. int32 ) then
            write (f,"(A)") "int32`."
         else
            write (f,"(A)") "?`."
         end if
         write (f,"(A)") ""
         write (f,"(A)") "    - (for integers only used as MPI variables)"
         write (f,"(A)") ""

         write (f,"(A)", advance="no") &
            "    - MPI signed integer precision `MPI_MP = "
         if ( MPI_MP .eq. MPI_INTEGER8 ) then
            write (f,"(A)") "MPI_INTEGER8`."
         else if ( MPI_MP .eq. MPI_INTEGER4 ) then
            write (f,"(A)") "MPI_INTEGER4`."
         else if ( MPI_MP .eq. MPI_INTEGER ) then
            write (f,"(A)") "MPI_INTEGER`."
         else
            write (f,"(A)") "?`."
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
      integer(MP) :: ierr, processor_length, world_rank
      integer(SP) :: i_dim
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
         do i_dim = 1_SP, s%n_dim
            write (label,"(A, I1, A)") "(", i_dim, ")"
            call write_string_table_entry( f, trim(label), &
               COORDS_COLUMN_WIDTH )
         end do
         call write_string_table_entry( f, "points", POINTS_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
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
         do i_dim = 1_SP, s%n_dim
            call write_rule_table_entry( f, COORDS_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED )
         end do
         call write_rule_table_entry( f, POINTS_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            call write_rule_table_entry( f, NX_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED, end_row=(i_dim .eq. s%n_dim) )
         end do
      end if

      do world_rank = 0_MP, s%world_size-1_MP
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            call write_integer_table_entry( f, int(s%world_rank,SP), &
               RANK_COLUMN_WIDTH )
            call write_string_table_entry(  f, trim(processor_name), &
               HOSTNAME_COLUMN_WIDTH )
            call write_integer_table_entry( f, s%ib, &
               IB_COLUMN_WIDTH )
            call write_integer_table_entry( f, int(s%block_rank,SP), &
               RANK_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               call write_integer_table_entry( f, &
               int(s%block_coords(i_dim),SP), COORDS_COLUMN_WIDTH )
            end do
            call write_integer_table_entry( f, product(s%nx), &
               POINTS_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               call write_integer_table_entry( f, s%nx(i_dim), &
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
      integer(MP) :: ierr, world_rank, neighbor
      integer(SP) :: i_dim, i_dir, j_dir, face_count=0
      integer(SP), dimension(N_DIR) :: dirs = (/ LOWER_DIR, UPPER_DIR /)
      character(len=STRING_LENGTH) :: label

      if ( s%world_rank .eq. WORLD_MASTER ) then
         write (*,"(A)") "## Process neighbors"
         write (*,"(A)") ""

         call write_string_table_entry( f, "`world_rank`", RANK_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
            do i_dir = 1_SP, N_DIR
               j_dir = dirs(i_dir)
               face_count = face_count + 1
               if ( j_dir .eq. LOWER_DIR ) then
                  write (label,"(I1, A)") i_dim, "L"
               else
                  write (label,"(I1, A)") i_dim, "U"
               end if
               call write_string_table_entry( f, trim(label), &
                  RANK_COLUMN_WIDTH, &
                  end_row=( face_count .eq. int(s%n_dim*N_DIR,MP) ) )
            end do
         end do

         face_count = 0_SP
         call write_rule_table_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            do i_dir = 1_SP, N_DIR
               face_count = face_count + 1_SP
               call write_rule_table_entry( f, RANK_COLUMN_WIDTH, &
                  alignment=RIGHT_ALIGNED, &
                  end_row=( face_count .eq. int(s%n_dim*N_DIR,MP) ) )
            end do
         end do
      end if

      do world_rank = 0_MP, s%world_size-1_MP
         face_count = 0_SP
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            call write_integer_table_entry( f, int(s%world_rank,SP), &
               RANK_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               do i_dir = 1_SP, N_DIR
                  j_dir = dirs(i_dir)
                  face_count = face_count + 1_SP
                  neighbor = s%neighbors(i_dim,j_dir)
                  if ( neighbor .eq. MPI_PROC_NULL ) then
                     call write_string_table_entry( f, "-", &
                        RANK_COLUMN_WIDTH, &
                        end_row=( face_count .eq. int(s%n_dim*N_DIR,MP) ) )
                  else
                     call write_integer_table_entry( f, int(neighbor,SP), &
                        RANK_COLUMN_WIDTH, &
                        end_row=( face_count .eq. int(s%n_dim*N_DIR,MP) ) )
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
