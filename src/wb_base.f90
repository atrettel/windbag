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
   use wb_exit
   use wb_text
   implicit none

   private

   public find_input_file, print_initial_information, wb_subdomain_construct, &
      wb_subdomain_destroy

   integer(MP), public, parameter :: BLOCK_MASTER = 0_MP
   integer(MP), public, parameter :: WORLD_MASTER = 0_MP

   integer(SP), public, parameter :: DEFAULT_BLOCK_NEIGHBOR = 1_SP
   integer(SP), public, parameter :: DEFAULT_BLOCK_NUMBER   = 1_SP
   integer(SP), public, parameter ::      NO_BLOCK_NEIGHBOR = 0_SP

   integer(SP), public, parameter :: DEFAULT_NUMBER_OF_BLOCKS = 1_SP

   integer(SP), public, parameter :: DEFAULT_NUMBER_OF_COMPONENTS = 1_SP

   integer(SP), public, parameter :: DEFAULT_NUMBER_OF_DIMENSIONS = 3_SP
   integer(SP), public, parameter ::     MIN_NUMBER_OF_DIMENSIONS = 1_SP
   integer(SP), public, parameter ::     MAX_NUMBER_OF_DIMENSIONS = 3_SP

   integer(SP), public, parameter :: NUMBER_OF_DIRECTIONS = 2_SP
   integer(SP), public, parameter ::     LOWER_DIRECTION  = 1_SP
   integer(SP), public, parameter ::     UPPER_DIRECTION  = 2_SP

   integer(SP), public, parameter :: NO_FIELD = 0_SP

   integer(SP), public, parameter :: DEFAULT_NUMBER_OF_GHOST_POINTS = 3_SP

   integer(SP), public, parameter :: NUMBER_OF_PHASES = 1_SP

   integer(SP), public, parameter :: DEFAULT_NX = 0_SP

   integer(MP), public, parameter :: DEFAULT_NP = 0_MP

   logical, public, parameter :: DEFAULT_REORDER = .false.

   character(len=*), public, parameter ::      PROGRAM_NAME = "windbag"
   character(len=*), public, parameter ::           VERSION = "0.0.0"
   character(len=*), public, parameter :: DEFAULT_CASE_NAME = "casename"

   type, private :: WB_Block
      private
      integer(MP), dimension(:), allocatable :: np
      integer(SP), dimension(:), allocatable :: nx
      integer(SP), dimension(:,:), allocatable :: neighbors
      logical, dimension(:), allocatable :: periods
      logical :: reorder
   end type WB_Block

   type, private :: WB_Process
      private
      integer(MP) :: block_rank
      integer(MP), dimension(:), allocatable :: block_coords
      integer(SP), dimension(:), allocatable :: nx
      integer(SP) :: ib
   end type WB_Process

   type, public :: WB_Subdomain
      private
      character(len=:), allocatable :: case_name
      integer(MP) :: block_rank, block_size
      integer(MP) :: world_rank, world_size
      integer(SP) :: ib, nb
      integer(SP) :: n_dim
      integer(SP) :: nc
      integer(SP) :: nf
      integer(SP) :: ng
      integer(SP) :: nv
      integer(SP) :: i_iter
      integer(SP), dimension(:), allocatable :: nx
      integer(MP), dimension(:), allocatable :: block_coords
      integer(MP), dimension(:,:), allocatable :: neighbors
      real(FP) :: t
      real(FP), dimension(:,:,:,:), allocatable :: f
      type(MPI_Comm) :: comm_block
      type(WB_Block) :: local_block
   end type WB_Subdomain

   interface neighbor
      module procedure wb_block_neighbor, wb_subdomain_neighbor
   end interface neighbor

   interface num_blocks
      module procedure wb_subdomain_total_blocks
   end interface num_blocks

   interface num_components
      module procedure wb_subdomain_components
   end interface num_components

   interface num_dimensions
      module procedure wb_subdomain_dimensions
   end interface num_dimensions

   interface num_dimensions_mp
      module procedure wb_subdomain_dimensions_mp
   end interface num_dimensions_mp

   interface num_ghost_points
      module procedure wb_subdomain_ghost_points
   end interface num_ghost_points

   interface num_variables
      module procedure wb_subdomain_variables
   end interface num_variables

   interface num_points
      module procedure wb_block_points, wb_subdomain_points
   end interface num_points

   interface total_points
      module procedure wb_block_total_points, wb_subdomain_total_points
   end interface total_points

   interface wb_subdomain_construct
      module procedure wb_subdomain_construct_namelist, &
         wb_subdomain_construct_variables
   end interface wb_subdomain_construct
contains
   subroutine assign_blocks( blocks, block_assignments, world_size )
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      integer(SP), dimension(:), allocatable, intent(inout) :: &
         block_assignments
      integer(MP) :: assigned_processes, world_rank, world_size
      integer(SP) :: ib

      ib = 1_SP
      assigned_processes = 0_MP
      do world_rank = 0_MP, world_size-1_MP
         block_assignments(world_rank) = ib
         assigned_processes = assigned_processes + 1_MP
         if ( assigned_processes .eq. wb_block_size( blocks(ib) ) ) then
            assigned_processes = 0_MP
            ib = ib + 1_SP
         end if
      end do
   end subroutine assign_blocks

   subroutine calculate_number_of_variables( sd )
      type(WB_Subdomain), intent(inout) :: sd

      sd%nv = phase_rule( num_components(sd), NUMBER_OF_PHASES ) + &
              dimension_rule( num_dimensions(sd) )
   end subroutine calculate_number_of_variables

   subroutine check_block_bounds( ib, nb )
      integer(SP), intent(in) :: ib, nb

      if ( ib .lt. 1_SP .or. ib .gt. nb ) then
         call wb_abort( "block N1 is out of acceptable range [N2, N3]", &
            EXIT_DATAERR, (/ ib, 1_SP, nb /) )
      end if
   end subroutine check_block_bounds

   subroutine check_block_dimension_arrays( blocks, nb, n_dim, ng )
      integer(SP), intent(in) :: nb, n_dim, ng
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      integer(SP) :: ib, i_dir, i_dim

      do ib = 1_SP, nb
         do i_dim = 1_SP, n_dim
            if ( wb_block_processes( blocks(ib), i_dim ) .lt. 1_MP ) then
               call wb_abort( "number of processes in direction N1 of &
                              &block N2 is less than 1", &
                              EXIT_DATAERR, (/ i_dim, ib /) )
            end if
            if ( num_points( blocks(ib), i_dim ) .lt. ng ) then
               call wb_abort( "number of points in direction N1 of block &
                              &N2 is less than number of ghost points N3", &
                              EXIT_DATAERR, (/ i_dim, ib, ng /) )
            end if
            do i_dir = 1_SP, NUMBER_OF_DIRECTIONS
               if ( neighbor( blocks(ib), i_dim, i_dir ) .lt. &
                  NO_BLOCK_NEIGHBOR ) then
                  call wb_abort( "neighbor to block N1 in direction N2 and &
                                 &dimension N3 is negative", &
                                 EXIT_DATAERR, (/ ib, i_dir, i_dim /) )
               end if
            end do
         end do
      end do
   end subroutine check_block_dimension_arrays

   subroutine check_block_neighbors( blocks, nb, n_dim )
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      integer(SP), intent(in) :: nb, n_dim
      integer(SP) :: ib, i_dim, j_dim, neighbor_l, neighbor_u

      ! Ensure that lower and upper pairs exist, and that their number of
      ! points and processes match.  Since this is just checking and not
      ! calculation, it must occur on the world master.  It is impossible to
      ! have each block master check this for their own blocks since the block
      ! communicators do not exist yet.  If that were possible, it would
      ! eliminate the loop over all blocks.
      do ib = 1_SP, nb
         do i_dim = 1_SP, n_dim
            neighbor_l = neighbor( blocks(ib), i_dim, LOWER_DIRECTION )
            if ( neighbor_l .ne. NO_BLOCK_NEIGHBOR ) then
               neighbor_u = neighbor( blocks(neighbor_l), i_dim, UPPER_DIRECTION )
               if ( ib .ne. neighbor_u ) then
                  call wb_abort( "lower face of block N1 does not neighbor &
                                 &upper face of block N2 in direction N3", &
                     EXIT_DATAERR, &
                     (/ ib, neighbor_l, i_dim /) )
               else
                  do j_dim = 1_SP, n_dim
                     if ( j_dim .ne. i_dim .and. &
                        wb_block_processes( blocks(ib),         j_dim ) .ne. &
                        wb_block_processes( blocks(neighbor_l), j_dim ) ) then
                        call wb_abort( "face in direction N1 shared by &
                                       &blocks N2 and N3 does not match &
                                       &processes in direction N4", &
                           EXIT_DATAERR, &
                           (/ i_dim, ib, neighbor_l, j_dim /) )
                     end if
                     if ( j_dim .ne. i_dim .and. &
                        num_points( blocks(ib),         j_dim ) .ne. &
                        num_points( blocks(neighbor_l), j_dim ) ) then
                        call wb_abort( "face in direction N1 shared by &
                                       &blocks N2 and N3 does not match &
                                       &points in direction N4", &
                           EXIT_DATAERR, &
                           (/ i_dim, ib, neighbor_l, j_dim /) )
                     end if
                  end do
               end if
            end if
         end do
      end do
   end subroutine check_block_neighbors

   subroutine check_block_total_points( sd )
      integer(SP) :: points_in_block, points_in_processes
      integer(MP) :: ierr
      type(WB_Subdomain), intent(in) :: sd
      type(MPI_Comm) :: comm_block
      type(WB_Block) :: local_block

      call wb_subdomain_block_communicator( sd, comm_block )
      call wb_subdomain_local_block( sd, local_block )

      points_in_block = total_points(local_block)
      points_in_processes = 0_SP
      call mpi_reduce( total_points(sd), points_in_processes, &
         num_dimensions_mp(sd), MPI_SP, MPI_SUM, BLOCK_MASTER, &
         comm_block, ierr )
      if ( wb_subdomain_is_block_master(sd) .and. &
         points_in_block .ne. points_in_processes ) then
         call wb_abort( "total points in block N1 (N2) does not match sum of &
                        &points in individual processes (N3)", &
            EXIT_FAILURE, &
            (/ wb_subdomain_block_number(sd), points_in_block, &
            points_in_processes /) )
      end if
   end subroutine check_block_total_points

   subroutine check_block_world_size( blocks, nb )
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      integer(SP), intent(in) :: nb
      integer(SP) :: ib
      integer(MP) :: ierr, world_size, world_size_from_blocks

      world_size_from_blocks = 0_MP
      do ib = 1_SP, nb
         world_size_from_blocks = world_size_from_blocks + &
            wb_block_size( blocks(ib) )
      end do
      call mpi_comm_size( MPI_COMM_WORLD, world_size, ierr )
      if ( world_size_from_blocks .ne. world_size ) then
         call wb_abort( &
            "size of block domain decomposition (N1) does not match &
            &world size (N2)", EXIT_DATAERR, &
            int( (/ world_size_from_blocks, world_size /), SP ) )
      end if
   end subroutine check_block_world_size

   subroutine check_general_variables( nb, nc, n_dim, ng, world_size )
      integer(SP), intent(in) :: nb, nc, n_dim, ng
      integer(MP), intent(in) :: world_size

      if ( nb .lt. 1_SP  .or. nb .gt. int(world_size,SP) ) then
         ! The second condition is a feature of the code and not a bug.  It
         ! allows the code to treat communication between blocks as the same
         ! as communication between processes, but that plan only works if
         ! each block uses at least one process.
         call wb_abort( &
            "number of blocks must be in interval [N1, N2]", &
            EXIT_DATAERR, (/ 1_SP, int(world_size,SP) /) )
      end if
      if ( nc .lt. 1_SP ) then
         call wb_abort( &
            "number of components must be at least 1", EXIT_DATAERR )
      end if
      if ( ng .lt. 1_SP ) then
         call wb_abort( "number of ghost points is less than 1", &
            EXIT_DATAERR )
      end if
      if ( n_dim .lt. MIN_NUMBER_OF_DIMENSIONS .or. &
           n_dim .gt. MAX_NUMBER_OF_DIMENSIONS ) then
         call wb_abort( &
            "number of dimensions must be in interval [N1, N2]", &
            EXIT_DATAERR, (/ MIN_NUMBER_OF_DIMENSIONS, &
                             MAX_NUMBER_OF_DIMENSIONS /) )
      end if
   end subroutine check_general_variables

   subroutine check_points( sd )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: i_dim
      integer(MP) :: ierr, world_rank

      do world_rank = 0_MP, wb_subdomain_world_size(sd)-1_MP
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         do i_dim = 1, num_dimensions(sd)
            if ( num_points(sd,i_dim) .lt. num_ghost_points(sd) .and. &
               wb_subdomain_world_rank(sd) .eq. world_rank ) then
            call wb_abort( &
               "number of points N1 in dimension N2 is less than the number &
               &of ghost points N3 for subdomain N4", &
               EXIT_FAILURE, (/ num_points(sd,i_dim), i_dim, &
               num_ghost_points(sd), int(wb_subdomain_world_rank(sd),SP) /) )
            end if
         end do
      end do
   end subroutine check_points

   subroutine decompose_domain( sd, blocks, processes )
      integer(MP) :: ierr, world_rank
      integer(SP) :: i_dim
      type(MPI_Comm) :: comm_split
      type(WB_Subdomain), intent(inout) :: sd
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      type(WB_Process), dimension(:), allocatable, intent(out) :: processes
      integer(SP), dimension(:), allocatable :: block_nx, block_assignments, &
         block_neighbors_l, block_neighbors_u
      integer(MP), dimension(:), allocatable :: block_np
      logical, dimension(:), allocatable :: block_periods

      allocate( block_assignments(0_MP:wb_subdomain_world_size(sd)-1_MP) )
      call assign_blocks( blocks, block_assignments, &
         wb_subdomain_world_size(sd) )

      sd%ib = block_assignments(wb_subdomain_world_rank(sd))
      allocate( block_neighbors_l(num_dimensions(sd)), &
                block_neighbors_u(num_dimensions(sd)), &
                         block_np(num_dimensions(sd)), &
                         block_nx(num_dimensions(sd)), &
                    block_periods(num_dimensions(sd)) )
      call wb_block_neighbors_vector( blocks(wb_subdomain_block_number(sd)), &
         block_neighbors_l, LOWER_DIRECTION )
      call wb_block_neighbors_vector( blocks(wb_subdomain_block_number(sd)), &
         block_neighbors_u, UPPER_DIRECTION )
      call wb_block_processes_vector( blocks(wb_subdomain_block_number(sd)), &
         block_np )
      call wb_block_points_vector( blocks(wb_subdomain_block_number(sd)), &
         block_nx )
      call wb_block_periods_vector( blocks(wb_subdomain_block_number(sd)), &
         block_periods )
      call wb_block_construct( sd%local_block, num_dimensions(sd), &
         block_np, block_nx, block_neighbors_l, block_neighbors_u, &
         wb_subdomain_block_number(sd) )

      call mpi_comm_split( MPI_COMM_WORLD, &
         wb_subdomain_block_number_mp(sd), 0_MP, comm_split, ierr )
      call mpi_cart_create( comm_split, num_dimensions_mp(sd), &
         block_np, block_periods, wb_block_reorder(sd%local_block), &
         sd%comm_block, ierr )
      call mpi_comm_free( comm_split, ierr )
      call mpi_comm_rank( sd%comm_block, sd%block_rank, ierr )
      call mpi_comm_size( sd%comm_block, sd%block_size, ierr )
      call mpi_cart_coords( sd%comm_block, wb_subdomain_block_rank(sd), &
         num_dimensions_mp(sd), sd%block_coords, ierr )

      sd%nx = block_nx / int(block_np,SP)
      do i_dim = 1_SP, num_dimensions(sd)
         if ( sd%block_coords(i_dim) .eq. &
            wb_block_processes( sd%local_block, i_dim ) - 1_MP ) then
            sd%nx(i_dim) = sd%nx(i_dim) + modulo( &
               num_points( sd%local_block, i_dim ), &
               int( wb_block_processes( sd%local_block, i_dim ), SP ) )
         end if
      end do

      allocate( processes(0_MP:wb_subdomain_world_size(sd)-1_MP) )
      do world_rank = 0_MP, wb_subdomain_world_size(sd)-1_MP
         if ( world_rank .eq. wb_subdomain_world_rank(sd) ) then
            call wb_process_construct( processes(world_rank), &
               num_dimensions(sd), block_assignments(world_rank), &
               world_rank, wb_subdomain_block_rank(sd), sd%block_coords, &
               sd%nx )
         else
            call wb_process_construct( processes(world_rank), &
               num_dimensions(sd), block_assignments(world_rank), &
               world_rank )
         end if
      end do

      deallocate( block_assignments, &
                  block_neighbors_l, &
                  block_neighbors_u, &
                  block_np,          &
                  block_nx,          &
                  block_periods )
   end subroutine decompose_domain

   pure function dimension_rule( n_dim ) result( f )
      integer(SP), intent(in) :: n_dim
      integer(SP) :: f

      f = 2_SP * n_dim + n_dim**2_SP
   end function dimension_rule

   subroutine find_input_file( filename )
      character(len=STRING_LENGTH), intent(out)  :: filename
      integer :: argc, filename_length, ierr
      integer(MP) :: ierr_mpi, world_rank
      logical :: file_exists

      call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr_mpi )
      if ( world_rank .eq. WORLD_MASTER ) then
         argc = command_argument_count()
         if ( argc .eq. 0 ) then
            call wb_abort( "no input file given", EXIT_USAGE )
         end if
         call get_command_argument( 1, filename, filename_length, ierr )
         inquire( file=filename, exist=file_exists, iostat=ierr )
         if ( file_exists .eqv. .false. ) then
            call wb_abort( "input file does not exist", EXIT_NOINPUT )
         end if
      end if
      call mpi_barrier( MPI_COMM_WORLD, ierr )
   end subroutine find_input_file

   subroutine identify_process_neighbors( sd, blocks, processes )
      integer(MP) :: ierr, world_rank
      integer(SP) :: block_neighbor, i_dir, i_dim
      integer(MP), dimension(NUMBER_OF_DIRECTIONS) :: block_ranks
      integer(MP), dimension(:), allocatable :: block_coords, &
         process_block_coords
      type(WB_Subdomain), intent(inout) :: sd
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      type(WB_Process), dimension(:), allocatable, intent(in) :: processes

      allocate( block_coords(num_dimensions(sd)), &
        process_block_coords(num_dimensions(sd)) )
      do i_dim = 1_SP, num_dimensions(sd)
         call mpi_cart_shift( sd%comm_block, int(i_dim,MP)-1_MP, 1_MP, &
            block_ranks(LOWER_DIRECTION), block_ranks(UPPER_DIRECTION), ierr )
         do i_dir = 1_SP, NUMBER_OF_DIRECTIONS
            if ( block_ranks(i_dir) .ne. MPI_PROC_NULL ) then
               do world_rank = 0_MP, wb_subdomain_world_size(sd)-1_MP
                  if ( wb_process_block_number( processes(world_rank) ) .eq. &
                     wb_subdomain_block_number(sd) .and. &
                     wb_process_block_rank( processes(world_rank) ) .eq. &
                     block_ranks(i_dir) ) then
                     exit
                  end if
               end do
               sd%neighbors(i_dim,i_dir) = world_rank
            else
               block_neighbor = neighbor( sd%local_block, i_dim, i_dir )
               if ( block_neighbor .eq. NO_BLOCK_NEIGHBOR ) then
                  sd%neighbors(i_dim,i_dir) = MPI_PROC_NULL
               else
                  ! This block neighbors another block.  This neighboring block
                  ! sits on the opposite side of the current block.  Both
                  ! processes share the same block coordinates except for the
                  ! current spatial dimension i_dim.  The block coordinates for
                  ! dimension i_dim would be opposites for each.
                  call wb_subdomain_block_coords_vector( sd, block_coords )
                  if ( i_dir .eq. LOWER_DIRECTION ) then
                     block_coords(i_dim) = &
                        wb_block_processes(blocks(block_neighbor),i_dim)-1_MP
                  else
                     block_coords(i_dim) = 0_MP
                  end if
                  do world_rank = 0_MP, wb_subdomain_world_size(sd)-1_MP
                     call wb_process_block_coords( processes(world_rank), &
                        process_block_coords )
                     if ( wb_process_block_number( processes(world_rank) ) &
                        .eq. block_neighbor .and. all( process_block_coords &
                        .eq. block_coords ) ) then
                        exit
                     end if
                  end do
                  sd%neighbors(i_dim,i_dir) = world_rank
               end if
            end if
         end do
      end do
      deallocate( block_coords, process_block_coords )
   end subroutine identify_process_neighbors

   pure function phase_rule( c, p ) result( f )
      integer(SP), intent(in) :: c, p
      integer(SP) :: f

      f = c - p + 2_SP
   end function phase_rule

   subroutine print_initial_information( sd )
      type(WB_Subdomain), intent(in) :: sd

      call write_environment_information( output_unit, sd )
      call write_scalar_variables(        output_unit, sd )
      call write_block_information(       output_unit, sd )
      call write_subdomain_information(   output_unit, sd )
      call write_subdomain_neighbors(     output_unit, sd )
   end subroutine print_initial_information

   subroutine read_block_namelists( filename, nb, n_dim, blocks )
      character(len=STRING_LENGTH), intent(in) :: filename
      type(WB_Block), dimension(:), allocatable, intent(out) :: blocks
      integer(SP), intent(in) :: nb, n_dim
      integer :: file_unit
      integer(MP) :: ierr, world_rank, world_size
      integer(SP) :: ib, jb
      integer(MP), dimension(:), allocatable :: np
      integer(SP), dimension(:), allocatable :: nx, neighbors_l, neighbors_u
      namelist /block/ ib, np, nx, neighbors_l, neighbors_u

      allocate( np(n_dim) )
      allocate( nx(n_dim) )
      allocate( neighbors_l(n_dim) )
      allocate( neighbors_u(n_dim) )

      ib             = DEFAULT_BLOCK_NUMBER
      np(:)          = DEFAULT_NP
      nx(:)          = DEFAULT_NX
      neighbors_l(:) = DEFAULT_BLOCK_NEIGHBOR
      neighbors_u(:) = DEFAULT_BLOCK_NEIGHBOR

      allocate( blocks(nb) )
      call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr )
      call mpi_comm_size( MPI_COMM_WORLD, world_size, ierr )
      if ( world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
      end if
      do jb = 1_SP, nb
         if ( world_rank .eq. WORLD_MASTER ) then
            read( unit=file_unit, nml=block )
            call check_block_bounds( ib, nb )
         end if

         call mpi_bcast( ib, 1_MP, MPI_INTEGER, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( np, int(n_dim,MP), MPI_INTEGER, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( nx, int(n_dim,MP), MPI_SP, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( neighbors_l, int(n_dim,MP), MPI_SP, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )
         call mpi_bcast( neighbors_u, int(n_dim,MP), MPI_SP, WORLD_MASTER, &
            MPI_COMM_WORLD, ierr )

         call wb_block_construct( blocks(ib), n_dim, np, nx, neighbors_l, &
            neighbors_u, ib )
      end do
      if ( world_rank .eq. WORLD_MASTER ) then
         close( unit=file_unit )
      end if

      deallocate( np )
      deallocate( nx )
      deallocate( neighbors_l )
      deallocate( neighbors_u )
   end subroutine read_block_namelists

   subroutine read_general_namelist( filename, case_name, nb, nc, n_dim, ng )
      character(len=STRING_LENGTH), intent(in) :: filename
      character(len=STRING_LENGTH), intent(out) :: case_name
      integer(SP), intent(out) :: nb, nc, n_dim, ng
      integer :: file_unit
      integer(MP) :: ierr, world_rank
      namelist /general/ case_name, nb, nc, ng, n_dim

      case_name = DEFAULT_CASE_NAME
      nb        = DEFAULT_NUMBER_OF_BLOCKS
      nc        = DEFAULT_NUMBER_OF_COMPONENTS
      n_dim     = DEFAULT_NUMBER_OF_DIMENSIONS
      ng        = DEFAULT_NUMBER_OF_GHOST_POINTS

      call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr )
      if ( world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         read( unit=file_unit, nml=general )
         close( unit=file_unit )
      end if

      call mpi_bcast( case_name, int(STRING_LENGTH,MP), MPI_CHARACTER, &
         WORLD_MASTER, MPI_COMM_WORLD, ierr )
      call mpi_bcast( nb,    1_MP, MPI_SP, WORLD_MASTER, MPI_COMM_WORLD, ierr )
      call mpi_bcast( nc,    1_MP, MPI_SP, WORLD_MASTER, MPI_COMM_WORLD, ierr )
      call mpi_bcast( n_dim, 1_MP, MPI_SP, WORLD_MASTER, MPI_COMM_WORLD, ierr )
      call mpi_bcast( ng,    1_MP, MPI_SP, WORLD_MASTER, MPI_COMM_WORLD, ierr )
   end subroutine read_general_namelist

   subroutine wb_block_construct( blk, n_dim, np, nx, neighbors_l, &
      neighbors_u, ib )
      type(WB_Block), intent(inout) :: blk
      integer(SP), intent(in) :: n_dim
      integer(MP), dimension(:), allocatable, optional, intent(in) :: np
      integer(SP), dimension(:), allocatable, optional, intent(in) :: nx, &
         neighbors_l, neighbors_u
      integer(SP), optional, intent(in) :: ib
      integer(SP) :: i_dim

      allocate( blk%np(n_dim) )
      allocate( blk%nx(n_dim) )
      allocate( blk%neighbors(n_dim,NUMBER_OF_DIRECTIONS) )
      allocate( blk%periods(n_dim) )

      blk%np(:)          = DEFAULT_NP
      blk%nx(:)          = DEFAULT_NX
      blk%neighbors(:,:) = DEFAULT_BLOCK_NEIGHBOR
      blk%periods(:)     = .false.

      if ( present(np) ) then
         blk%np         = np
      end if
      if ( present(nx) ) then
         blk%nx = nx
      end if
      if ( present(neighbors_l) ) then
         blk%neighbors(:,LOWER_DIRECTION) = neighbors_l
      end if
      if ( present(neighbors_u) ) then
         blk%neighbors(:,UPPER_DIRECTION) = neighbors_u
      end if
      blk%reorder = DEFAULT_REORDER

      if ( present(ib) ) then
         do i_dim = 1_SP, n_dim
            if ( blk%neighbors(i_dim,LOWER_DIRECTION) .eq. ib .and. &
                 blk%neighbors(i_dim,UPPER_DIRECTION) .eq. ib ) then
               blk%periods(i_dim) = .true.
            else
               blk%periods(i_dim) = .false.
            end if
         end do
      end if
   end subroutine wb_block_construct

   subroutine wb_block_destroy( blk )
      type(WB_Block), intent(inout) :: blk

      deallocate( blk%np )
      deallocate( blk%nx )
      deallocate( blk%neighbors )
      deallocate( blk%periods )
   end subroutine wb_block_destroy

   function wb_block_neighbor( blk, i_dim, i_dir ) result( neighbor )
      type(WB_Block), intent(in) :: blk
      integer(SP), intent(in) :: i_dim, i_dir
      integer(SP) :: neighbor

      neighbor = blk%neighbors(i_dim,i_dir)
   end function wb_block_neighbor

   subroutine wb_block_neighbors_vector( blk, neighbors, i_dir )
      type(WB_Block), intent(in) :: blk
      integer(SP), intent(in) :: i_dir
      integer(SP), dimension(:), allocatable, intent(inout) :: neighbors

      neighbors = blk%neighbors(:,i_dir)
   end subroutine wb_block_neighbors_vector

   subroutine wb_block_periods_vector( blk, periods )
      type(WB_Block), intent(in) :: blk
      logical, dimension(:), allocatable, intent(inout) :: periods

      periods = blk%periods
   end subroutine wb_block_periods_vector

   function wb_block_points( blk, i_dim ) result( nx )
      type(WB_Block), intent(in) :: blk
      integer(SP), intent(in) :: i_dim
      integer(SP) :: nx

      nx = blk%nx(i_dim)
   end function wb_block_points

   subroutine wb_block_points_vector( blk, nx )
      type(WB_Block), intent(in) :: blk
      integer(SP), dimension(:), allocatable, intent(inout) :: nx

      nx = blk%nx
   end subroutine wb_block_points_vector

   function wb_block_processes( blk, i_dim ) result( np )
      type(WB_Block), intent(in) :: blk
      integer(SP), intent(in) :: i_dim
      integer(MP) :: np

      np = blk%np(i_dim)
   end function wb_block_processes

   subroutine wb_block_processes_vector( blk, np )
      type(WB_Block), intent(in) :: blk
      integer(MP), dimension(:), allocatable, intent(inout) :: np

      np = blk%np
   end subroutine wb_block_processes_vector

   function wb_block_reorder( blk ) result( reorder )
      type(WB_Block), intent(in) :: blk
      logical :: reorder

      reorder = blk%reorder
   end function wb_block_reorder

   function wb_block_size( blk ) result( block_size )
      type(WB_Block), intent(in) :: blk
      integer(MP) :: block_size

      block_size = product(blk%np)
   end function wb_block_size

   function wb_block_total_points( blk ) result( points_in_block )
      type(WB_Block), intent(in) :: blk
      integer(SP) :: points_in_block

      points_in_block = product(blk%nx)
   end function wb_block_total_points

   function wb_process_block_number( process ) result( block_number )
      type(WB_Process), intent(in) :: process
      integer(SP) :: block_number

      block_number = process%ib
   end function wb_process_block_number

   function wb_process_block_rank( process ) result( block_rank )
      type(WB_Process), intent(in) :: process
      integer(MP) :: block_rank

      block_rank = process%block_rank
   end function wb_process_block_rank

   subroutine wb_process_block_coords( process, block_coords )
      type(WB_Process), intent(in) :: process
      integer(MP), dimension(:), allocatable, intent(out) :: block_coords

      block_coords = process%block_coords
   end subroutine wb_process_block_coords

   subroutine wb_process_construct( process, n_dim, ib, world_rank, &
      block_rank, block_coords, nx  )
      type(WB_Process), intent(inout) :: process
      integer(SP), intent(in) :: n_dim, ib
      integer(MP), intent(in) :: world_rank
      integer(MP), optional, intent(in) :: block_rank
      integer(MP), dimension(:), allocatable, optional, intent(in) :: &
         block_coords
      integer(SP), dimension(:), allocatable, optional, intent(in) :: nx
      integer(MP) :: ierr

      allocate( process%block_coords(n_dim) )
      allocate( process%nx(n_dim) )

      process%ib = ib
      if ( present(block_rank) ) then
         process%block_rank = block_rank
      end if
      if ( present(block_coords) ) then
         process%block_coords = block_coords
      end if
      if ( present(nx) ) then
         process%nx = nx
      end if

      call mpi_bcast( process%block_rank, 1_MP, MPI_INTEGER, world_rank, &
         MPI_COMM_WORLD, ierr )
      call mpi_bcast( process%block_coords, int(n_dim,MP), MPI_INTEGER, &
         world_rank, MPI_COMM_WORLD, ierr )
      call mpi_bcast( process%nx, int(n_dim,MP), MPI_SP, world_rank, &
         MPI_COMM_WORLD, ierr )
   end subroutine wb_process_construct

   subroutine wb_process_destroy( process )
      type(WB_Process), intent(inout) :: process

      deallocate( process%block_coords )
      deallocate( process%nx )
   end subroutine wb_process_destroy

   subroutine wb_subdomain_block_communicator( sd, comm_block )
      type(WB_Subdomain), intent(in) :: sd
      type(MPI_Comm), intent(inout) :: comm_block

      comm_block = sd%comm_block
   end subroutine wb_subdomain_block_communicator

   function wb_subdomain_block_coord( sd, i_dim ) result( block_coord )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: i_dim
      integer(MP) :: block_coord

      block_coord = sd%block_coords(i_dim)
   end function wb_subdomain_block_coord

   subroutine wb_subdomain_block_coords_vector( sd, block_coords )
      type(WB_Subdomain), intent(in) :: sd
      integer(MP), dimension(:), allocatable, intent(inout) :: block_coords

      block_coords = sd%block_coords
   end subroutine wb_subdomain_block_coords_vector

   function wb_subdomain_block_number( sd ) result( ib )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: ib

      ib = sd%ib
   end function wb_subdomain_block_number

   function wb_subdomain_block_number_mp( sd ) result( ib )
      type(WB_Subdomain), intent(in) :: sd
      integer(MP) :: ib

      ib = int(wb_subdomain_block_number(sd),MP)
   end function wb_subdomain_block_number_mp

   function wb_subdomain_block_rank( sd ) result( block_rank )
      type(WB_Subdomain), intent(in) :: sd
      integer(MP) :: block_rank

      block_rank = sd%block_rank
   end function

   subroutine wb_subdomain_case_name( sd, case_name )
      type(WB_Subdomain), intent(in) :: sd
      character(len=:), allocatable, intent(inout) :: case_name

      case_name = sd%case_name
   end subroutine wb_subdomain_case_name

   subroutine wb_subdomain_construct_namelist( sd, filename )
      type(WB_Subdomain), intent(inout) :: sd
      character(len=STRING_LENGTH), intent(in) :: filename
      character(len=STRING_LENGTH) :: case_name
      integer(SP) :: ib, nb, nc, n_dim, ng
      integer(MP) :: ierr, world_rank, world_size
      type(WB_Block), dimension(:), allocatable :: blocks

      call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr )
      call mpi_comm_size( MPI_COMM_WORLD, world_size, ierr )
      call read_general_namelist( filename, case_name, nb, nc, n_dim, ng )
      if ( world_rank .eq. WORLD_MASTER ) then
         call check_general_variables( nb, nc, n_dim, ng, world_size )
      end if
      call mpi_barrier( MPI_COMM_WORLD, ierr )
      call read_block_namelists( filename, nb, n_dim, blocks )
      if ( world_rank .eq. WORLD_MASTER ) then
         call check_block_dimension_arrays( blocks, nb, n_dim, ng )
         call check_block_neighbors( blocks, nb, n_dim )
         call check_block_world_size( blocks, nb )
      end if
      call mpi_barrier( MPI_COMM_WORLD, ierr )
      call wb_subdomain_construct_variables( sd, nb, nc, n_dim, ng, blocks, &
         case_name )

      do ib = 1_SP, nb
         call wb_block_destroy( blocks(ib) )
      end do
      deallocate( blocks )
   end subroutine wb_subdomain_construct_namelist

   subroutine wb_subdomain_construct_variables( sd, nb, nc, n_dim, ng, &
      blocks, case_name )
      type(WB_Subdomain), intent(inout) :: sd
      character(len=STRING_LENGTH), optional, intent(in) :: case_name
      integer(SP), intent(in) :: nb, nc, n_dim, ng
      integer(MP) :: ierr, world_rank
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      type(WB_Process), dimension(:), allocatable :: processes

      if ( present(case_name) ) then
         sd%case_name = trim(case_name)
      else
         sd%case_name = trim(DEFAULT_CASE_NAME)
      end if

      sd%nb    = nb
      sd%nc    = nc
      sd%n_dim = n_dim
      sd%ng    = ng

      call mpi_comm_rank( MPI_COMM_WORLD, sd%world_rank, ierr )
      call mpi_comm_size( MPI_COMM_WORLD, sd%world_size, ierr )

      allocate( sd%nx(sd%n_dim) )
      allocate( sd%block_coords(sd%n_dim) )
      allocate( sd%neighbors(sd%n_dim,NUMBER_OF_DIRECTIONS) )

      call decompose_domain( sd, blocks, processes )
      call check_points( sd )
      call check_block_total_points( sd )
      call identify_process_neighbors( sd, blocks, processes )
      do world_rank = 0_MP, sd%world_size-1_MP
         call wb_process_destroy( processes(world_rank) )
      end do
      deallocate( processes )

      call calculate_number_of_variables( sd )
   end subroutine wb_subdomain_construct_variables

   subroutine wb_subdomain_destroy( sd )
      integer(MP) :: ierr
      type(WB_Subdomain), intent(inout) :: sd

      deallocate( sd%case_name )
      deallocate( sd%nx )
      deallocate( sd%block_coords )
      deallocate( sd%neighbors )
      call mpi_comm_free( sd%comm_block, ierr )
      call wb_block_destroy( sd%local_block )
   end subroutine wb_subdomain_destroy

   function wb_subdomain_dimensions( sd ) result( n_dim )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: n_dim

      n_dim = sd%n_dim
   end function wb_subdomain_dimensions

   function wb_subdomain_dimensions_mp( sd ) result( n_dim )
      type(WB_Subdomain), intent(in) :: sd
      integer(MP) :: n_dim

      n_dim = int(wb_subdomain_dimensions(sd),MP)
   end function wb_subdomain_dimensions_mp

   function wb_subdomain_faces( sd ) result( n_faces )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: n_faces

      n_faces = sd%n_dim * NUMBER_OF_DIRECTIONS
   end function wb_subdomain_faces

   function wb_subdomain_ghost_points( sd ) result( ng )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: ng

      ng = sd%ng
   end function wb_subdomain_ghost_points

   function wb_subdomain_is_block_master( sd ) result( is_block_master )
      type(WB_Subdomain), intent(in) :: sd
      logical :: is_block_master

      is_block_master = sd%block_rank .eq. BLOCK_MASTER
   end function wb_subdomain_is_block_master

   function wb_subdomain_is_world_master( sd ) result( is_world_master )
      type(WB_Subdomain), intent(in) :: sd
      logical :: is_world_master

      is_world_master = sd%world_rank .eq. WORLD_MASTER
   end function wb_subdomain_is_world_master

   subroutine wb_subdomain_local_block( sd, local_block )
      type(WB_Subdomain), intent(in) :: sd
      type(WB_Block), intent(inout) :: local_block

      local_block = sd%local_block
   end subroutine wb_subdomain_local_block

   function wb_subdomain_neighbor( sd, i_dim, i_dir ) result( neighbor )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP), intent(in) :: i_dim, i_dir
      integer(MP) :: neighbor

      neighbor = sd%neighbors(i_dim,i_dir)
   end function wb_subdomain_neighbor

   function wb_subdomain_points( sd, i_dim ) result( points )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP), intent(in) :: i_dim
      integer(SP) :: points

      if ( i_dim .lt. MIN_NUMBER_OF_DIMENSIONS .or. &
           i_dim .gt. MAX_NUMBER_OF_DIMENSIONS ) then
         points = 0_SP
      else if ( i_dim .gt. sd%n_dim ) then
         points = 1_SP
      else
         points = sd%nx(i_dim)
      end if
   end function wb_subdomain_points

   function wb_subdomain_total_blocks( sd ) result( total_blocks )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: total_blocks

      total_blocks = sd%nb
   end function wb_subdomain_total_blocks

   function wb_subdomain_components( sd ) result( components )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: components

      components = sd%nc
   end function wb_subdomain_components

   function wb_subdomain_total_points( sd ) result( points_in_state )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: points_in_state

      points_in_state = product(sd%nx)
   end function wb_subdomain_total_points

   function wb_subdomain_variables( sd ) result( nv )
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: nv

      nv = sd%nv
   end function wb_subdomain_variables

   function wb_subdomain_world_rank( sd ) result( world_rank )
      type(WB_Subdomain), intent(in) :: sd
      integer(MP) :: world_rank

      world_rank = sd%world_rank
   end function wb_subdomain_world_rank

   function wb_subdomain_world_size( sd ) result( world_size )
      type(WB_Subdomain), intent(in) :: sd
      integer(MP) :: world_size

      world_size = sd%world_size
   end function wb_subdomain_world_size

   subroutine write_block_information( f, sd )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: sd
      integer(SP) :: ib, i_dim
      integer(MP) :: ierr
      character(len=STRING_LENGTH) :: label
      type(WB_Block) :: local_block

      if ( wb_subdomain_is_world_master(sd) ) then
         call write_log_heading( f, "Block information", level=2_SP )

         call write_table_entry( f, "`ib`", IB_COLUMN_WIDTH )
         call write_table_entry( f, "`block_size`", SIZE_COLUMN_WIDTH )
         do i_dim = 1_SP, num_dimensions(sd)
            write (label,"(A, I1, A)") "`np(", i_dim, ")`"
            call write_table_entry( f, label, NP_COLUMN_WIDTH )
         end do
         call write_table_entry( f, "points", POINTS_COLUMN_WIDTH )
         do i_dim = 1_SP, num_dimensions(sd)
            write (label,"(A, I1, A)") "`nx(", i_dim, ")`"
            call write_table_entry( f, label, NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. num_dimensions(sd)) )
         end do

         call write_table_rule_entry( f, IB_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_table_rule_entry( f, SIZE_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, num_dimensions(sd)
            call write_table_rule_entry( f, NP_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED )
         end do
         call write_table_rule_entry( f, POINTS_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, num_dimensions(sd)
            call write_table_rule_entry( f, NX_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED, end_row=(i_dim .eq. num_dimensions(sd)) )
         end do
      end if

      do ib = 1_SP, num_blocks(sd)
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( ib .eq. wb_subdomain_block_number(sd) .and. &
            wb_subdomain_is_block_master(sd) ) then
            call wb_subdomain_local_block( sd, local_block )
            call write_table_entry( f, wb_subdomain_block_number(sd), &
               IB_COLUMN_WIDTH )
            call write_table_entry( f, int(wb_block_size(local_block),SP), &
               SIZE_COLUMN_WIDTH )
            do i_dim = 1_SP, num_dimensions(sd)
               call write_table_entry( f, &
               int(wb_block_processes(local_block,i_dim),SP), &
                  NP_COLUMN_WIDTH )
            end do
            call write_table_entry( f, total_points(local_block), &
               POINTS_COLUMN_WIDTH )
            do i_dim = 1_SP, num_dimensions(sd)
               call write_table_entry( f, num_points(local_block,i_dim), &
                  NX_COLUMN_WIDTH, end_row=(i_dim .eq. num_dimensions(sd)) )
            end do
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( wb_subdomain_is_world_master(sd) ) then
         call write_blank_line( f )
      end if
   end subroutine write_block_information

   subroutine write_environment_information( f, sd )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: sd
      integer(MP) :: ierr, mpi_major_version_number, &
         mpi_minor_version_number, version_length
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: lib_version

      if ( wb_subdomain_is_world_master(sd) ) then
         call write_log_heading( f, PROGRAM_NAME )

         call write_log_heading( f, "Environment parameters", level=2_SP )
         call write_table_entry( f, "Question", QUESTION_COLUMN_WIDTH )
         call write_table_entry( f, "Answer", ANSWER_COLUMN_WIDTH, &
            end_row=.true. )
         call write_table_rule_entry( f, QUESTION_COLUMN_WIDTH )
         call write_table_rule_entry( f, ANSWER_COLUMN_WIDTH, end_row=.true. )
         call write_table_entry( f, "Are MPI subarrays supported?", &
            QUESTION_COLUMN_WIDTH )
         call write_table_entry( f, MPI_SUBARRAYS_SUPPORTED, &
            ANSWER_COLUMN_WIDTH, end_row=.true. )
         call write_table_entry( f, "Is it big-endian?", &
            QUESTION_COLUMN_WIDTH )
         call write_table_entry( f, ARCH_IS_BIG_ENDIAN, ANSWER_COLUMN_WIDTH, &
            end_row=.true. )
         call write_blank_line( f )

         call write_log_heading( f, "Data precision", level=2_SP )
         call write_table_entry( f, "Data type", DATA_TYPE_COLUMN_WIDTH )
         call write_table_entry( f, "Variable", VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, "Fortran", VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, "MPI", VARIABLE_COLUMN_WIDTH, &
            end_row=.true. )
         call write_table_rule_entry( f, DATA_TYPE_COLUMN_WIDTH )
         call write_table_rule_entry( f, VARIABLE_COLUMN_WIDTH )
         call write_table_rule_entry( f, VARIABLE_COLUMN_WIDTH )
         call write_table_rule_entry( f, VARIABLE_COLUMN_WIDTH, &
            end_row=.true. )
         call write_table_entry( f, "Floating point", DATA_TYPE_COLUMN_WIDTH )
         call write_table_entry( f, "`FP`", VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, identify_real_precision(FP), &
            VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, identify_real_precision(MPI_FP), &
            VARIABLE_COLUMN_WIDTH, end_row=.true. )
         call write_table_entry( f, "Integer", DATA_TYPE_COLUMN_WIDTH )
         call write_table_entry( f, "`SP`", VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, identify_integer_precision(SP), &
            VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, identify_integer_precision(MPI_SP), &
            VARIABLE_COLUMN_WIDTH, end_row=.true. )
         call write_table_entry( f, "Integer (MPI only)", &
            DATA_TYPE_COLUMN_WIDTH )
         call write_table_entry( f, "`MP`", VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, identify_integer_precision(MP), &
            VARIABLE_COLUMN_WIDTH )
         call write_table_entry( f, identify_integer_precision(MPI_MP), &
            VARIABLE_COLUMN_WIDTH, end_row=.true. )
         call write_blank_line( f )

         write (f,"(A, A)") "- Program version: ", VERSION
         call write_blank_line( f )

         call mpi_get_version( mpi_major_version_number, &
            mpi_minor_version_number, ierr )
         write (f,"(A, I1, A, I1)") "- MPI version: ", &
            mpi_major_version_number, ".", mpi_minor_version_number
         call write_blank_line( f )

         call mpi_get_library_version( lib_version, version_length, ierr )
         write (f,"(A, A)") "- MPI library version: ", trim(lib_version)
         call write_blank_line( f )

         write (f,"(A, A, A, A, A)") "- Compiled using ", compiler_version(), &
            " using the following options: `", compiler_options(), "`"
         call write_blank_line( f )
      end if
   end subroutine write_environment_information

   subroutine write_subdomain_information( f, sd )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: sd
      integer(MP) :: ierr, processor_length, world_rank
      integer(SP) :: i_dim
      character(len=STRING_LENGTH) :: label
      character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

      if ( wb_subdomain_is_world_master(sd) ) then
         call write_log_heading( f, "Subdomain information", level=2_SP )

         call write_table_entry( f, "`world_rank`", RANK_COLUMN_WIDTH )
         call write_table_entry( f, "hostname", HOSTNAME_COLUMN_WIDTH )
         call write_table_entry( f, "`ib`",           IB_COLUMN_WIDTH )
         call write_table_entry( f, "`block_rank`", RANK_COLUMN_WIDTH )
         do i_dim = 1_SP, num_dimensions(sd)
            write (label,"(A, I1, A)") "(", i_dim, ")"
            call write_table_entry( f, label, COORDS_COLUMN_WIDTH )
         end do
         call write_table_entry( f, "points", POINTS_COLUMN_WIDTH )
         do i_dim = 1_SP, num_dimensions(sd)
            write (label,"(A, I1, A)") "`nx(", i_dim, ")`"
            call write_table_entry( f, label, NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. num_dimensions(sd)) )
         end do

         call write_table_rule_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_table_rule_entry( f, HOSTNAME_COLUMN_WIDTH, &
            alignment=LEFT_ALIGNED  )
         call write_table_rule_entry( f, IB_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_table_rule_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, num_dimensions(sd)
            call write_table_rule_entry( f, COORDS_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED )
         end do
         call write_table_rule_entry( f, POINTS_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, num_dimensions(sd)
            call write_table_rule_entry( f, NX_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED, end_row=(i_dim .eq. num_dimensions(sd)) )
         end do
      end if

      call mpi_get_processor_name( processor_name, processor_length, ierr )
      do world_rank = 0_MP, wb_subdomain_world_size(sd)-1_MP
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( wb_subdomain_world_rank(sd) .eq. world_rank ) then
            call write_table_entry( f, int(wb_subdomain_world_rank(sd),SP), &
               RANK_COLUMN_WIDTH )
            call write_table_entry(  f, processor_name, &
               HOSTNAME_COLUMN_WIDTH )
            call write_table_entry( f, wb_subdomain_block_number(sd), &
               IB_COLUMN_WIDTH )
            call write_table_entry( f, int(wb_subdomain_block_rank(sd),SP), &
               RANK_COLUMN_WIDTH )
            do i_dim = 1_SP, num_dimensions(sd)
               call write_table_entry( f, &
               int(wb_subdomain_block_coord(sd,i_dim),SP), &
               COORDS_COLUMN_WIDTH )
            end do
            call write_table_entry( f, total_points(sd), &
               POINTS_COLUMN_WIDTH )
            do i_dim = 1_SP, num_dimensions(sd)
               call write_table_entry( f, num_points(sd,i_dim), &
                  NX_COLUMN_WIDTH, end_row=(i_dim .eq. num_dimensions(sd)) )
            end do
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( wb_subdomain_is_world_master(sd) ) then
         call write_blank_line( f )
      end if
   end subroutine write_subdomain_information

   subroutine write_subdomain_neighbors( f, sd )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: sd
      integer(MP) :: ierr, world_rank, sd_neighbor
      integer(SP) :: i_dim, i_dir, j_dir, face_count
      integer(SP), dimension(NUMBER_OF_DIRECTIONS) :: dirs
      character(len=STRING_LENGTH) :: label

      face_count = 0_SP
      dirs = (/ LOWER_DIRECTION, UPPER_DIRECTION /)

      if ( wb_subdomain_is_world_master(sd) ) then
         call write_log_heading( f, "Subdomain neighbors", level=2_SP )

         call write_table_entry( f, "`world_rank`", RANK_COLUMN_WIDTH )
         do i_dim = 1_SP, num_dimensions(sd)
            do i_dir = 1_SP, NUMBER_OF_DIRECTIONS
               j_dir = dirs(i_dir)
               face_count = face_count + 1_SP
               if ( j_dir .eq. LOWER_DIRECTION ) then
                  write (label,"(I1, A)") i_dim, "L"
               else
                  write (label,"(I1, A)") i_dim, "U"
               end if
               call write_table_entry( f, label, RANK_COLUMN_WIDTH, &
                  end_row=( face_count .eq. wb_subdomain_faces(sd) ) )
            end do
         end do

         face_count = 0_SP
         call write_table_rule_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, num_dimensions(sd)
            do i_dir = 1_SP, NUMBER_OF_DIRECTIONS
               face_count = face_count + 1_SP
               call write_table_rule_entry( f, RANK_COLUMN_WIDTH, &
                  alignment=RIGHT_ALIGNED, &
                  end_row=( face_count .eq. wb_subdomain_faces(sd) ) )
            end do
         end do
      end if

      do world_rank = 0_MP, wb_subdomain_world_size(sd)-1_MP
         face_count = 0_SP
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( wb_subdomain_world_rank(sd) .eq. world_rank ) then
            call write_table_entry( f, int(wb_subdomain_world_rank(sd),SP), &
               RANK_COLUMN_WIDTH )
            do i_dim = 1_SP, num_dimensions(sd)
               do i_dir = 1_SP, NUMBER_OF_DIRECTIONS
                  j_dir = dirs(i_dir)
                  face_count = face_count + 1_SP
                  sd_neighbor = neighbor( sd, i_dim, j_dir )
                  if ( sd_neighbor .eq. MPI_PROC_NULL ) then
                     call write_table_entry( f, "-", &
                        RANK_COLUMN_WIDTH, &
                        end_row=( face_count .eq. wb_subdomain_faces(sd) ) )
                  else
                     call write_table_entry( f, int(sd_neighbor,SP), &
                        RANK_COLUMN_WIDTH, &
                        end_row=( face_count .eq. wb_subdomain_faces(sd) ) )
                  end if
               end do
            end do
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( wb_subdomain_is_world_master(sd) ) then
         call write_blank_line( f )
      end if
   end subroutine write_subdomain_neighbors

   subroutine write_scalar_variables( f, sd )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: sd
      character(len=:), allocatable :: case_name

      call wb_subdomain_case_name( sd, case_name )

      if ( wb_subdomain_is_world_master(sd) ) then
         call write_log_heading( f, "State scalar variables", level=2_SP )
         call write_table_entry( f, "Property", PROPERTY_COLUMN_WIDTH )
         call write_table_entry( f, "Value", VALUE_COLUMN_WIDTH, &
            end_row=.true. )
         call write_table_rule_entry( f, PROPERTY_COLUMN_WIDTH, &
            alignment=LEFT_ALIGNED )
         call write_table_rule_entry( f, VALUE_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED, end_row=.true. )
         call write_table_entry( f, "Case name", PROPERTY_COLUMN_WIDTH )
         call write_table_entry( f, case_name, VALUE_COLUMN_WIDTH, &
            end_row=.true. )
         call write_table_entry( f, "Number of blocks", PROPERTY_COLUMN_WIDTH )
         call write_table_entry( f, num_blocks(sd), &
            VALUE_COLUMN_WIDTH, end_row=.true. )
         call write_table_entry( f, "Number of components", PROPERTY_COLUMN_WIDTH )
         call write_table_entry( f, num_components(sd), &
            VALUE_COLUMN_WIDTH, end_row=.true. )
         call write_table_entry( f, "Number of dimensions", &
            PROPERTY_COLUMN_WIDTH )
         call write_table_entry( f, num_dimensions(sd), &
            VALUE_COLUMN_WIDTH, end_row=.true. )
         call write_table_entry( f, "Number of ghost points", &
            PROPERTY_COLUMN_WIDTH )
         call write_table_entry( f, num_ghost_points(sd), &
            VALUE_COLUMN_WIDTH, end_row=.true. )
         call write_table_entry( f, "Number of variables", &
            PROPERTY_COLUMN_WIDTH )
         call write_table_entry( f, num_variables(sd), &
            VALUE_COLUMN_WIDTH, end_row=.true. )
         call write_blank_line( f )
      end if

      deallocate( case_name )
   end subroutine write_scalar_variables
end module wb_base
