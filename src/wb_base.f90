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

   integer(SP), public, parameter :: NO_BLOCK_NEIGHBOR = 0_SP

   integer(SP), public, parameter :: MIN_N_DIM = 1_SP
   integer(SP), public, parameter :: MAX_N_DIM = 3_SP

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

   interface total_points
      module procedure wb_block_total_points, wb_subdomain_total_points
   end interface total_points

   interface wb_subdomain_construct
      module procedure wb_subdomain_construct_namelist, &
         wb_subdomain_construct_variables
   end interface wb_subdomain_construct
contains
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
            if ( wb_block_points( blocks(ib), i_dim ) .lt. ng ) then
               call wb_abort( "number of points in direction N1 of block &
                              &N2 is less than number of ghost points N3", &
                              EXIT_DATAERR, (/ i_dim, ib, ng /) )
            end if
            do i_dir = 1_SP, N_DIR
               if ( wb_block_neighbor( blocks(ib), i_dim, i_dir ) .lt. &
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
            neighbor_l = wb_block_neighbor( blocks(ib), i_dim, LOWER_DIR )
            if ( neighbor_l .ne. NO_BLOCK_NEIGHBOR ) then
               neighbor_u = wb_block_neighbor( blocks(neighbor_l), i_dim, &
                  UPPER_DIR )
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
                        wb_block_points( blocks(ib),         j_dim ) .ne. &
                        wb_block_points( blocks(neighbor_l), j_dim ) ) then
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

   subroutine check_block_total_points( s, blocks )
      integer(SP) :: points_in_block, points_in_processes
      integer(MP) :: ierr
      type(WB_Subdomain), intent(in) :: s
      type(WB_Block), dimension(:), allocatable :: blocks

      points_in_block = total_points(blocks(s%ib))
      points_in_processes = 0_SP
      call mpi_reduce( total_points(s), points_in_processes, &
         int(s%n_dim,MP), MPI_SP, MPI_SUM, BLOCK_MASTER, s%comm_block, ierr )
      if ( s%block_rank .eq. BLOCK_MASTER .and. &
         points_in_block .ne. points_in_processes ) then
         call wb_abort( "total points in block N1 (N2) does not match sum of &
                        &points in individual processes (N3)", &
            EXIT_FAILURE, &
            (/ s%ib, points_in_block, points_in_processes /) )
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

   subroutine check_general_variables( nb, n_dim, ng, world_size )
      integer(SP), intent(in) :: nb, n_dim, ng
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
      if ( ng .lt. 1_SP ) then
         call wb_abort( "number of ghost points is less than 1", &
            EXIT_DATAERR )
      end if
      if ( n_dim .lt. MIN_N_DIM .or. n_dim .gt. MAX_N_DIM ) then
         call wb_abort( &
            "number of dimensions must be in interval [N1, N2]", &
            EXIT_DATAERR, (/ MIN_N_DIM, MAX_N_DIM /) )
      end if
   end subroutine check_general_variables

   subroutine decompose_domain( s, blocks, processes )
      integer(MP) :: assigned_processes, ierr, world_rank
      integer(SP) :: ib, i_dim
      type(MPI_Comm) :: comm_split
      type(WB_Subdomain), intent(inout) :: s
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      type(WB_Process), dimension(:), allocatable, intent(out) :: processes
      integer(SP), dimension(:), allocatable :: nx
      integer(MP), dimension(:), allocatable :: np
      logical, dimension(:), allocatable :: periods

      ib = 1_SP
      assigned_processes = 0_MP
      allocate( processes(0_MP:s%world_size-1_MP) )
      do world_rank = 0_MP, s%world_size-1_MP
         call wb_process_construct( processes(world_rank), s%n_dim, ib )
         assigned_processes = assigned_processes + 1_MP
         if ( assigned_processes .eq. wb_block_size( blocks(ib) ) ) then
            assigned_processes = 0_MP
            ib = ib + 1_SP
         end if
      end do

      s%ib = processes(s%world_rank)%ib
      allocate( periods(s%n_dim), nx(s%n_dim), np(s%n_dim) )
      call wb_block_periods_vector(   blocks(s%ib), periods )
      call wb_block_points_vector(    blocks(s%ib), nx      )
      call wb_block_processes_vector( blocks(s%ib), np      )
      call mpi_comm_split( MPI_COMM_WORLD, int(s%ib,MP), 0_MP, comm_split, &
         ierr )
      call mpi_cart_create( comm_split, int(s%n_dim,MP), np, periods, &
         wb_block_reorder( blocks(s%ib) ), s%comm_block, ierr )
      call mpi_comm_free( comm_split, ierr )
      call mpi_comm_rank( s%comm_block, s%block_rank, ierr )
      call mpi_comm_size( s%comm_block, s%block_size, ierr )
      call mpi_cart_coords( s%comm_block, s%block_rank, int(s%n_dim,MP), &
         s%block_coords, ierr )

      s%nx = nx / int(np,SP)
      do i_dim = 1_SP, s%n_dim
         if ( s%block_coords(i_dim) .eq. &
            wb_block_processes( blocks(s%ib), i_dim ) - 1_MP ) then
            s%nx(i_dim) = s%nx(i_dim) + modulo( &
               wb_block_points( blocks(s%ib), i_dim ), &
               int( wb_block_processes( blocks(s%ib), i_dim ), SP ) )
         end if
      end do

      processes(s%world_rank)%block_rank   = s%block_rank
      processes(s%world_rank)%block_coords = s%block_coords
      processes(s%world_rank)%nx           = s%nx
      do world_rank = 0_MP, s%world_size-1_MP
         call mpi_bcast( processes(world_rank)%block_rank, 1_MP, &
            MPI_INTEGER, world_rank, MPI_COMM_WORLD, ierr )
         call mpi_bcast( processes(world_rank)%block_coords, int(s%n_dim,MP), &
            MPI_INTEGER, world_rank, MPI_COMM_WORLD, ierr )
         call mpi_bcast( processes(world_rank)%nx, int(s%n_dim,MP), &
            MPI_SP, world_rank, MPI_COMM_WORLD, ierr )
      end do

      deallocate( periods, nx, np )
   end subroutine decompose_domain

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

   subroutine identify_process_neighbors( s, blocks, processes )
      integer(MP) :: ierr, world_rank
      integer(SP) :: block_neighbor, i_dir, i_dim
      integer(MP), dimension(N_DIR) :: block_ranks
      integer(MP), dimension(:), allocatable :: block_coords
      type(WB_Subdomain), intent(inout) :: s
      type(WB_Block), dimension(:), allocatable, intent(in) :: blocks
      type(WB_Process), dimension(:), allocatable, intent(in) :: processes

      allocate( block_coords(s%n_dim) )
      do i_dim = 1_SP, s%n_dim
         call mpi_cart_shift( s%comm_block, int(i_dim,MP)-1_MP, 1_MP, &
            block_ranks(LOWER_DIR), block_ranks(UPPER_DIR), ierr )
         do i_dir = 1_SP, N_DIR
            if ( block_ranks(i_dir) .ne. MPI_PROC_NULL ) then
               do world_rank = 0_MP, s%world_size-1_MP
                  if ( processes(world_rank)%ib .eq. s%ib .and. &
                     processes(world_rank)%block_rank .eq. &
                     block_ranks(i_dir) ) then
                     exit
                  end if
               end do
               s%neighbors(i_dim,i_dir) = world_rank
            else
               block_neighbor = blocks(s%ib)%neighbors(i_dim,i_dir)
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
                        blocks(block_neighbor)%np(i_dim)-1_MP
                  else
                     block_coords(i_dim) = 0_MP
                  end if
                  do world_rank = 0_MP, s%world_size-1_MP
                     if ( processes(world_rank)%ib .eq. block_neighbor &
                        .and. all( processes(world_rank)%block_coords &
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

   subroutine print_initial_information( s )
      type(WB_Subdomain), intent(in) :: s

      call write_environment_information( output_unit, s )
      call write_scalar_variables(        output_unit, s )
      call write_block_information(       output_unit, s )
      call write_process_information(     output_unit, s )
      call write_process_neighbors(       output_unit, s )
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

      ib             = DEFAULT_IB
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

   subroutine read_general_namelist( filename, case_name, nb, n_dim, ng )
      character(len=STRING_LENGTH), intent(in) :: filename
      character(len=STRING_LENGTH), intent(out) :: case_name
      integer(SP), intent(out) :: nb, n_dim, ng
      integer :: file_unit
      integer(MP) :: ierr, world_rank
      namelist /general/ case_name, nb, ng, n_dim

      case_name = DEFAULT_CASE_NAME
      nb        = DEFAULT_NB
      n_dim     = DEFAULT_N_DIM
      ng        = DEFAULT_NG

      call mpi_comm_rank( MPI_COMM_WORLD, world_rank, ierr )
      if ( world_rank .eq. WORLD_MASTER ) then
         open( newunit=file_unit, file=filename, form="formatted", &
            action="read" )
         read( unit=file_unit, nml=general )
         close( unit=file_unit )

      end if

      call mpi_bcast( case_name, int(STRING_LENGTH,MP), MPI_CHARACTER, &
         WORLD_MASTER, MPI_COMM_WORLD, ierr )
      call mpi_bcast( nb, 1_MP, MPI_SP, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
      call mpi_bcast( n_dim, 1_MP, MPI_SP, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
      call mpi_bcast( ng, 1_MP, MPI_SP, WORLD_MASTER, &
         MPI_COMM_WORLD, ierr )
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
      allocate( blk%neighbors(n_dim,N_DIR) )
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
         blk%neighbors(:,LOWER_DIR) = neighbors_l
      end if
      if ( present(neighbors_u) ) then
         blk%neighbors(:,UPPER_DIR) = neighbors_u
      end if
      blk%reorder = DEFAULT_REORDER

      if ( present(ib) ) then
         do i_dim = 1_SP, n_dim
            if ( blk%neighbors(i_dim,LOWER_DIR) .eq. ib .and. &
                 blk%neighbors(i_dim,UPPER_DIR) .eq. ib ) then
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

   subroutine wb_block_periods_vector( blk, periods )
      type(WB_Block), intent(in) :: blk
      logical, dimension(:), allocatable, intent(out) :: periods

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
      integer(SP), dimension(:), allocatable, intent(out) :: nx

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
      integer(MP), dimension(:), allocatable, intent(out) :: np

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

   subroutine wb_process_construct( process, n_dim, ib )
      type(WB_Process), intent(inout) :: process
      integer(SP), intent(in) :: n_dim, ib

      allocate( process%block_coords(n_dim) )
      allocate( process%nx(n_dim) )

      process%ib = ib
   end subroutine wb_process_construct

   subroutine wb_process_destroy( process )
      type(WB_Process), intent(inout) :: process

      deallocate( process%block_coords )
      deallocate( process%nx )
   end subroutine wb_process_destroy

   subroutine wb_subdomain_construct_namelist( s, filename )
      type(WB_Subdomain), intent(inout) :: s
      character(len=STRING_LENGTH), intent(in) :: filename
      character(len=STRING_LENGTH) :: case_name
      integer(SP) :: ib, nb, n_dim, ng
      integer(MP) :: world_rank
      type(WB_Block), dimension(:), allocatable :: blocks
      type(WB_Process), dimension(:), allocatable :: processes

      call read_general_namelist( filename, case_name, nb, n_dim, ng )
      call wb_subdomain_construct_variables( s, nb, n_dim, ng, case_name )
      call read_block_namelists( filename, nb, n_dim, blocks )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         call check_block_dimension_arrays( blocks, nb, n_dim, ng )
         call check_block_neighbors( blocks, nb, n_dim )
         call check_block_world_size( blocks, nb )
      end if
      call decompose_domain( s, blocks, processes )
      call check_block_total_points( s, blocks )
      call identify_process_neighbors( s, blocks, processes )
      s%local_block = blocks(s%ib)

      do ib = 1_SP, nb
         call wb_block_destroy( blocks(ib) )
      end do
      deallocate( blocks )

      do world_rank = 0_MP, s%world_size-1_MP
         call wb_process_destroy( processes(world_rank) )
      end do
      deallocate( processes )
   end subroutine wb_subdomain_construct_namelist

   subroutine wb_subdomain_construct_variables( s, nb, n_dim, ng, case_name )
      type(WB_Subdomain), intent(inout) :: s
      character(len=STRING_LENGTH), optional, intent(in) :: case_name
      integer(SP), optional, intent(in) :: nb, n_dim, ng
      integer(MP) :: ierr

      if ( present(case_name) ) then
         s%case_name = trim(case_name)
      else
         s%case_name = trim(DEFAULT_CASE_NAME)
      end if

      if ( present(nb) ) then
         s%nb = nb
      else
         s%nb = DEFAULT_NB
      end if

      if ( present(n_dim) ) then
         s%n_dim = n_dim
      else
         s%n_dim = DEFAULT_N_DIM
      end if

      if ( present(ng) ) then
         s%ng = ng
      else
         s%ng = DEFAULT_NG
      end if

      call mpi_comm_rank( MPI_COMM_WORLD, s%world_rank, ierr )
      call mpi_comm_size( MPI_COMM_WORLD, s%world_size, ierr )
      call check_general_variables( s%nb, s%n_dim, s%ng, s%world_size )

      allocate( s%nx(s%n_dim) )
      allocate( s%block_coords(s%n_dim) )
      allocate( s%neighbors(s%n_dim,N_DIR) )
      call wb_block_construct( s%local_block, s%n_dim )
   end subroutine wb_subdomain_construct_variables

   subroutine wb_subdomain_destroy( s )
      integer(MP) :: ierr
      type(WB_Subdomain), intent(inout) :: s

      deallocate( s%case_name )
      deallocate( s%nx )
      deallocate( s%block_coords )
      deallocate( s%neighbors )
      call mpi_comm_free( s%comm_block, ierr )
      call wb_block_destroy( s%local_block )
   end subroutine wb_subdomain_destroy

   function wb_subdomain_total_points( s ) result( points_in_state )
      type(WB_Subdomain), intent(in) :: s
      integer(SP) :: points_in_state

      points_in_state = product(s%nx)
   end function wb_subdomain_total_points

   subroutine write_block_information( f, s )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: s
      integer(SP) :: ib, i_dim
      integer(MP) :: ierr
      character(len=STRING_LENGTH) :: label

      if ( s%world_rank .eq. WORLD_MASTER ) then
         call write_log_heading( f, "Block information", level=2_SP )

         call write_table_entry( f, "`ib`", IB_COLUMN_WIDTH )
         call write_table_entry( f, "`block_size`", SIZE_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
            write (label,"(A, I1, A)") "`np(", i_dim, ")`"
            call write_table_entry( f, label, NP_COLUMN_WIDTH )
         end do
         call write_table_entry( f, "points", POINTS_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
            write (label,"(A, I1, A)") "`nx(", i_dim, ")`"
            call write_table_entry( f, label, NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. s%n_dim) )
         end do

         call write_table_rule_entry( f, IB_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_table_rule_entry( f, SIZE_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            call write_table_rule_entry( f, NP_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED )
         end do
         call write_table_rule_entry( f, POINTS_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            call write_table_rule_entry( f, NX_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED, end_row=(i_dim .eq. s%n_dim) )
         end do
      end if

      do ib = 1_SP, s%nb
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( ib .eq. s%ib .and. s%block_rank .eq. BLOCK_MASTER ) then
            call write_table_entry( f, s%ib, IB_COLUMN_WIDTH )
            call write_table_entry( f, int(wb_block_size(s%local_block),SP), &
               SIZE_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               call write_table_entry( f, int(s%local_block%np(i_dim),SP), &
                  NP_COLUMN_WIDTH )
            end do
            call write_table_entry( f, total_points(s%local_block), &
               POINTS_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               call write_table_entry( f, s%local_block%nx(i_dim), &
                  NX_COLUMN_WIDTH, end_row=(i_dim .eq. s%n_dim) )
            end do
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         call write_blank_line( f )
      end if
   end subroutine write_block_information

   subroutine write_environment_information( f, s )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: s
      integer(MP) :: ierr, mpi_major_version_number, &
         mpi_minor_version_number, version_length
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: lib_version

      if ( s%world_rank .eq. WORLD_MASTER ) then
         call write_log_heading( f, PROGRAM_NAME )

         call write_log_heading( f, "Environment parameters", level=2_SP )
         call write_table_entry( f, "Question", 30_SP )
         call write_table_entry( f, "Answer", 10_SP, end_row=.true. )
         call write_table_rule_entry( f, 30_SP )
         call write_table_rule_entry( f, 10_SP, end_row=.true. )
         call write_table_entry( f, "Are MPI subarrays supported?", 30_SP )
         call write_table_entry( f, MPI_SUBARRAYS_SUPPORTED, 10_SP, &
            end_row=.true. )
         call write_table_entry( f, "Is it big-endian?", 30_SP )
         call write_table_entry( f, ARCH_IS_BIG_ENDIAN, 10_SP, &
            end_row=.true. )
         call write_blank_line( f )

         call write_log_heading( f, "Data precision", level=2_SP )
         call write_table_entry( f, "Data type", 20_SP )
         call write_table_entry( f, "Variable", 10_SP )
         call write_table_entry( f, "Fortran", 15_SP )
         call write_table_entry( f, "MPI", 15_SP, end_row=.true. )
         call write_table_rule_entry( f, 20_SP )
         call write_table_rule_entry( f, 10_SP )
         call write_table_rule_entry( f, 15_SP )
         call write_table_rule_entry( f, 15_SP, end_row=.true. )
         call write_table_entry( f, "Floating point", 20_SP )
         call write_table_entry( f, "`FP`", 10_SP )
         call write_table_entry( f, identify_real_precision(FP), 15_SP )
         call write_table_entry( f, identify_real_precision(MPI_FP), 15_SP, &
            end_row=.true. )
         call write_table_entry( f, "Integer", 20_SP )
         call write_table_entry( f, "`SP`", 10_SP )
         call write_table_entry( f, identify_integer_precision(SP), 15_SP )
         call write_table_entry( f, identify_integer_precision(MPI_SP), &
            15_SP, end_row=.true. )
         call write_table_entry( f, "Integer (MPI only)", 20_SP )
         call write_table_entry( f, "`MP`", 10_SP )
         call write_table_entry( f, identify_integer_precision(MP), 15_SP )
         call write_table_entry( f, identify_integer_precision(MPI_MP), &
            15_SP, end_row=.true. )
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

   subroutine write_process_information( f, s )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: s
      integer(MP) :: ierr, processor_length, world_rank
      integer(SP) :: i_dim
      character(len=STRING_LENGTH) :: label
      character(len=MPI_MAX_PROCESSOR_NAME) :: processor_name

      if ( s%world_rank .eq. WORLD_MASTER ) then
         call write_log_heading( f, "Process information", level=2_SP )

         call write_table_entry( f, "`world_rank`", RANK_COLUMN_WIDTH )
         call write_table_entry( f, "hostname", HOSTNAME_COLUMN_WIDTH )
         call write_table_entry( f, "`ib`",           IB_COLUMN_WIDTH )
         call write_table_entry( f, "`block_rank`", RANK_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
            write (label,"(A, I1, A)") "(", i_dim, ")"
            call write_table_entry( f, label, COORDS_COLUMN_WIDTH )
         end do
         call write_table_entry( f, "points", POINTS_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
            write (label,"(A, I1, A)") "`nx(", i_dim, ")`"
            call write_table_entry( f, label, NX_COLUMN_WIDTH, &
               end_row=(i_dim .eq. s%n_dim) )
         end do

         call write_table_rule_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_table_rule_entry( f, HOSTNAME_COLUMN_WIDTH, &
            alignment=LEFT_ALIGNED  )
         call write_table_rule_entry( f, IB_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         call write_table_rule_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            call write_table_rule_entry( f, COORDS_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED )
         end do
         call write_table_rule_entry( f, POINTS_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            call write_table_rule_entry( f, NX_COLUMN_WIDTH, &
               alignment=RIGHT_ALIGNED, end_row=(i_dim .eq. s%n_dim) )
         end do
      end if

      call mpi_get_processor_name( processor_name, processor_length, ierr )
      do world_rank = 0_MP, s%world_size-1_MP
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            call write_table_entry( f, int(s%world_rank,SP), &
               RANK_COLUMN_WIDTH )
            call write_table_entry(  f, processor_name, &
               HOSTNAME_COLUMN_WIDTH )
            call write_table_entry( f, s%ib, &
               IB_COLUMN_WIDTH )
            call write_table_entry( f, int(s%block_rank,SP), &
               RANK_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               call write_table_entry( f, &
               int(s%block_coords(i_dim),SP), COORDS_COLUMN_WIDTH )
            end do
            call write_table_entry( f, total_points(s), &
               POINTS_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               call write_table_entry( f, s%nx(i_dim), &
                  NX_COLUMN_WIDTH, end_row=(i_dim .eq. s%n_dim) )
            end do
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         call write_blank_line( f )
      end if
   end subroutine write_process_information

   subroutine write_process_neighbors( f, s )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: s
      integer(MP) :: ierr, world_rank, neighbor
      integer(SP) :: i_dim, i_dir, j_dir, face_count
      integer(SP), dimension(N_DIR) :: dirs
      character(len=STRING_LENGTH) :: label

      face_count = 0_SP
      dirs = (/ LOWER_DIR, UPPER_DIR /)

      if ( s%world_rank .eq. WORLD_MASTER ) then
         call write_log_heading( f, "Process neighbors", level=2_SP )

         call write_table_entry( f, "`world_rank`", RANK_COLUMN_WIDTH )
         do i_dim = 1_SP, s%n_dim
            do i_dir = 1_SP, N_DIR
               j_dir = dirs(i_dir)
               face_count = face_count + 1_SP
               if ( j_dir .eq. LOWER_DIR ) then
                  write (label,"(I1, A)") i_dim, "L"
               else
                  write (label,"(I1, A)") i_dim, "U"
               end if
               call write_table_entry( f, label, RANK_COLUMN_WIDTH, &
                  end_row=( face_count .eq. int(s%n_dim*N_DIR,MP) ) )
            end do
         end do

         face_count = 0_SP
         call write_table_rule_entry( f, RANK_COLUMN_WIDTH, &
            alignment=RIGHT_ALIGNED )
         do i_dim = 1_SP, s%n_dim
            do i_dir = 1_SP, N_DIR
               face_count = face_count + 1_SP
               call write_table_rule_entry( f, RANK_COLUMN_WIDTH, &
                  alignment=RIGHT_ALIGNED, &
                  end_row=( face_count .eq. int(s%n_dim*N_DIR,MP) ) )
            end do
         end do
      end if

      do world_rank = 0_MP, s%world_size-1_MP
         face_count = 0_SP
         call mpi_barrier( MPI_COMM_WORLD, ierr )
         if ( s%world_rank .eq. world_rank ) then
            call write_table_entry( f, int(s%world_rank,SP), &
               RANK_COLUMN_WIDTH )
            do i_dim = 1_SP, s%n_dim
               do i_dir = 1_SP, N_DIR
                  j_dir = dirs(i_dir)
                  face_count = face_count + 1_SP
                  neighbor = s%neighbors(i_dim,j_dir)
                  if ( neighbor .eq. MPI_PROC_NULL ) then
                     call write_table_entry( f, "-", &
                        RANK_COLUMN_WIDTH, &
                        end_row=( face_count .eq. int(s%n_dim*N_DIR,MP) ) )
                  else
                     call write_table_entry( f, int(neighbor,SP), &
                        RANK_COLUMN_WIDTH, &
                        end_row=( face_count .eq. int(s%n_dim*N_DIR,MP) ) )
                  end if
               end do
            end do
         end if
      end do

      call mpi_barrier( MPI_COMM_WORLD, ierr )
      if ( s%world_rank .eq. WORLD_MASTER ) then
         call write_blank_line( f )
      end if
   end subroutine write_process_neighbors

   subroutine write_scalar_variables( f, s )
      integer, intent(in) :: f
      type(WB_Subdomain), intent(in) :: s

      if ( s%world_rank .eq. WORLD_MASTER ) then
         call write_log_heading( f, "State scalar variables", level=2_SP )
         call write_table_entry( f, "Variable", 30_SP )
         call write_table_entry( f, "Value", 20_SP, end_row=.true. )
         call write_table_rule_entry( f, 30_SP, alignment=LEFT_ALIGNED )
         call write_table_rule_entry( f, 20_SP, alignment=RIGHT_ALIGNED, &
            end_row=.true. )
         call write_table_entry( f, "Case name", 30_SP )
         call write_table_entry( f, s%case_name, 20_SP, &
            end_row=.true. )
         call write_table_entry( f, "Number of dimensions", 30_SP )
         call write_table_entry( f, s%n_dim, 20_SP, &
            end_row=.true. )
         call write_table_entry( f, "Number of blocks", 30_SP )
         call write_table_entry( f, s%nb, 20_SP, &
            end_row=.true. )
         call write_table_entry( f, "Number of ghost points", 30_SP )
         call write_table_entry( f, s%ng, 20_SP, &
            end_row=.true. )
         call write_blank_line( f )
      end if
   end subroutine write_scalar_variables
end module wb_base
