!-----------------------------------------------------------------------
!
! MODULE: sweep module
!> @brief
!> This module performs a sweep of the grid in the new parallelization
!> scheme based on the OpenCL branch.
!
!-----------------------------------------------------------------------


MODULE omp_sweep_module

USE global_module, ONLY: i_knd, r_knd

USE geom_module, ONLY: nx, ny_gl, nz_gl

USE data_module, ONLY: ng

USE sn_module, ONLY: nang, noct

  IMPLICIT NONE

  PRIVATE

  TYPE cell
    INTEGER(i_knd) :: i, j, k
  END TYPE cell

  TYPE plane
    INTEGER(i_knd) :: num_cells
    INTEGER(i_knd) :: index
    TYPE(cell), ALLOCATABLE, DIMENSION(:) :: cells
  END TYPE plane

  ! Flux arrays
  REAL(r_knd), DIMENSION(:,:,:,:,:,:), POINTER :: flux_in, flux_out
  REAL(r_knd), DIMENSION(:,:,:,:), POINTER :: flux_i, flux_j, flux_k

  PUBLIC :: omp_sweep_alloc, omp_sweep_dealloc, omp_sweep

CONTAINS


  SUBROUTINE omp_sweep_alloc
    ALLOCATE ( flux_in(nang,ng,nx,ny_gl,nz_gl,noct) )
    ALLOCATE ( flux_out(nang,ng,nx,ny_gl,nz_gl,noct) )
    ALLOCATE ( flux_i(nang,ng,ny_gl,nz_gl) )
    ALLOCATE ( flux_j(nang,ng,nx,nz_gl) )
    ALLOCATE ( flux_k(nang,ng,nx,ny_gl) )
  END SUBROUTINE omp_sweep_alloc


  SUBROUTINE omp_sweep_dealloc
    DEALLOCATE ( flux_in, flux_out )
    DEALLOCATE ( flux_i, flux_j, flux_k )
  END SUBROUTINE omp_sweep_dealloc


  SUBROUTINE zero_edge_flux
    flux_i = 0.0
    flux_j = 0.0
    flux_k = 0.0
  END SUBROUTINE zero_edge_flux



  SUBROUTINE omp_sweep
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: o, p, c, g, a, nplanes
    INTEGER(i_knd) :: xhi, yhi, zhi, istep, jstep, kstep
    INTEGER(i_knd) :: i, j, k
    TYPE(plane), ALLOCATABLE, DIMENSION(:) :: planes
!_______________________________________________________________________
!
!   Set up the sweep order
!_______________________________________________________________________

    nplanes = nx + ny_gl + nz_gl - 2

    ALLOCATE ( planes(nplanes) )

    CALL compute_sweep_order ( planes, nplanes )

!_______________________________________________________________________
!
!   Sweep over the grid
!_______________________________________________________________________

    CALL zero_edge_flux

    ! Loop over octants
    octants: DO o = 1, noct

      ! Calculate the corner co-ordinates for this octant
      SELECT CASE (o)
        CASE (1)
          xhi = nx
          yhi = ny_gl
          zhi = nz_gl
        CASE (2)
          xhi = nx
          yhi = ny_gl
          zhi = 1
        CASE (3)
          xhi = nx
          yhi = 1
          zhi = nz_gl
        CASE (4)
          xhi = nx
          yhi = 1
          zhi = 1
        CASE (5)
         xhi = 1
         yhi = ny_gl
         zhi = nz_gl
        CASE (6)
          xhi = 1
          yhi = ny_gl
          zhi = 1
        CASE (7)
          xhi = 1
          yhi = 1
          zhi = nz_gl
        CASE (8)
          xhi = 1
          yhi = 1
          zhi = 1
      END SELECT

      ! Set the order to traverse each axis
      IF ( xhi == nx ) THEN
        istep = -1
      ELSE
        istep = 1
      END IF

      IF ( yhi == ny_gl ) THEN
        jstep = -1
      ELSE
        jstep = 1
      END IF

      IF ( zhi == nz_gl ) THEN
        kstep = -1
      ELSE
        kstep = 1
      END IF

      ! Loop over wavefronts
      wavefront: DO p = 1, nplanes
        ! Loop over cells in the wavefront
        cells: DO c = 1, planes(p)%num_cells

          ! Get the cell index
          IF ( istep > 0 ) THEN
            i = planes(p)%cells(c)%i
          ELSE
            i = nx - planes(p)%cells(c)%i + 1
          END IF
          IF ( jstep > 0) THEN
            j = planes(p)%cells(c)%j
          ELSE
            j = nx - planes(p)%cells(c)%j + 1
          END IF
          IF ( kstep > 0) THEN
            k = planes(p)%cells(c)%k
          ELSE
            k = nx - planes(p)%cells(c)%k + 1
          END IF

          ! Loop over the energy groups
          groups: DO g = 1, ng
            ! Loop over the angles in the octant
            angles: DO a = 1, nang

              !stuff

            END DO angles
          END DO groups
        END DO cells
      END DO wavefront

      CALL zero_edge_flux

    END DO octants

    ! Deallocate the plane list
    DO a = p, nplanes
      DEALLOCATE ( planes(p)%cells )
    END DO
    DEALLOCATE ( planes )


  END SUBROUTINE omp_sweep


  SUBROUTINE compute_sweep_order ( planes, nplanes )

    TYPE(plane), DIMENSION(:) :: planes
    INTEGER(i_knd) :: nplanes, i, j, k, n, idx

    DO i = 1, nplanes
      planes%num_cells = 0
    END DO

    ! Cells in each plane have equal co-ordinate sum
    DO k = 1, nz_gl
      DO j = 1, ny_gl
        DO i = 1, nx
          planes(i+j+k-2)%num_cells = planes(i+j+k-2)%num_cells + 1
        END DO
      END DO
    END DO
    
    ! Allocate the memory for each plane
    DO i = 1, nplanes
      ALLOCATE ( planes(i)%cells(planes(i)%num_cells) )
      planes(i)%index = 1
    END DO

    ! Store the indexes in the plane array
    DO k = 1, nz_gl
      DO j = 1, ny_gl
        DO i = 1, nx
          n = i + j + k - 2
          idx = planes(n)%index
          planes(n)%cells(idx)%i = i
          planes(n)%cells(idx)%j = j
          planes(n)%cells(idx)%k = k
          planes(n)%index = planes(n)%index + 1
        END DO
      END DO
    END DO


  END SUBROUTINE compute_sweep_order    


END MODULE omp_sweep_module 
