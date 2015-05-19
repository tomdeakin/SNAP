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

  USE geom_module, ONLY: nx, ny_gl, nz_gl, hi, hj, hk, dinv

  USE data_module, ONLY: ng, vdelt

  USE sn_module, ONLY: nang, noct, ec, cmom, mu

  USE solvar_module, ONLY: qtot

  USE time_module, ONLY: wtime

  IMPLICIT NONE

  TYPE cell
    INTEGER(i_knd) :: i, j, k
  END TYPE cell

  TYPE plane
    INTEGER(i_knd) :: num_cells
    INTEGER(i_knd) :: index
    TYPE(cell), ALLOCATABLE, DIMENSION(:) :: cells
  END TYPE plane

  PUBLIC

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
  !$omp workshare
    flux_i(:,:,:,:) = 0.0
    flux_j(:,:,:,:) = 0.0
    flux_k(:,:,:,:) = 0.0
  !$omp end workshare
  END SUBROUTINE zero_edge_flux



  SUBROUTINE omp_sweep
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: o, p, c, g, a, nplanes
    INTEGER(i_knd) :: xhi, yhi, zhi, istep, jstep, kstep
    INTEGER(i_knd) :: i, j, k, l
    REAL(r_knd) :: source, psi
    TYPE(plane), ALLOCATABLE, DIMENSION(:) :: planes
    REAL(r_knd) :: t1, t2

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

    
    !$omp workshare
    flux_in(:,:,:,:,:,:) = 0.0
    flux_out(:,:,:,:,:,:) = 0.0
    !$omp end workshare

    CALL wtime ( t1 )

    ! Loop over octants
    octants: DO o = 1, noct
      CALL zero_edge_flux

      ! Calculate the corner co-ordinates for this octant
      SELECT CASE (o)
        CASE (1)
          xhi = nx
          yhi = ny_gl
          zhi = nz_gl
        CASE (2)
          xhi = 1
          yhi = ny_gl
          zhi = nz_gl
        CASE (3)
          xhi = nx
          yhi = 1
          zhi = nz_gl
        CASE (4)
          xhi = 1
          yhi = 1
          zhi = nz_gl
        CASE (5)
          xhi = nx
          yhi = ny_gl
          zhi = 1
        CASE (6)
          xhi = 1
          yhi = ny_gl
          zhi = 1
        CASE (7)
          xhi = nx
          yhi = 1
          zhi = 1
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

      CALL omp_sweep_c( p, ng, nang, nx, ny_gl, nz_gl, istep, jstep, kstep, o, noct, cmom, &
                         qtot, ec, mu, hi, hj, hk, vdelt, dinv, &
                         flux_i, flux_j, flux_k, flux_in, flux_out )

      ! CALL THE C!!!!!!

        ! Loop over cells in the wavefront
        !!$omp parallel
        !!dir$ ivdep
        !!$omp do private(i,j,k,psi,source)
        ! cells: DO c = 1, planes(p)%num_cells

        !   ! Get the cell index
        !   IF ( istep > 0 ) THEN
        !     i = planes(p)%cells(c)%i
        !   ELSE
        !     i = nx - planes(p)%cells(c)%i + 1
        !   END IF
        !   IF ( jstep > 0) THEN
        !     j = planes(p)%cells(c)%j
        !   ELSE
        !     j = ny_gl - planes(p)%cells(c)%j + 1
        !   END IF
        !   IF ( kstep > 0) THEN
        !     k = planes(p)%cells(c)%k
        !   ELSE
        !     k = nz_gl - planes(p)%cells(c)%k + 1
        !   END IF

        !   ! Loop over the energy groups
        !   !dir$ ivdep
        !   groups: DO g = 1, ng
        !     ! Loop over the angles in the octant
        !     !dir$ ivdep
        !     angles: DO a = 1, nang
        !       ! Compute the angular flux
        !       source = qtot(1,i,j,k,g)
        !       !dir$ novector
        !       DO l = 2, cmom
        !         source = source + ec(a,l,o) * qtot(l,i,j,k,g)
        !       END DO

        !       psi = source + (flux_i(a,g,j,k) * mu(a) * hi) + (flux_j(a,g,i,k) * hj(a)) + (flux_k(a,g,i,j) * hk(a))

        !       IF (vdelt(g) /= 0.0) psi = psi + vdelt(g) * flux_in(a,g,i,j,k,o)

        !       psi = psi * dinv(a,i,j,k,g)

        !       ! TODO: fixup

        !       flux_i(a,g,j,k) = 2.0 * psi - flux_i(a,g,j,k)
        !       flux_j(a,g,i,k) = 2.0 * psi - flux_j(a,g,i,k)
        !       flux_k(a,g,i,j) = 2.0 * psi - flux_k(a,g,i,j)

        !       IF (vdelt(g) /= 0.0) psi = (2.0 * psi) - flux_in(a,g,i,j,k,o)

        !       flux_out(a,g,i,j,k,o) = psi

        !     END DO angles

        !     ! Reduce to the scalar flux
        !     ! TODO

        !   END DO groups
        ! END DO cells
      !!$omp end do
      !!$omp end parallel
      END DO wavefront


    END DO octants

    CALL wtime ( t2 )

    ! Deallocate the plane list
    DO a = p, nplanes
      DEALLOCATE ( planes(p)%cells )
    END DO
    DEALLOCATE ( planes )

    PRINT *, "Sweep took", t2-t1

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

  SUBROUTINE transpose_angular_flux ( flux_in, flux_out )

    REAL(r_knd), DIMENSION(:,:,:,:,:,:), POINTER :: flux_in, flux_out
    INTEGER(i_knd) :: o, g, a

    DO o = 1, noct
      DO g = 1, ng
        DO a = 1, nang
          flux_out(a,:,:,:,o,g) = flux_in(a,g,:,:,:,o)
        END DO
      END DO
    END DO

  END SUBROUTINE transpose_angular_flux


END MODULE omp_sweep_module 
