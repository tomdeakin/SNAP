!-----------------------------------------------------------------------
!
! MODULE: solvar_module
!> @brief
!> This module contains several variables that are used in the solution
!> process, including their allocation and deallocation. Also includes
!> initialization of sweep parameters.
! 
!-----------------------------------------------------------------------

MODULE solvar_module

  USE global_module, ONLY: i_knd, r_knd, zero

  USE plib_module, ONLY: ichunk

  USE geom_module, ONLY: nx, ny, nz, nc

  USE sn_module, ONLY: nang, noct, nmom, cmom

  USE data_module, ONLY: ng

  USE control_module, ONLY: timedep, angcpy

  USE ISO_C_BINDING

  IMPLICIT NONE

!-----------------------------------------------------------------------
! C function to allocate to 2MB boundaries
!-----------------------------------------------------------------------
  INTERFACE
    TYPE(C_PTR) FUNCTION ALLOC(len) BIND(C,name="alloc")
      IMPORT :: C_PTR
      IMPLICIT NONE
      INTEGER :: len
    END FUNCTION
  END INTERFACE

  PUBLIC

  SAVE
!_______________________________________________________________________
!
! Module variables
!
! ptr_in(nang,nx,ny,nz,noct,ng)   - Incoming time-edge flux pointer
! ptr_out(nang,nx,ny,nz,noct,ng)  - Outgoing time-edge flux pointer
!
! flux0(nx,ny,nz,ng)         - Scalar flux moments array
! flux0po(nx,ny,nz,ng)       - Previous outer copy of scalar flux array
! flux0pi(nx,ny,nz,ng)       - Previous inner copy of scalar flux array
! fluxm(cmom-1,nx,ny,nz,ng)  - Flux moments array
!
! q2grp0(nx,ny,nz,ng)        - Isotropic out-of-group + fixed sources
! q2grpm(cmom-1,nx,ny,nz,ng) - Anisotropic out-of-group + fixed sources
! qtot(cmom,ichunk,ny,nz,nc,ng) - Total source: q2grp0 + q2grpm +
!                                 within-group source
!
! t_xs(nx,ny,nz,ng)       - Total cross section on mesh
! a_xs(nx,ny,nz,ng)       - Absorption cross section on mesh
! s_xs(nx,ny,nz,nmom,ng)  - In-group scattering cross section on mesh
!
! psii(nang,ny,nz,ng)     - Working psi_x array
! psij(nang,ichunk,nz,ng) - Working psi_y array
! psik(nang,ichunk,ny,ng) - Working psi_z array
!
! jb_in(nang,ichunk,nz,ng)  - y-dir boundary flux in from comm
! jb_out(nang,ichunk,nz,ng) - y-dir boundary flux out to comm
! kb_in(nang,ichunk,ny,ng)  - z-dir boundary flux in from comm
! kb_out(nang,ichunk,ny,ng) - z-dir boundary flux out to comm
!
! flkx(nx+1,ny,nz,ng)     - x-dir leakage array
! flky(nx,ny+1,nz,ng)     - y-dir leakage array
! flkz(nx,ny,nz+1,ng)     - z-dir leakage array
!
! fmin(ng)       - dummy flux min
! fmax(ng)       - dummy flux max
!
! pop(ng)        - particle population spectrum
!
!_______________________________________________________________________

  REAL(r_knd), ALLOCATABLE, DIMENSION(:) :: fmin, fmax, pop

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:,:,:) :: flux0, flux0po,      &
    flux0pi, q2grp0, t_xs, a_xs, psii, psij, psik, jb_in, jb_out,      &
    kb_in, kb_out, flkx, flky, flkz

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: q2grpm, fluxm, s_xs

  REAL(r_knd), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: qtot

  REAL(r_knd), DIMENSION(:,:,:,:,:,:), POINTER :: ptr_in, ptr_out


  CONTAINS


  SUBROUTINE solvar_allocate ( ierr )

!-----------------------------------------------------------------------
!
! Allocate solution arrays.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(OUT) :: ierr
    TYPE(C_PTR) :: tmp_ptr
!_______________________________________________________________________
!
!   Allocate ptr_in/out if needed. If angcpy=1, only allocate for one
!   copy and point the other pointer at it. If angcpy=2, allocate for
!   two copies. Provide an initial condition of zero. This may be
!   changed in the future if necessary.
!_______________________________________________________________________

    ierr = 0

    NULLIFY( ptr_in, ptr_out )

    IF ( timedep == 1 ) THEN
      IF ( angcpy == 1 ) THEN
        tmp_ptr = ALLOC( nang*nx*ny*nz*noct*ng )
        CALL C_F_POINTER( tmp_ptr, ptr_in, (/ nang,nx,ny,nz,noct,ng /) )
        !ALLOCATE( ptr_in(nang,nx,ny,nz,noct,ng), STAT=ierr )
      ELSE
        ALLOCATE( ptr_in(nang,nx,ny,nz,noct,ng),                       &
          ptr_out(nang,nx,ny,nz,noct,ng), STAT=ierr )
      END IF
    ELSE
      ALLOCATE( ptr_in(0,0,0,0,0,0), ptr_out(0,0,0,0,0,0), STAT=ierr )
    END IF
    IF ( ierr /= 0 ) RETURN

    IF ( timedep == 1 ) THEN
      ptr_in = zero
      IF ( angcpy == 1 ) THEN
        ptr_out => ptr_in
      ELSE
        ptr_out = zero
      END IF
    END IF
!_______________________________________________________________________
!
!   Allocate the flux moments arrays. Keep an old copy. If isotropic,
!   allocate fluxm as a dummy array to make passing contiguous pieces of
!   it in argument lists possible even in debug mode. There are better
!   ways to do this, but want to keep data structures simple for others
!   to change as they want easily.
!_______________________________________________________________________

    ALLOCATE( flux0(nx,ny,nz,ng), flux0po(nx,ny,nz,ng),                &
      flux0pi(nx,ny,nz,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN

    IF ( cmom > 1 ) THEN
      ALLOCATE( fluxm(cmom-1,nx,ny,nz,ng), STAT=ierr )
      IF ( ierr /= 0 ) RETURN
    ELSE
      ALLOCATE( fluxm(0:0,nx,ny,nz,ng), STAT=ierr )
      IF ( ierr /= 0 ) RETURN
    END IF

    flux0   = zero
    flux0po = zero
    flux0pi = zero
    fluxm   = zero
!_______________________________________________________________________
!
!   Allocate the source arrays. Do the same thing for q2grpm as was done
!   to fluxm above.
!_______________________________________________________________________

    ALLOCATE( q2grp0(nx,ny,nz,ng), qtot(cmom,ichunk,ny,nz,nc,ng),      &
      STAT=ierr )
    IF ( ierr /= 0 ) RETURN

    IF ( cmom > 1 ) THEN
      ALLOCATE( q2grpm(cmom-1,nx,ny,nz,ng), STAT=ierr )
      IF ( ierr /= 0 ) RETURN
    ELSE
      ALLOCATE( q2grpm(0:0,nx,ny,nz,ng), STAT=ierr )
      IF ( ierr /= 0 ) RETURN
    END IF

    q2grp0 = zero
    q2grpm = zero
    qtot = zero
!_______________________________________________________________________
!
!   Allocate the cross section expanded to spatial mesh arrays
!_______________________________________________________________________

    ALLOCATE( t_xs(nx,ny,nz,ng), a_xs(nx,ny,nz,ng),                    &
      s_xs(nx,ny,nz,nmom,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN

    t_xs = zero
    a_xs = zero
    s_xs = zero
!_______________________________________________________________________
!
!   Working arrays
!_______________________________________________________________________

    tmp_ptr = ALLOC( nang*ny*nz*ng )
    CALL C_F_POINTER( tmp_ptr, psii, (/ nang, ny, nz, ng /) )

    tmp_ptr = ALLOC( nang*ichunk*nz*ng )
    CALL C_F_POINTER( tmp_ptr, psij, (/ nang, ichunk, nz, ng /) )

    tmp_ptr = ALLOC( nang*ichunk*ny*ng )
    CALL C_F_POINTER( tmp_ptr, psik, (/ nang, ichunk, ny, ng /) )

    !ALLOCATE( psii(nang,ny,nz,ng), psij(nang,ichunk,nz,ng),            &
    !  psik(nang,ichunk,ny,ng), STAT=ierr )
    !IF ( ierr /= 0 ) RETURN

    psii = zero
    psij = zero
    psik = zero
!_______________________________________________________________________
!
!   PE boundary flux arrays
!_______________________________________________________________________

    ALLOCATE( jb_in(nang,ichunk,nz,ng), jb_out(nang,ichunk,nz,ng),     &
      kb_in(nang,ichunk,ny,ng), kb_out(nang,ichunk,ny,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN

    jb_in  = zero
    jb_out = zero
    kb_in  = zero
    kb_out = zero
!_______________________________________________________________________
!
!   Leakage arrays
!_______________________________________________________________________

    ALLOCATE( flkx(nx+1,ny,nz,ng), flky(nx,ny+1,nz,ng),                &
      flkz(nx,ny,nz+1,ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN

    flkx = zero
    flky = zero
    flkz = zero
!_______________________________________________________________________
!
!   Flux extremes and particle population spectrum
!_______________________________________________________________________

    ALLOCATE( fmin(ng), fmax(ng), pop(ng), STAT=ierr )
    IF ( ierr /= 0 ) RETURN

    fmin = zero
    fmax = zero
    pop = zero
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE solvar_allocate


  SUBROUTINE solvar_deallocate

!-----------------------------------------------------------------------
!
! Deallocate solve_module arrays.
!
!-----------------------------------------------------------------------
!_______________________________________________________________________

    DEALLOCATE( ptr_in )
    IF ( angcpy==2 .OR. timedep==0 ) DEALLOCATE( ptr_out )
    DEALLOCATE( flux0, flux0po, flux0pi, fluxm )
    DEALLOCATE( q2grp0, q2grpm, qtot )
    DEALLOCATE( t_xs, a_xs, s_xs )
    DEALLOCATE( psii, psij, psik )
    DEALLOCATE( jb_in, jb_out, kb_in, kb_out )
    DEALLOCATE( flkx, flky, flkz )
    DEALLOCATE( fmin, fmax, pop )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE solvar_deallocate


END MODULE solvar_module
