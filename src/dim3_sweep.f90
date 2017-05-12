!-----------------------------------------------------------------------
!
! MODULE: dim3_sweep_module
!> @brief
!> This module contains the 2D and 3D mesh sweep logic.
!
!-----------------------------------------------------------------------

MODULE dim3_sweep_module

  USE global_module, ONLY: i_knd, l_knd, r_knd, zero, two, one, half

  USE plib_module, ONLY: ichunk, firsty, lasty, firstz, lastz

  USE geom_module, ONLY: nx, hi, hj, hk, ndimen, ny, nz, nc

  USE sn_module, ONLY: cmom, nang, mu, eta, xi, w

  USE data_module, ONLY: src_opt, qim

  USE control_module, ONLY: fixup, tolr, last_oct, update_ptr

  USE thrd_comm_module, ONLY: sweep_recv_bdry, sweep_send_bdry

  IMPLICIT NONE

  PUBLIC :: dim3_sweep


  CONTAINS


  SUBROUTINE dim3_sweep ( ich, id, d1, d2, d3, d4, jd, kd, oct, g, t,  &
    iop, reqs, szreq, psii, psij, psik, qtot, ec, vdelt, ptr_in,       &
    ptr_out, dinv, flux0, fluxm, jb_in, jb_out, kb_in, kb_out, wmu,    &
    weta, wxi, flkx, flky, flkz, t_xs, fmin, fmax )

!-----------------------------------------------------------------------
!
! 3-D slab mesh sweeper.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: ich, id, d1, d2, d3, d4, jd, kd, oct,&
      g, t, iop, szreq

    INTEGER(i_knd), DIMENSION(szreq), INTENT(INOUT) :: reqs

    REAL(r_knd), INTENT(IN) :: vdelt

    REAL(r_knd), INTENT(INOUT) :: fmin, fmax

    REAL(r_knd), DIMENSION(nang), INTENT(IN) :: wmu, weta, wxi

    REAL(r_knd), DIMENSION(nang,cmom), INTENT(IN) :: ec

    REAL(r_knd), DIMENSION(nang,ny,nz), INTENT(INOUT) :: psii

    REAL(r_knd), DIMENSION(nang,ichunk,nz), INTENT(INOUT) :: psij,     &
      jb_in, jb_out

    REAL(r_knd), DIMENSION(nang,ichunk,ny), INTENT(INOUT) :: psik,     &
      kb_in, kb_out

    REAL(r_knd), DIMENSION(nx,ny,nz), INTENT(IN) :: t_xs

    REAL(r_knd), DIMENSION(nx,ny,nz), INTENT(INOUT) :: flux0

    REAL(r_knd), DIMENSION(nx+1,ny,nz), INTENT(INOUT) :: flkx

    REAL(r_knd), DIMENSION(nx,ny+1,nz), INTENT(INOUT) :: flky

    REAL(r_knd), DIMENSION(nx,ny,nz+1), INTENT(INOUT) :: flkz

    REAL(r_knd), DIMENSION(nang,ichunk,ny,nz), INTENT(IN) :: dinv

    REAL(r_knd), DIMENSION(cmom-1,nx,ny,nz), INTENT(INOUT) :: fluxm

    REAL(r_knd), DIMENSION(cmom,ichunk,ny,nz), INTENT(IN) :: qtot

    REAL(r_knd), DIMENSION(d1,d2,d3,d4), INTENT(IN) :: ptr_in

    REAL(r_knd), DIMENSION(d1,d2,d3,d4), INTENT(OUT) :: ptr_out
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: ist, iclo, ichi, jst, jlo, jhi, kst, klo, khi, k,&
      j, ic, i, l, ibl, ibr, ibb, ibt, ibf, ibk, id1, a

    LOGICAL(l_knd) :: receive

    REAL(r_knd), DIMENSION(nang) :: sum_hv

    REAL(r_knd), DIMENSION(nang) :: psi, pc, den

    REAL(r_knd), DIMENSION(nang,4) :: hv, fxhv

! Only required because OpenMP cannot cope
! "Subobjects are not allowed in this OpenMP* clause; a named variable must be
! specified"
    REAL(r_knd) :: cell_flux0
    REAL(r_knd), DIMENSION(cmom-1) :: cell_fluxm
!_______________________________________________________________________
!
!   Set the receive control flag to true. Will set to false within the
!   loop for the first cell. Pushing the receive to this point helps
!   overlap communication and computation.
!_______________________________________________________________________

    receive = .TRUE.
!_______________________________________________________________________
!
!   Set up the sweep order given octant info.
!_______________________________________________________________________

    IF ( id == 1 ) THEN
      ist = -1; iclo = ichunk; ichi = 1
    ELSE
      ist = +1; iclo = 1; ichi = ichunk
    END IF

    IF ( jd == 1 ) THEN
      jst = -1; jlo = ny; jhi = 1
    ELSE
      jst = +1; jlo = 1; jhi = ny
    END IF

    IF ( kd == 1 ) THEN
      kst = -1; klo = nz; khi = 1
    ELSE
      kst = +1; klo = 1; khi = nz
    END IF
!_______________________________________________________________________
!
!   Initialize the fixup counter
!_______________________________________________________________________

    fxhv = zero
!_______________________________________________________________________
!
!   Loop over the cells using bounds/stride above
!_______________________________________________________________________

    k_loop:  DO k = klo, khi, kst
    j_loop:  DO j = jlo, jhi, jst
    ic_loop: DO ic = iclo, ichi, ist
!_______________________________________________________________________
!
!     Index ic refers to location in current ichunk, ich. Compute the
!     global i index for x-direction.
!_______________________________________________________________________

      i = (ich-1)*ichunk + ic

!_______________________________________________________________________
!
!     Dummy operation to match real transport code where nx is not
!     always evenly divided by ichunk
!_______________________________________________________________________

      IF ( i > nx ) THEN
        jb_out(:,ic,k) = zero
        kb_out(:,ic,j) = zero
        CYCLE ic_loop
      END IF
!_______________________________________________________________________
!
! Set the boundary conditions ahead of the solve - don't want MPI calls
! in the middle of the kernel.
!_______________________________________________________________________

!_______________________________________________________________________
!
!     Left/right boundary conditions, always vacuum. Dummy operations
!     may be ignored.
!_______________________________________________________________________

      ibl = 0
      ibr = 0
      IF ( i==nx .AND. ist==-1 ) THEN
        psii(:,j,k) = zero
      ELSE IF ( i==1 .AND. ist==1 ) THEN
        SELECT CASE ( ibl )
          CASE ( 0 )
            psii(:,j,k) = zero
          CASE ( 1 )
            psii(:,j,k) = zero
        END SELECT
      END IF
!_______________________________________________________________________
!
!     Process boundaries in y- and z-directions come from preset
!     boundary conditions (vacuum) or from upstream neighbors. Must call
!     to receive any data. If upstream is boundary, no communication.
!     Only do this at first cell of chunk but receive all incoming data
!     for chunk. Resetting to .FALSE. skips unneeded communication call
!     for all other cells.
!_______________________________________________________________________

      IF ( receive ) THEN
        CALL sweep_recv_bdry ( g, jd, kd, iop, t, reqs, szreq, nc,     &
          nang, ichunk, ny, nz, jb_in, kb_in )
        receive = .FALSE.
      END IF
!_______________________________________________________________________
!
!     Top/bottom boundary conditions for chunk. Either received from
!     upstream neighbor or set to zero (always vacuum) if no usptream
!     neighbor.
!_______________________________________________________________________

      ibb = 0
      ibt = 0
      IF ( j == jlo ) THEN
        IF ( jd==1 .AND. lasty ) THEN
          psij(:,ic,k) = zero
        ELSE IF ( jd==2 .AND. firsty ) THEN
          SELECT CASE ( ibb )
            CASE ( 0 )
              psij(:,ic,k) = zero
            CASE ( 1 )
              psij(:,ic,k) = zero
          END SELECT
        ELSE
          psij(:,ic,k) = jb_in(:,ic,k)
        END IF
      END IF
!_______________________________________________________________________
!
!     Front/back boundary condtions. Either received from upstream
!     neighbor or set to zero (always vacuum) if no upstream neighbor.
!_______________________________________________________________________

      ibf = 0
      ibk = 0
      IF ( k == klo ) THEN
        IF ( (kd==1 .AND. lastz) .OR. ndimen<3 ) THEN
          psik(:,ic,j) = zero
        ELSE IF ( kd==2 .AND. firstz ) THEN
          SELECT CASE ( ibf )
            CASE ( 0 )
              psik(:,ic,j) = zero
            CASE ( 1 )
              psik(:,ic,j) = zero
          END SELECT
        ELSE
          psik(:,ic,j) = kb_in(:,ic,j)
        END IF
      END IF
!_______________________________________________________________________
!
! Zero the temp scalar flux variables for this cell.
!_______________________________________________________________________
      cell_flux0 = zero
      cell_fluxm = zero
!_______________________________________________________________________
!
! Loop over the angles in this octant. All compute is done in this loop
! to good vectorisation and cache reuse. This also removes the need for
! Fortran vector notation, and ensures the loops are fused which is
! not the case using the vector notation.
!_______________________________________________________________________

  !$OMP SIMD REDUCTION(+:cell_flux0, cell_fluxm) ALIGNED(ptr_out:64)
  !DIR$ VECTOR NONTEMPORAL(ptr_out)
      a_loop: DO a = 1, nang
!_______________________________________________________________________
!
!     Compute the angular source. Add the MMS contribution if necessary.
!_______________________________________________________________________

        psi(a) = qtot(1,ic,j,k)

        DO l = 2, cmom
          psi(a) = psi(a) + ec(a,l)*qtot(l,ic,j,k)
        END DO

        IF ( src_opt == 3 ) psi(a) = psi(a) + qim(a,i,j,k,oct,g)
!_______________________________________________________________________
!
!     Compute initial solution
!_______________________________________________________________________

        IF ( vdelt /= zero ) THEN
          pc(a) = ( psi(a) + psii(a,j,k)*mu(a)*hi + psij(a,ic,k)*eta(a)*hj +         &
            psik(a,ic,j)*xi(a)*hk + ptr_in(a,i,j,k)*vdelt ) * dinv(a,ic,j,k)
        ELSE
          pc(a) = ( psi(a) + psii(a,j,k)*mu(a)*hi + psij(a,ic,k)*eta(a)*hj +         &
            psik(a,ic,j)*xi(a)*hk ) * dinv(a,ic,j,k)
        END IF
!_______________________________________________________________________
!
!     Compute outgoing edges with diamond difference, no negative flux
!     fixup
!_______________________________________________________________________

        IF ( fixup == 0 ) THEN
  
          psi(a) = pc(a)
  
          psii(a,j,k) = two*psi(a) - psii(a,j,k)
          psij(a,ic,k) = two*psi(a) - psij(a,ic,k)
          IF ( ndimen == 3 ) psik(a,ic,j) = two*psi(a) - psik(a,ic,j)
          IF ( vdelt/=zero .AND. update_ptr ) THEN
            ptr_out(a,i,j,k) = two*psi(a) - ptr_in(a,i,j,k)
          END IF
  
        ELSE
!_______________________________________________________________________
!
!       Use negative flux fixup. Compute outgoing edges. If negative,
!       set to zero and rebalance. Multi-pass until all negativities
!       eliminated. Initialize counters and determine angles that will
!       need fixup first.
!_______________________________________________________________________

          hv(a,:) = one
          sum_hv(a) = SUM( hv(a,:) )
  
          fixup_loop: DO
  
            fxhv(a,1) = two*pc(a) - psii(a,j,k)
            fxhv(a,2) = two*pc(a) - psij(a,ic,k)
            IF ( ndimen == 3 ) fxhv(a,3) = two*pc(a) - psik(a,ic,j)
            IF ( vdelt /= zero ) fxhv(a,4) = two*pc(a) - ptr_in(a,i,j,k)
  
            WHERE ( fxhv(a,:) < zero ) hv(a,:) = zero
!_______________________________________________________________________
!
!         Exit loop when all angles are fixed up, i.e., no change in hv
!_______________________________________________________________________

            IF ( sum_hv(a) == SUM( hv(a,:) ) ) EXIT fixup_loop
            sum_hv(a) = SUM( hv(a,:) )
!_______________________________________________________________________
!
!         Recompute balance equation numerator and denominator and get
!         new cell average flux
!_______________________________________________________________________

            IF ( vdelt /= zero ) THEN
              pc(a) = psi(a) + half * ( psii(a,j,k)*mu(a)*hi*(one+hv(a,1)) +      &
                                  psij(a,ic,k)*eta(a)*hj*(one+hv(a,2)) +    &
                                  psik(a,ic,j)*xi(a)*hk*(one+hv(a,3)) +     &
                                  ptr_in(a,i,j,k)*vdelt*(one+hv(a,4)) )
              den(a) = t_xs(i,j,k) + mu(a)*hi*hv(a,1) + eta(a)*hj*hv(a,2) +       &
                xi(a)*hk*hv(a,3) + vdelt*hv(a,4)
            ELSE
              pc(a) = psi(a) + half * ( psii(a,j,k)*mu(a)*hi*(one+hv(a,1)) +      &
                                  psij(a,ic,k)*eta(a)*hj*(one+hv(a,2)) +    &
                                  psik(a,ic,j)*xi(a)*hk*(one+hv(a,3)) )
              den(a) = t_xs(i,j,k) + mu(a)*hi*hv(a,1) + eta(a)*hj*hv(a,2) +       &
                xi(a)*hk*hv(a,3)
            END IF
  
            IF ( pc(a) <= zero ) den(a) = zero
  
            IF( den(a) < tolr ) THEN
              pc(a) = zero
              den(a) = one
            END IF
  
            pc(a) = pc(a) / den(a)
  
          END DO fixup_loop
!_______________________________________________________________________
!
!       Fixup done, compute edges with resolved center
!_______________________________________________________________________
  
          psi(a) = pc(a)
  
          psii(a,j,k) = fxhv(a,1) * hv(a,1)
          psij(a,ic,k) = fxhv(a,2) * hv(a,2)
          IF ( ndimen == 3 ) psik(a,ic,j) = fxhv(a,3) * hv(a,3)
          IF ( vdelt/=zero .AND. update_ptr ) THEN
            ptr_out(a,i,j,k) = fxhv(a,4) * hv(a,4)
          END IF
  
        END IF
!_______________________________________________________________________
!
!     Compute the flux moments
!_______________________________________________________________________

        psi(a) = w(a)*psi(a)
  
        cell_flux0 = cell_flux0 + psi(a)
        DO l = 1, cmom-1
          cell_fluxm(l) = cell_fluxm(l) + ( ec(a,l+1)*psi(a) )
        END DO
!_______________________________________________________________________
!
! End the angle loop
! The remaining statements need all the angles
!_______________________________________________________________________

      END DO a_loop
!_______________________________________________________________________
!
! Save the scalar flux
!_______________________________________________________________________

      IF ( oct == 1 ) THEN
        flux0(i,j,k) = cell_flux0
        DO l = 1, cmom-1
          fluxm(l,i,j,k) = cell_fluxm(l)
        END DO
      ELSE
        flux0(i,j,k) = flux0(i,j,k) + cell_flux0
        DO l = 1, cmom-1
          fluxm(l,i,j,k) = fluxm(l,i,j,k) + cell_fluxm(l)
        END DO
      END IF
!_______________________________________________________________________
!
!     Save edge fluxes (dummy if checks for unused non-vacuum BCs)
!_______________________________________________________________________
 
      IF ( j == jhi ) THEN
        IF ( jd==2 .AND. lasty ) THEN
          CONTINUE
        ELSE IF ( jd==1 .AND. firsty ) THEN
          IF ( ibb == 1 ) CONTINUE
        ELSE
          jb_out(:,ic,k) = psij(:,ic,k)
        END IF
      END IF

      IF ( k == khi ) THEN
        IF ( kd==2 .AND. lastz ) THEN
          CONTINUE
        ELSE IF ( kd==1 .AND. firstz ) THEN
          IF ( ibf == 1 ) CONTINUE
        ELSE
          kb_out(:,ic,j) = psik(:,ic,j)
        END IF
      END IF
!_______________________________________________________________________
!
!     Compute dummy leakages (not used elsewhere currently)
!_______________________________________________________________________

      flkx(i+id-1,j,k) = flkx(i+id-1,j,k) + ist*SUM( wmu*psii(:,j,k) )
      flky(i,j+jd-1,k) = flky(i,j+jd-1,k) + jst*SUM( weta*psij(:,ic,k) )
      flkz(i,j,k+kd-1) = flkz(i,j,k+kd-1) + kst*SUM( wxi*psik(:,ic,j) )

!_______________________________________________________________________
!
!     Calculate dummy min and max scalar fluxes (not used elsewhere
!     currently)
!_______________________________________________________________________

      IF ( oct == last_oct ) THEN
        fmin = MIN( fmin, flux0(i,j,k) )
        fmax = MAX( fmax, flux0(i,j,k) )
      END IF
!_______________________________________________________________________
!
!   Finish the loops
!_______________________________________________________________________

    END DO ic_loop
    END DO j_loop
    END DO k_loop
!_______________________________________________________________________
!
!   Send data to downstream neighbors
!_______________________________________________________________________

    CALL sweep_send_bdry ( g, jd, kd, iop, t, reqs, szreq, nc, nang,   &
      ichunk, ny, nz, jb_out, kb_out )
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE dim3_sweep


END MODULE dim3_sweep_module
