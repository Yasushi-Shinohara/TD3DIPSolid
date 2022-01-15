Subroutine Writeout_OMPinfo(Nk)
  Use Omp_lib
  Implicit none
  Integer, intent(in) :: Nk
  Integer :: NOMP

  !$ NOMP = omp_get_max_threads()
  !$ Write(*,'("# ++++++++++++++++++++++++++++")') 
  !$ Write(*,'("# Number of OpenMP thread =",i8)')  NOMP
  !$ If (Nk < NOMP) then
  !$   Write(*,'("# CAUTION: OMP parallelization is not efficient because of Nk < NOMP.")')
  !$ End If
  !$ Write(*,'("# ++++++++++++++++++++++++++++")') 
 
  Return
End Subroutine Writeout_OMPinfo
!==========================================================================================
Subroutine uGbk_forward_RK4(uGbk, hGGk, NG, Nocc, Nk, dt)
  Implicit none
  Complex(kind(0d0)), parameter :: zI=(0.0d0, 1.0d0)
  Integer, intent(in) :: NG, Nocc, Nk
  Double precision, intent(in) :: dt
  Complex(kind(0d0)), intent(inout) :: uGbk(1:Nk,1:Nocc,1:NG)
  Complex(kind(0d0)), intent(in) :: hGGk(1:Nk,1:NG,1:NG)
  Integer :: ik, ig, ib
  Complex(kind(0d0)) :: k1(1:Nocc,1:NG), k2(1:Nocc,1:NG), k3(1:Nocc,1:NG), k4(1:Nocc,1:NG)

  !$omp parallel
  !$omp do private(ik,k1,k2,k3,k4)
  Do ik = 1, Nk
    k1 = u_h2hu(uGbk(ik,:,:)              , hGGk(ik,:,:))/zI
    k2 = u_h2hu(uGbk(ik,:,:) + 0.5d0*dt*k1, hGGk(ik,:,:))/zI
    k3 = u_h2hu(uGbk(ik,:,:) + 0.5d0*dt*k2, hGGk(ik,:,:))/zI
    k4 = u_h2hu(uGbk(ik,:,:) + dt*k3      , hGGk(ik,:,:))/zI
    uGbk(ik,:,:) = uGbk(ik,:,:) + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)*dt/6.0d0
  End Do
  !$omp end do
  !$omp end parallel
  Return
Contains
  Function u_h2hu(u, h) result(hu)
    Implicit none
    Complex(kind(0d0)), intent(in) :: u(1:Nocc, 1:NG), h(1:NG, 1:NG)
    Complex(kind(0d0)) :: hu(1:Nocc, 1:NG)
    Integer :: i

!    hu(1:Nocc, 1:NG) = 0.d0
!    Do ig = 1,NG
!      Do ib = 1,Nocc
!        Do i = 1,NG
!          hu(ib,ig) = hu(ib,ig) + h(i,ig)*u(ib,i)
!        End Do
!      End Do 
!    End Do
    hu = matmul(u,h)
  End Function u_h2hu
End Subroutine uGbk_forward_RK4
!==========================================================================================
Subroutine uGbk_forward_RK4FFT(uGbk, tGGk, vx, NG, Nocc, Nk, dt)
  Implicit none
  Complex(kind(0d0)), parameter :: zI=(0.0d0, 1.0d0)
  Integer, intent(in) :: NG, Nocc, Nk
  Double precision, intent(in) :: dt
  Complex(kind(0d0)), intent(inout) :: uGbk(1:Nk,1:Nocc,1:NG)
  Complex(kind(0d0)), intent(in) :: tGGk(1:Nk,1:NG,1:NG)
  Double precision, intent(in) :: vx(1:NG)
  Integer :: ik, ig, ib
  Complex(kind(0d0)) :: k1(1:Nocc,1:NG), k2(1:Nocc,1:NG), k3(1:Nocc,1:NG), k4(1:Nocc,1:NG)

  !$omp parallel
  !$omp do private(ik,k1,k2,k3,k4)
  Do ik = 1, Nk
    k1 = u_t_v2hu_FFT(uGbk(ik,:,:)              , tGGk(ik,:,:), vx)/zI
    k2 = u_t_v2hu_FFT(uGbk(ik,:,:) + 0.5d0*dt*k1, tGGk(ik,:,:), vx)/zI
    k3 = u_t_v2hu_FFT(uGbk(ik,:,:) + 0.5d0*dt*k2, tGGk(ik,:,:), vx)/zI
    k4 = u_t_v2hu_FFT(uGbk(ik,:,:) + dt*k3      , tGGk(ik,:,:), vx)/zI
    uGbk(ik,:,:) = uGbk(ik,:,:) + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)*dt/6.0d0
  End Do
  !$omp end do
  !$omp end parallel
  Return
Contains
  Function u_t_v2hu_FFT(u, t, vx) result(hu)
    Use, intrinsic :: iso_c_binding
    Implicit none
    Complex(kind(0d0)), intent(in) :: u(1:Nocc, 1:NG), t(1:NG, 1:NG)
    Double precision, intent(in) :: vx(1:NG)
    Complex(kind(0d0)) :: hu(1:Nocc, 1:NG), tu(1:Nocc, 1:NG), ux(1:Nocc, 1:NG), vux(1:Nocc, 1:NG), vu(1:Nocc, 1:NG)
    Integer :: i, ib
!For FFTW
    Type(C_PTR) :: planf, planb
    Complex(C_DOUBLE_COMPLEX) :: uwork(1:Nocc, 1:NG)
    Include 'fftw3.f03'

    planf = fftw_plan_dft_1d(NG, uwork(ib,:), ux(ib,:), FFTW_FORWARD, FFTW_ESTIMATE)
    planb = fftw_plan_dft_1d(NG, ux(ib,:), uwork(ib,:), FFTW_BACKWARD, FFTW_ESTIMATE)

    uwork(:,:) = u(:,:) !INTENT = IN variable can not be substituted directly in the fftw_execute_dft without clear reasons
    Do ig = 1,NG
      tu(:,ig) = t(ig,ig)*u(:,ig)
    End Do
    Do ib = 1,Nocc
      Call fftw_execute_dft(planb,uwork(ib,:),ux(ib,:)) !iFFT ux <= u
      vux(ib,:) = vx(:)*ux(ib,:)
      Call fftw_execute_dft(planf,vux(ib,:),vu(ib,:)) !FFT vux => vu
    End Do
    vu = vu/float(NG)
    
    hu = tu + vu
    Call fftw_destroy_plan(planf)
    Call fftw_destroy_plan(planb)
  End Function u_t_v2hu_FFT
  !==
  Function u_t_v2hu_FFT_old(u, t, vx) result(hu)
    Implicit none
    Complex(kind(0d0)), intent(in) :: u(1:Nocc, 1:NG), t(1:NG, 1:NG)
    Double precision, intent(in) :: vx(1:NG)
    Complex(kind(0d0)) :: hu(1:Nocc, 1:NG), tu(1:Nocc, 1:NG), ux(1:Nocc, 1:NG), vux(1:Nocc, 1:NG), vu(1:Nocc, 1:NG)
    Integer :: i, ib
!For FFTW
    Integer(8) :: planf, planb
    Complex(kind(0d0)) :: work1(1:NG), work2(1:NG)
    Include 'fftw3.f'

    Call dfftw_plan_dft_1d(planf, NG, work1, work2, FFTW_FORWARD, FFTW_ESTIMATE)
    Call dfftw_plan_dft_1d(planb, NG, work2, work1, FFTW_BACKWARD, FFTW_ESTIMATE)
    Do ig = 1,NG
      tu(:,ig) = t(ig,ig)*u(:,ig)
    End Do
    Do ib = 1,Nocc
      Call dfftw_execute_dft(planb,u(ib,:),ux(ib,:)) !iFFT ux <= u
      vux(ib,:) = vx(:)*ux(ib,:)
      Call dfftw_execute_dft(planf,vux(ib,:),vu(ib,:)) !FFT vux => vu
    End Do
    vu = vu/float(NG)
    
    hu = tu + vu
    Call dfftw_destroy_plan(planf)
    Call dfftw_destroy_plan(planb)
  End Function u_t_v2hu_FFT_old
End Subroutine uGbk_forward_RK4FFT
!==========================================================================================
Subroutine uGbk_forward_exp(uGbk, hGGk, NG, Nocc, Nk, dt)
  Implicit none
  Complex(kind(0d0)), parameter :: zI=(0.0d0, 1.0d0)
  Integer, intent(in) :: NG, Nocc, Nk
  Double precision, intent(in) :: dt
  Complex(kind(0d0)), intent(inout) :: uGbk(1:Nk,1:Nocc,1:NG)
  Complex(kind(0d0)), intent(in) :: hGGk(1:Nk,1:NG,1:NG)
  Integer :: ik, ig, ib
  Complex(kind(0d0)) :: U(1:NG,1:NG)

  !$omp parallel
  !$omp do private(ik,U)
  Do ik = 1, Nk
    U = hdt2U(hGGk(ik,:,:)*dt)
    uGbk(ik,:,:) = U_u2Uu(U, uGbk(ik,:,:))
  End Do
  !$omp end do
  !$omp end parallel
  Return
Contains
  !==
  Function hdt2U(hdt) result(U)
    Implicit none
    Complex(kind(0d0)), intent(in) :: hdt(1:NG, 1:NG)
    Complex(kind(0d0)) :: U(1:NG, 1:NG)
    Complex(kind(0d0)) :: coef(1:NG, 1:NG)
    Double precision :: eigs(1:NG)
    Integer :: i,j,k

    U = (0.d0, 0.d0)
    Call eigh(NG,transpose(hdt),eigs,coef) !Transpose is for the conversion from Row-major to Column-major
    Do i = 1,NG
      Do J = 1,NG
        Do k = 1,NG
          U(j,k) = U(j,k) + exp(-zI*eigs(i))*coef(j,i)*conjg(coef(k,i))
        End Do
      End Do
    End Do
    U = transpose(U)                      !Transpose is for the conversion from Column-major to Row-major 
  End Function hdt2U
  !==
  Function U_u2Uu(Unitary, uorb) result(Unitaryuorb)
    Implicit none
    Complex(kind(0d0)), intent(in) :: uorb(1:Nocc, 1:NG), Unitary(1:NG, 1:NG)
    Complex(kind(0d0)) :: Unitaryuorb(1:Nocc, 1:NG)
    Integer :: i 

!    Uu(1:Nocc, 1:NG) = 0.d0
!    Do ig = 1,NG
!      Do ib = 1,Nocc
!        Do i = 1,NG
!          Uu(ib,jg) = Uu(ib,ig) + U(i,ig)*u(ib,i)
!        End Do
!      End Do 
!    End Do
    Unitaryuorb = matmul(uorb,Unitary)
  End Function U_u2Uu
End Subroutine uGbk_forward_exp
!==========================================================================================
Subroutine eigh(N,H,E,V)
  Implicit none
  Integer, intent(in) :: N
  Complex(kind(0d0)), intent(in) :: H(1:N, 1:N)
  Double precision, intent(out) :: E(1:N)
  Complex(kind(0d0)), intent(out) :: V(1:N, 1:N)
!For ZHEEV
  Integer :: LWORK_EV 
  Complex(kind(0d0)) :: WORK_EV(2*(2*N-1)) !The argument is just LWORK_EV
  Double precision :: RWORK_EV(3*N - 2)
  Integer :: INFO_EV
  LWORK_EV = 2*(2*N-1)

  V = H
  Call ZHEEV('V','U',N,V,N,E,WORK_EV,LWORK_EV,RWORK_EV,INFO_EV)
    
  Return
End Subroutine eigh
