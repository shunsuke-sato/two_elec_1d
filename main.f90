module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)

  real(8),parameter :: lambda_int = 1.0d0
  
! grid
  integer :: nx
  real(8) :: length_x, dx
  real(8),allocatable :: xn(:)

! wavefunction
  real(8),allocatable :: wfn(:,:),wfn_t(:,:),hwfn_t(:,:)
  real(8),allocatable :: v_ext(:), v_int(:,:), v_all(:,:)

  real(8) :: eps_gs
  real(8),allocatable :: rho_gs(:)


contains
  function erf_x(x) result(y)
    real(8),intent(in) :: x
    real(8),parameter :: epsilon_s = 1d-3
    real(8) :: y
    
    if(abs(x) > epsilon_s)then
      y = erf(x)/x
    else
      y = 2d0/sqrt(pi)*( 1d0 - x**2/3d0 + x**4/10d0 - x**6/42d0 + x**8/216d0)
    end if
    
  end function erf_x


end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none


  call input
  call initialize

  call cg_gs
  call calc_hartree_fock
!  call ks_inversion
!  call scf_ks

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  nx = 600
  length_x = 60d0
  dx = length_x/nx
  

end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: ix, iy
  real(8) :: r

  allocate(wfn(0:nx,0:nx))
  allocate(hwfn_t(0:nx,0:nx))
  allocate(wfn_t(0-4:nx+4,0-4:nx+4))
  wfn_t = 0d0

  allocate(rho_gs(0:nx))
  allocate(v_ext(0:nx),v_int(0:nx,0:nx),v_all(0:nx,0:nx))


  allocate(xn(0:nx))
  do ix = 0,nx
    xn(ix) = -0.5d0*length_x + dx*ix
  end do

! external potential
  do ix = 0,nx
    v_ext(ix) = -2d0/sqrt(0.5d0**2+xn(ix)**2)
!    v_ext(ix) = -1d0/sqrt(2d0+xn(ix)**2)
  end do

! e-e interaction
  do iy = 0,nx
    do ix = 0,nx
      v_int(ix,iy) = lambda_int/sqrt(0.5d0**2+(xn(ix)-xn(iy))**2)
!      v_int(ix,iy) = 0d0 
    end do
  end do


  do iy = 0,nx
    do ix = 0,nx
      v_all(ix,iy) = v_ext(ix) + v_ext(iy) + v_int(ix,iy)
    end do
  end do


end subroutine initialize
!-------------------------------------------------------------------------------
subroutine laplacian(wfn_in,Lwfn_out,fact)
  use global_variables
  real(8),intent(in) :: wfn_in(0:nx,0:nx), fact
  real(8),intent(out) :: Lwfn_out(0:nx,0:nx)
  real(8),parameter :: ct0 = -2d0, ct1 = 1d0
  real(8) :: c0,c1
  integer :: ix,iy

  c0 = 2d0*fact*ct0/dx**2
  c1 = fact*ct1/dx**2

  wfn_t(0:nx,0:nx) = wfn_in(0:nx,0:nx)

  do iy = 0,nx
    do ix = 0,nx

      Lwfn_out(ix,iy) = c0*wfn_t(ix,iy) &
        +c1*( wfn_t(ix+1,iy) + wfn_t(ix-1,iy) + wfn_t(ix,iy+1) + wfn_t(ix,iy-1) )

    end do
  end do

end subroutine laplacian
!-------------------------------------------------------------------------------
subroutine hpsi(wfn_in,hwfn_out)
  use global_variables
  real(8),intent(in) :: wfn_in(0:nx,0:nx)
  real(8),intent(out) :: hwfn_out(0:nx,0:nx)

  call laplacian(wfn_in,hwfn_out,-0.5d0)
  hwfn_out = hwfn_out + v_all*wfn_in

end subroutine hpsi
!-------------------------------------------------------------------------------
subroutine cg_gs
  use global_variables
  implicit none
  integer :: ncg, icg
  real(8),parameter :: res_epsilon = 1d-12
  real(8),allocatable :: xvec(:,:),gvec(:,:),pvec(:,:)
  real(8),allocatable :: hxvec(:,:),hpvec(:,:)
  real(8) :: gg0, gg, xx, xhx, pp, php, xp, xhp
  real(8) :: alpha, beta, aa, bb, cc, ss, lambda, res
  integer :: ix,iy

  ncg = 300

  allocate(xvec(0:nx,0:nx),pvec(0:nx,0:nx),gvec(0:nx,0:nx))
  allocate(hxvec(0:nx,0:nx),hpvec(0:nx,0:nx))

  do iy = 0,nx
    do ix = 0,nx
      wfn(ix,iy) = exp(-0.5d0*(xn(ix)-1d0)**2)*exp(-0.5d0*(xn(iy)+2d0)**2)
    end do
  end do

  do iy = 0,nx
    do ix = 0,nx
      xvec(ix,iy) = wfn(ix,iy) + wfn(iy,ix)
    end do
  end do
  ss = sum(xvec**2)*dx**2
  xvec = xvec/sqrt(ss)

  xx = sum(xvec**2)*dx**2
  call hpsi(xvec,hxvec)
  xhx = sum(xvec*hxvec)*dx**2
  lambda = xhx/xx

  res = sum(abs(hxvec-lambda*xvec)**2)*dx**2/xx

  do icg = 1,ncg

    gvec = 2d0*(hxvec - lambda*xvec)/xx
    gg0  = sum(gvec**2)*dx**2
    select case(icg)
    case(1)
      pvec = -gvec
    case default
      beta = gg0/gg
      pvec = -gvec + beta*pvec
    end select
    gg = gg0
          
    call hpsi(pvec, hpvec)
    pp  = sum(pvec**2)*dx**2
    php = sum(pvec*hpvec)*dx**2
    xp  = sum(pvec*xvec)*dx**2
    xhp = sum(hxvec*pvec)*dx**2
    
    aa=php*xp-xhp*pp
    bb=php*xx-xhx*pp
    cc=xhp*xx-xhx*xp
    ss=bb**2-4d0*aa*cc
    if(ss > 0d0)then
      alpha=(-bb+sqrt(ss))/(2d0*aa)
    else
      exit
    end if
    
    xvec = xvec + alpha*pvec
    call hpsi(xvec, hxvec)
    xx  = sum(xvec**2)*dx**2
    xhx = sum(xvec*hxvec)*dx**2
    lambda = xhx/xx

    res = sum(abs(hxvec-lambda*xvec)**2)*dx**2/xx
    if(res<res_epsilon)exit
    
  end do


! symmetrize
  wfn = xvec
  do iy = 0,nx
    do ix = 0,nx
      xvec(ix,iy) = wfn(ix,iy) + wfn(iy,ix)
    end do
  end do
  ss = sum(xvec**2)*dx**2
  xvec = xvec/sqrt(ss)


  call hpsi(xvec, hxvec)
  xx  = sum(xvec**2)*dx**2
  xhx = sum(xvec*hxvec)*dx**2
  lambda = xhx/xx
  res = sum(abs(hxvec-lambda*xvec)**2)*dx**2/xx
  xvec = xvec/xx
  write(*,*)"icg =",icg
  write(*,"(A,2x,2e26.16e3)")"lambda,res=",lambda,res

  wfn = xvec
  eps_gs = lambda

  do ix = 0,nx
    rho_gs(ix) = sum(wfn(:,ix)**2)*dx
  end do

  ss = sum(rho_gs)*dx
  rho_gs = rho_gs/sqrt(ss)

  open(20,file='rho_gs.out')
  do ix = 0,nx
    write(20,"(999e26.16e3)")xn(ix),rho_gs(ix)
  end do
  close(20)


end subroutine cg_gs
!-------------------------------------------------------------------------------
subroutine calc_hartree_fock
  use global_variables
  implicit none
  real(8) :: wfn_hf(0:nx), hwfn_1p_t(0:nx)
  real(8) :: wfn_ks(0:nx)
  real(8) :: v_hf(0:nx), rho_1p(0:nx)
  real(8) :: ss
  integer :: ix,iy
  real(8) :: eps, res
  integer :: iscf, nscf
  real(8) :: dt
  real(8) :: ss_se_hf, ss_se_ks, ss_hf_ks

  nscf = 2000
  dt = 0.01d0

  do ix = 0,nx
     wfn_hf(ix) = sqrt(rho_gs(ix))
  end do
  ss = sum(wfn_hf**2)*dx
  wfn_hf = wfn_hf/sqrt(ss)
  wfn_ks = wfn_hf

  rho_1p = wfn_hf**2
  v_hf = matmul(v_int, rho_1p)*dx + v_ext

  do iscf = 1, nscf
     call laplacian_1p(wfn_hf,hwfn_1p_t, -0.5d0)
     hwfn_1p_t = hwfn_1p_t + wfn_hf*v_hf

     eps = sum(wfn_hf*hwfn_1p_t)*dx
     res = sum((hwfn_1p_t - eps*wfn_hf)**2)*dx
     write(*,"(I7,2x,2e16.6e3)")iscf,eps,res

     wfn_hf  = wfn_hf  - dt*hwfn_1p_t


     ss = sum(wfn_hf**2)*dx
     wfn_hf = wfn_hf/sqrt(ss)

     rho_1p = wfn_hf**2
     v_hf = matmul(v_int, rho_1p)*dx + v_ext

  end do

  ss_se_hf = 0d0
  ss_se_ks = 0d0
  ss_hf_ks = 0d0

  do iy = 0,nx
     do ix = 0,nx

        ss_se_hf = ss_se_hf + wfn(ix,iy)*wfn_hf(ix)*wfn_hf(iy)
        ss_se_ks = ss_se_ks + wfn(ix,iy)*wfn_ks(ix)*wfn_ks(iy)
        ss_hf_ks = ss_hf_ks + wfn_hf(ix)*wfn_hf(iy)*wfn_ks(ix)*wfn_ks(iy)

     end do
  end do

  ss_se_hf = (ss_se_hf*dx**2)**2
  ss_se_ks = (ss_se_ks*dx**2)**2
  ss_hf_ks = (ss_hf_ks*dx**2)**2

  write(*,"(999e26.16e3)")lambda_int,ss_se_hf,ss_se_ks,ss_hf_ks


  

end subroutine calc_hartree_fock
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine laplacian_1p(wfn_in,Lwfn_out,fact)
  use global_variables
  real(8),intent(in) :: wfn_in(0:nx), fact
  real(8),intent(out) :: Lwfn_out(0:nx)
  real(8) :: wfn_tmp(-1:nx+1)
  real(8),parameter :: ct0 = -2d0, ct1 = 1d0
  real(8) :: c0,c1
  integer :: ix,iy

  c0 = fact*ct0/dx**2
  c1 = fact*ct1/dx**2

  wfn_tmp(-1) = 0d0
  wfn_tmp(0:nx) = wfn_in(0:nx)
  wfn_tmp(nx+1) = 0d0

    do ix = 0,nx
      Lwfn_out(ix) = c0*wfn_tmp(ix) +c1*( wfn_tmp(ix+1) + wfn_tmp(ix-1))
    end do

  end subroutine laplacian_1p
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
