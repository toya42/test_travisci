module constants
   use,intrinsic :: iso_fortran_env
   implicit none
!   real(real128),parameter :: PI = 4.0q0*atan(1q0)
   real(real64),parameter :: PI = 4.0d0*atan(1d0)
   real(real64),parameter :: NU = 2.0d-6
!   complex(8),parameter :: IU = (0.0d0,1.0d0)

!! IM>=3*KM+1, JM>=3*LM+1
!   integer(int32),parameter :: KM = 341
!   integer(int32),parameter :: LM = 341
   integer(int32),parameter :: IM = 4*2**6
   integer(int32),parameter :: JM = 4*2**6
!   integer(int32),parameter :: KL = (2*KM+1)*(2*LM+1)
   real(real64),parameter :: DT = 2.0d-2
   real(real64),parameter :: TM = 200.0d0
!   integer(int32),parameter :: NMAX = 10000   ! TM/DT
!   integer(int32),parameter :: NOUT = 50
   integer(int32),parameter :: NMAX = 10   ! TM/DT
   integer(int32),parameter :: NOUT = 10
end module constants
!===
module arrays
   use,intrinsic :: iso_fortran_env
   use constants
   implicit none
   real(real64),dimension(2,0:IM/2-1,-JM/2:JM/2-1) :: aij  ! spectral coefficients of vorticity
end module arrays
!===
module fft_mkl
   use,intrinsic :: iso_fortran_env
   use mkl_dfti
   implicit none
   private
   public fft_initialize,fft_finalize,fft_execute_forward,fft_execute_backward
   contains
   subroutine fft_initialize(imax,jmax,des_r2c,des_c2r)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      type(dfti_descriptor),pointer :: des_r2c,des_c2r
      integer(int32),dimension(12) :: dftistatus
      integer(int32),dimension(2) :: length
      integer(int32),dimension(3) :: stride_i,stride_o
      integer(int32) :: i
      real(real64) :: ni

      ni = 1.0d0/dble(imax*jmax)
      length(1) = imax
      length(2) = jmax
      stride_i(1) = 0
      stride_i(2) = 1
      stride_i(3) = imax
      stride_o(1) = 0
      stride_o(2) = 1
      stride_o(3) = imax/2+1

!r2c
      dftistatus( 1) = dfticreatedescriptor(des_r2c,dfti_double,dfti_real,2,length)
      dftistatus( 2) = dftisetvalue        (des_r2c,dfti_forward_scale,ni)
      dftistatus( 3) = dftisetvalue        (des_r2c,dfti_placement,dfti_not_inplace)
      dftistatus( 4) = dftisetvalue        (des_r2c,dfti_conjugate_even_storage,dfti_complex_complex)
      dftistatus( 5) = dftisetvalue        (des_r2c,dfti_packed_format,dfti_cce_format)
      dftistatus( 6) = dftisetvalue        (des_r2c,dfti_input_strides, stride_i)
      dftistatus( 7) = dftisetvalue        (des_r2c,dfti_output_strides,stride_o)
      dftistatus( 8) = dfticommitdescriptor(des_r2c)
!c2r
      dftistatus( 9) = dfticopydescriptor  (des_r2c,des_c2r)
      dftistatus(10) = dftisetvalue        (des_c2r,dfti_input_strides, stride_o)
      dftistatus(11) = dftisetvalue        (des_c2r,dfti_output_strides,stride_i)
      dftistatus(12) = dfticommitdescriptor(des_c2r)

#ifdef debug
      print *,'fft_initialize) IMAX   =',imax
      print *,'fft_initialize) JMAX   =',jmax
      do i=1,12
         if(dftistatus(i).ne.0) then
            print *,'fft_initialize) dft setting error:',i,dftistatus(i)
            stop
         end if
      end do
      print *,'fft_initialize) dft settings are completed'
#endif
  

      return
   end subroutine fft_initialize
   subroutine fft_finalize(des_r2c,des_c2r)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_r2c,des_c2r
      integer(int32),dimension(2) :: dftistatus
      integer(int32) :: i

      dftistatus(1)=DftiFreeDescriptor(Des_r2c)
      dftistatus(2)=DftiFreeDescriptor(Des_c2r)
!#ifdef debug
      do i=1,2
         if(dftistatus(i).ne.0) then
            print *,'fft_finalize) dft finalization error:',i
            stop
         end if
      end do
      print *,'fft_finalize) dft finalization is completed'
!#endif
      return
   end subroutine fft_finalize
   subroutine fft_execute_forward(des_r2c,length,real2d,complex2d)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_r2c
      integer(int32),dimension(2),intent(in) :: length
      real(real64),dimension(length(1)),intent(in)  :: real2d
      real(real64),dimension(length(2)),intent(out) :: complex2d
      integer(int32) :: dftistatus

      dftistatus = dfticomputeforward(des_r2c,real2d,complex2d)
#ifdef debug
      if(dftistatus.ne.0) then
         print *,'fft_execute_forward) dft finalization error:',dftistatus
         stop
      end if
#endif
      return
   end subroutine fft_execute_forward
   subroutine fft_execute_backward(des_c2r,length,complex2d,real2d)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_c2r
      integer(int32),dimension(2),intent(in) :: length
      real(real64),dimension(length(2)),intent(in)  :: complex2d
      real(real64),dimension(length(1)),intent(out) :: real2d
      integer(int32) :: dftistatus

      dftistatus = dfticomputebackward(des_c2r,complex2d,real2d)
#ifdef debug
      if(dftistatus.ne.0) then
         print *,'fft_execute_forward) dft finalization error:',dftistatus
         stop
      end if
#endif
      return
   end subroutine fft_execute_backward
end module fft_mkl

!===
!module 2dfft
!   use,intrinsic :: iso_fortran_env
!   use fft_mkl,only : fft_execute_forward,fft_execute_backward
!   use constants,only : IU,PI
!   implicit none
!   private
!   public physical_to_transform,transform_to_physical
!   contains
!! physical space(physical values) to transform space(spectral coefficients)
!   subroutine physical_to_transform(imax,jmax,kmax,lmax,zeta,aij)
!      use,intrinsic :: iso_fortran_env
!      implicit none
!      integer(int32),intent(in) :: imax,jmax,kmax,lmax
!      real(real64),dimension(0:imax-1,0:jmax-1),  intent(in)  :: zeta
!      complex(8),dimension(0:kmax,-lmax:lmax),intent(out) :: aij
!      complex(8),dimension(0:imax/2,0:jmax-1) :: c
!      complex(8),dimension(0:jmax-1) :: x
!      complex(8),dimension(0:jmax-1,0:kmax) :: g,d
!      integer(int32) :: i,j,k,l
!
!      aij = 0.0d0
!      c   = 0.0d0
!      x   = 0.0d0
!      g   = 0.0d0
!      d   = 0.0d0
!
!      do j=0,jmax-1
!!         call fft1d_r_forward(imax,zeta(0,j),itr,c(0,j),wci,wrif,y(0,j))
!         call fft1d_r_forward(imax,zeta(0,j),itr,c(0,j),wci,wrif)
!!         c(0,j) = real(c(0,j))
!      end do
!      
!!      do i=0,imax/2-1
!!         do j=0,jmax-1
!!            write(30,*) j+1,i,real(c(i,j)),imag(c(i,j))
!!!            write(51,*) j+1,i,real(y(i,j)),imag(y(i,j))
!!!            write(52,*) j+1,i,real(c(i,j)),imag(c(i,j))
!!         end do
!!      end do
!!
!!
!!      do k=0,kmax
!!         do j=0,jmax-1
!!            write(21,*) j,k,real(c(k,j)),imag(c(k,j))
!!         end do
!!      end do
!
!!      do j=0,jmax-1
!!         do i=0,kmax
!!            d(j,i) = c(i,j)
!!         end do
!!      end do
!      d = transpose(c)
!
!
!      do k=0,kmax
!         call fft1d_c_forward(jmax,d(0,k),itc,g(0,k),wcj)
!      end do
!
!!      do k=0,kmax
!!         do j=0,jmax-1
!!            write(22,*) j,k,real(g(j,k)),imag(g(j,k))
!!         end do
!!      end do
!
!      do k=1,kmax
!         do l=1,lmax
!            aij(k, l) = g(     l,k)
!            aij(k,-l) = g(jmax-l,k)
!         end do
!         aij(k,0) = g(0,k)
!      end do
!
!      do l=1,lmax
!         aij(0, l) = g(l,0)
!         aij(0,-l) = conjg(aij(0,l))
!!         if(aij(0,l).ne.aij(0,-l)) then
!!            write(*,*) l
!!            write(*,*) aij(0,l)
!!            write(*,*) aij(0,-l)
!!         end if
!      end do
!
!      aij(0,0) = g(0,0)
!
!!      do l=-lmax,lmax
!!         do k=0,kmax
!!            write(10,*) k,l,real(aij(k,l)),imag(aij(k,l))
!!         end do
!!      end do
!
!
!      return
!   end subroutine physical_to_transform
!
!! transform space(spectral coefficients) to physical space(physical values)
!   subroutine transformto_physical(imax,jmax,kmax,lmax,zeta,aij)
!      use,intrinsic :: iso_fortran_env
!      implicit none
!      integer(int32),intent(in) :: imax,jmax,kmax,lmax
!      real(real64),dimension(0:imax-1,0:jmax-1),        intent(out) :: zeta
!      complex(8),dimension(0:kmax,-lmax:lmax),intent(in)  :: aij
!!      complex(8),dimension(0:jmax-1) :: a,x,a2
!!      real(real64),   dimension(0:jmax-1) :: x2
!!      complex(8),dimension(0:kmax,  0:jmax-1) :: phitkn
!!      complex(8),dimension(0:imax-1,0:jmax-1) :: phimn
!      complex(8),dimension(0:jmax-1,0:imax/2) :: g,d
!      complex(8),dimension(0:imax/2,0:jmax-1) :: c
!      integer(int32) :: i,j,k,l
!
!      zeta=0.0d0
!      g = 0.0d0
!      d = 0.0d0
!      c = 0.0d0
!
!      do k=1,kmax
!         do l=1,lmax
!            g(     l,k) = aij(k, l)
!            g(jmax-l,k) = aij(k,-l)
!         end do
!!         do l=lmax+1,jmax-lmax-1
!!            g(l,k) = 0.0d0
!!         end do
!      end do
!
!      do k=1,kmax
!         g(0,k) = aij(k,0)
!      end do
!      
!      do l=1,lmax
!         g(     l,0) = aij(0,l)
!         g(jmax-l,0) = conjg(aij(0,l))
!!         g(jmax-l,0) = conjg(aij(0,-l))
!!         g(jmax-l,0) = (aij(0,-l))
!!         g(jmax-l,0) = (aij(0,l))
!      end do
!!      do l=lmax+1,jmax-lmax-1
!!         g(l,0) = 0.0d0
!!      end do
!
!      g(0,0) = real(aij(0,0))
!
!!      do i=kmax+1,imax/2-1
!!         do j=0,jmax-1
!!            g(j,i) = 0.0d0
!!         end do
!!      end do
!      
!      do k=0,kmax
!         call fft1d_c_backward(jmax,g(0,k),itc,d(0,k),wcj)
!      end do
!
!      c = transpose(d)
!      do j=0,jmax-1
!!!         do i=0,kmax
!!         do i=0,imax/2
!!            c(i,j) = d(j,i)
!!         end do
!         do i=kmax+1,imax/2
!            c(i,j) = 0.0d0
!         end do
!      end do
!
!      do j=0,jmax-1
!         call fft1d_r_backward(imax,c(0,j),itr,zeta(0,j),wci,wrib)
!      end do
!
!      return
!   end subroutine transformto_physical
!
!end module 2dfft
!!===
module flowfield
   use,intrinsic :: iso_fortran_env
   use constants, only : PI
   use fft_mkl,only : fft_execute_forward,fft_execute_backward
   implicit none
   private
   public flowfield_initialize,wrtd
   public aij_to_workc,workc_to_aij
!   public n_length
   real(real64),dimension(:),allocatable :: xi,yj
!org   real(real64),dimension(:,:),allocatable :: zeta_org
!   real(real64),dimension(:),allocatable :: workc
   integer(int32),dimension(2) :: n_length
   contains
   subroutine flowfield_initialize(imax,jmax,aij,des_n_r2c)
!   subroutine flowfield_initialize(imax,jmax,aij,des_n_r2c,des_n_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(out) :: aij
      type(dfti_descriptor),pointer :: des_n_r2c
!org      type(dfti_descriptor),pointer :: des_n_c2r
      real(real64),dimension(imax,jmax) :: zeta
      real(real64),dimension(imax+2,jmax) :: workc
      integer(int32) :: i,j
      real(real64) :: x,x1,x2,y,y1,y2,sigma,sgm2i,c
!org      real(real64) :: error
      
!      print *,'flowfield_initialize) #1'
      allocate(xi(imax))
      allocate(yj(jmax))
      n_length(1) = imax*jmax
      n_length(2) = (imax/2+1)*2*jmax
!      allocate(workc(n_length(2)))

!org      allocate(zeta_org(imax,jmax))

      do i=1,imax
         xi(i)  = 2.00d0*PI*dble(i-1)/dble(imax)
      end do
      do j=1,jmax
         yj(j)  = 2.00d0*PI*dble(j-1)/dble(jmax)
      end do
      do j=1,jmax
         do i=1,imax
            x  = xi(i)
            y  = yj(j)
            x1 = 0.80d0*PI
            y1 = x1
            x2 = 1.20d0*PI
            y2 = x2
            sigma = 0.10d0*PI
            sgm2i = 1.0d0/(sigma*sigma)
            c  = 0.0322447d0
            zeta(i,j) = dexp((dcos(x-x1)+dcos(y-y1)-2.0d0)*sgm2i)  &
&                      +dexp((dcos(x-x2)+dcos(y-y2)-2.0d0)*sgm2i)-c
!org            zeta_org(i,j) = zeta(i,j)
         end do
      end do


!      print *,'flowfield_initialize) #2'
      call fft_execute_forward(des_n_r2c,n_length,zeta,workc)
      aij = workc_to_aij(imax,jmax,workc)

!org      call fft_execute_backward(des_n_c2r,n_length,workc,zeta)
!org      error = 0.0d0
!org      do j=1,jmax
!org         do i=1,imax
!org            error = error+dabs(zeta(i,j)-zeta_org(i,j))
!org         end do
!org      end do
!org      write(*,*) 'wrtd) error:', error/dble(imax*jmax)



!      print *,'flowfield_initialize) #3'

      return
   end subroutine flowfield_initialize

   subroutine wrtd(imax,jmax,loop,aij,des_n_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      implicit none
      integer(int32),intent(in) :: imax,jmax
      integer(int32),intent(in) :: loop
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
      type(dfti_descriptor),pointer :: des_n_c2r
      real(real64),dimension(imax,jmax) :: zeta

      integer(int32) :: i,j
      real(real64),dimension(imax+2,jmax) :: workc
      real(real64),dimension(:,:,:,:),allocatable :: buf
      character(len=50) :: filename
!org      real(real64) :: error

      if(loop.eq.0) then
         allocate(buf(imax,jmax,1,3))
         do j=1,jmax
            do i=1,imax
               buf(i,j,1,1) = xi(i)
               buf(i,j,1,2) = yj(j)
               buf(i,j,1,3) = 0.0d0
            end do
         end do
         open(10,file='output/grid.xyz',form='unformatted',access='stream',status='replace')
!         open(10,file='grid.xyz',form='unformatted',status='replace')
         i=imax
         j=jmax
         write(10) i,j,1
         write(10) buf
         close(10)
         deallocate(buf)
         deallocate(xi,yj)
      end if

!      allocate(workc(n_length(2)))
      workc = aij_to_workc(imax,jmax,aij)
      call fft_execute_backward(des_n_c2r,n_length,workc,zeta)

!org      error = 0.0d0
!org      do j=1,jmax
!org         do i=1,imax
!org            error = error+dabs(zeta(i,j)-zeta_org(i,j))
!org         end do
!org      end do
!org      write(*,*) 'wrtd) error:', error/dble(imax*jmax)
!org      deallocate(zeta_org)
!org!      stop


!      allocate(buf(imax,jmax,1,1))
!      do j=1,jmax
!         do i=1,imax
!            buf(i,j,1,1)=zeta(i,j)
!         end do
!      end do
!      write(filename,'("output/zeta_",i5.5,".fun")') loop
!      open(10,file=filename,form='unformatted',access='stream',status='replace')
!      i=imax
!      j=jmax
!      write(10) i,j,1,1
!      write(10) buf
!      close(10)
!      deallocate(buf)

      allocate(buf(imax,jmax,1,5))
      do j=1,jmax
         do i=1,imax
            buf(i,j,1,1)=zeta(i,j)
            buf(i,j,1,2)=0.0d0
            buf(i,j,1,3)=0.0d0
            buf(i,j,1,4)=0.0d0
            buf(i,j,1,5)=zeta(i,j)
         end do
      end do
      write(filename,'("output/zeta_",i5.5,".q")') loop
      open(11,file=filename,form='unformatted',access='stream',status='replace')
      i=imax
      j=jmax
      write(11) i,j,1
      write(11) buf
      close(11)
      deallocate(buf)


      write(6,*) 'wrtd) flowfield output,loop=',loop

      return
   end subroutine wrtd
   
   pure function workc_to_aij(imax,jmax,workc) result(aij)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),dimension(imax+2,jmax),intent(in) :: workc
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: aij
      integer(int32) :: i,j,i2

      do j=0,jmax/2-1
         do i=0,imax/2-1
            i2 = 2*i
            aij(1,i,j) = workc(i2+1,j+1)
            aij(2,i,j) = workc(i2+2,j+1)
         end do
      end do
      do j=jmax/2,jmax-1
         do i=0,imax/2-1
            i2 = 2*i
!            aij(1,i,jmax/2-1-j) = workc(i2+1,j+1)
!            aij(2,i,jmax/2-1-j) = workc(i2+2,j+1)
            aij(1,i,-jmax+j) = workc(i2+1,j+1) 
            aij(2,i,-jmax+j) = workc(i2+2,j+1) 

         end do
      end do


!      do j=0,jmax/2-1
!         do i=0,imax/2-1
!            aij(1,i,j) = workc(2*i+1,j+1)
!            aij(2,i,j) = workc(2*i+2,j+1)
!         end do
!      end do
!      do j=-jmax/2,-1
!         do i=0,imax/2-1
!            aij(1,i,j) = workc(2*i+1,jmax/2-j)
!            aij(2,i,j) = workc(2*i+2,jmax/2-j)
!         end do
!      end do
      
   end function workc_to_aij
   pure function aij_to_workc(imax,jmax,aij) result(workc)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
      real(real64),dimension(imax+2,jmax) :: workc
      integer(int32) :: i,j,i2

      do j=0,jmax/2-1
         do i=0,imax/2-1
            i2 = 2*i
            workc(i2+1,j+1) = aij(1,i,j)
            workc(i2+2,j+1) = aij(2,i,j)
         end do
         workc(imax+1,j+1) = 0.0d0
         workc(imax+2,j+1) = 0.0d0
      end do
      do j=jmax/2,jmax-1
         do i=0,imax/2-1
            i2 = 2*i
!            workc(i2+1,j+1) = aij(1,i,jmax/2-1-j)
!            workc(i2+2,j+1) = aij(2,i,jmax/2-1-j)
            workc(i2+1,j+1) = aij(1,i,-jmax+j)
            workc(i2+2,j+1) = aij(2,i,-jmax+j)

         end do
         workc(imax+1,j+1) = 0.0d0
         workc(imax+2,j+1) = 0.0d0
      end do


!!      workc = 0.0d0
!      do j=0,jmax/2-1
!         do i=0,imax/2-1
!            workc(2*i+1,j+1) = aij(1,i,j)
!            workc(2*i+2,j+1) = aij(2,i,j)
!         end do
!         workc(imax+1,j+1) = 0.0d0
!         workc(imax+2,j+1) = 0.0d0
!      end do
!      do j=-jmax/2,-1
!         jj = jmax/2-1-j
!         do i=0,imax/2-1
!            workc(2*i+1,jj) = aij(1,i,j)
!            workc(2*i+2,jj) = aij(2,i,j)
!         end do
!         workc(imax+1,jj) = 0.0d0
!         workc(imax+2,jj) = 0.0d0
!      end do
      
   end function aij_to_workc


end module flowfield
!!===
module time_integration
   use,intrinsic :: iso_fortran_env
   use constants, only : NU
   implicit none
   private
   public timeintegration_initialize,rk4n_al
   real(real64),dimension(:,:),allocatable :: nuij,nlij,ztij
   integer(int32),dimension(2) :: p_length
   contains
   subroutine timeintegration_initialize(imax,jmax,dt)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),intent(in) :: dt
      integer(int32) :: i,j,ip,jp

      ip = imax*3/2
      jp = jmax*3/2
      p_length(1) = ip*jp
      p_length(2) = (ip/2+1)*2*jp   ! (ip+2)*jp

      allocate(nuij(0:imax/2-1,-jmax/2:jmax/2-1))
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            nuij(i,j) = dexp(-0.50d0*dt*NU*dble(i*i+j*j))
         end do
      end do
      
      allocate(nlij(0:imax/2-1,-jmax/2:jmax/2-1))
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               nlij(i,j) = -dble(i*i+j*j)
         end do
      end do
      nlij(0,0) = 1.0d0
      nlij = 1.0d0/nlij

      allocate(ztij(0:imax/2-1,-jmax/2:jmax/2-1))
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               ztij(i,j) = dble(i*i+j*j)
         end do
      end do


      return
   end subroutine timeintegration_initialize

   subroutine rk4n_al(imax,jmax,dt,w0,des_p_r2c,des_p_c2r)
!   subroutine rk4n_al(imax,jmax,dt,w0,des_p_r2c,des_p_c2r,des_n_r2c,des_n_c2r)
! time integration
!     RK4 for Non-linear term
!     Analytical solution for Linear term
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),intent(in) :: dt
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
!      type(dfti_descriptor),pointer :: des_n_r2c,des_n_c2r
      real(real64),parameter :: c16 = 1.0d0/6.0d0
      real(real64),parameter :: c13 = 1.0d0/3.0d0
      real(real64),dimension(imax*jmax),intent(inout) :: w0
      real(real64),dimension(imax*jmax) :: w1,w2,w3

      call intg_nonlinear(imax,jmax,w0,w1,des_p_r2c,des_p_c2r)
!      call intg_nonlinear(imax,jmax,w0,w1,des_p_r2c,des_p_c2r,des_n_r2c,des_n_c2r)

      w2 = w0+0.5d0*dt*w1
      w3 = w0+  c16*dt*w1

      call intg_linear(imax,jmax,w0)
      call intg_linear(imax,jmax,w2)
      call intg_linear(imax,jmax,w3)
      call intg_nonlinear(imax,jmax,w2,w1,des_p_r2c,des_p_c2r)
!      call intg_nonlinear(imax,jmax,w2,w1,des_p_r2c,des_p_c2r,des_n_r2c,des_n_c2r)
      
      w2 = w0+0.5d0*dt*w1
      w3 = w3+  c13*dt*w1

      call intg_nonlinear(imax,jmax,w2,w1,des_p_r2c,des_p_c2r)
!      call intg_nonlinear(imax,jmax,w2,w1,des_p_r2c,des_p_c2r,des_n_r2c,des_n_c2r)

      w2 = w0+      dt*w1
      w3 = w3+  c13*dt*w1

      call intg_linear(imax,jmax,w2)
      call intg_linear(imax,jmax,w3)
      call intg_nonlinear(imax,jmax,w2,w1,des_p_r2c,des_p_c2r)
!      call intg_nonlinear(imax,jmax,w2,w1,des_p_r2c,des_p_c2r,des_n_r2c,des_n_c2r)

      w0 = w3+  c16*dt*w1

      return
   end subroutine rk4n_al
   
   subroutine intg_nonlinear(imax,jmax,w,dw,des_p_r2c,des_p_c2r)
!   subroutine intg_nonlinear(imax,jmax,w,dw,des_p_r2c,des_p_c2r,des_n_r2c,des_n_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      use fft_mkl, only : fft_execute_forward,fft_execute_backward
      use flowfield, only : aij_to_workc,workc_to_aij!     ,n_length
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: w
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(out) :: dw
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
!      type(dfti_descriptor),pointer :: des_n_r2c,des_n_c2r
      real(real64),dimension(2,0:imax*3/4-1,-jmax*3/4:jmax*3/4-1) :: ht_p  ! 0 padding
!      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: ht    
!      real(real64),dimension(imax+2,jmax) :: work
      real(real64),dimension(imax*3/2+2,jmax*3/2) :: workc
      real(real64)   ,dimension(0:imax*3/2-1,0:jmax*3/2-1)   :: u,v,dzdx,dzdy,f_p
!      real(real64)   ,dimension(0:imax-1,0:jmax-1)   :: un
      integer(int32) :: i,j,ip,jp
      real(real64) :: di,dj

      ip = imax*3/2
      jp = jmax*3/2

      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            dw(1,i,j) = w(1,i,j)*nlij(i,j)
            dw(2,i,j) = w(2,i,j)*nlij(i,j)
         end do
      end do

! calc. u
      ht_p = 0.0d0
      do j=-jmax/2,jmax/2-1
         dj = dble(j)
         do i=0,imax/2-1
            ht_p(1,i,j) =  dj*dw(2,i,j)
            ht_p(2,i,j) = -dj*dw(1,i,j)
!            ht(1,i,j) =  dj*dw(2,i,j)
!            ht(2,i,j) = -dj*dw(1,i,j)
         end do
      end do
      workc = aij_to_workc(ip,jp,ht_p)
!      work = aij_to_workc(imax,jmax,ht)
!      call fft_execute_backward(des_n_c2r,n_length,work,un)
      call fft_execute_backward(des_p_c2r,p_length,workc,u)
!      do j=0,jp-1
!         do i=0,ip-1
!            write(12,*) i,j,u(i,j)
!         end do
!         write(12,*)
!      end do
!      do j=0,jmax-1
!         do i=0,imax-1
!            write(13,*) i,j,un(i,j)
!         end do
!         write(13,*)
!      end do
!
!      stop
! calc. v
      ht_p = 0.0d0
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            di = dble(i)
            ht_p(1,i,j) = -di*dw(2,i,j)
            ht_p(2,i,j) =  di*dw(1,i,j)
         end do
      end do
      workc = aij_to_workc(ip,jp,ht_p)
      call fft_execute_backward(des_p_c2r,p_length,workc,v)

! calc. dzdx
      ht_p = 0.0d0
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            di = dble(i)
            ht_p(1,i,j) =  di*ztij(i,j)*dw(2,i,j)
            ht_p(2,i,j) = -di*ztij(i,j)*dw(1,i,j)
         end do
      end do
      workc = aij_to_workc(ip,jp,ht_p)
      call fft_execute_backward(des_p_c2r,p_length,workc,dzdx)
! calc. dzdy
      ht_p = 0.0d0
      do j=-jmax/2,jmax/2-1
         dj = dble(j)
         do i=0,imax/2-1
            ht_p(1,i,j) =  dj*ztij(i,j)*dw(2,i,j)
            ht_p(2,i,j) = -dj*ztij(i,j)*dw(1,i,j)
         end do
      end do
      workc = aij_to_workc(ip,jp,ht_p)
      call fft_execute_backward(des_p_c2r,p_length,workc,dzdy)

! calc. Fij
      f_p = u*dzdx+v*dzdy

      call fft_execute_forward(des_p_r2c,p_length,f_p,workc)
      ht_p = workc_to_aij(ip,jp,workc)

       do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            dw(1,i,j) = -ht_p(1,i,j)
            dw(2,i,j) = -ht_p(2,i,j)
         end do
      end do
     
!      do l=-lmax,lmax
!         do k=0,kmax
!            dw(k,l) = dble(k*k-l*l)*uvht(k,l)+dble(k*l)*v2u2ht(k,l)
!!            if(k.eq.0) then
!!               write(*,*) k,l,dw(k,l)
!!            end if
!!            write(42,*) k,l,real(dw(k,l)),imag(dw(k,l))
!!            write(44,*) k,l,real(uvht(k,l)),imag(uvht(k,l))
!!            write(45,*) k,l,real(v2u2ht(k,l)),imag(v2u2ht(k,l))
!         end do
!      end do
!!      stop
      return

   end subroutine intg_nonlinear

   subroutine intg_linear(imax,jmax,w)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(inout) :: w
      integer(int32) :: i,j

      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            w(1,i,j) = w(1,i,j)*nuij(i,j)
            w(2,i,j) = w(2,i,j)*nuij(i,j)
         end do
      end do

!      do l=-lmax,lmax
!         do k=0,kmax
!            write(26,*) k,l,real(w(k,l)),imag(w(k,l))
!         end do
!      end do
!      stop

      return

   end subroutine intg_linear

end module time_integration
!===
program main
   use,intrinsic :: iso_fortran_env
!$   use omp_lib
   use constants
   use arrays
   use mkl_dfti
   use fft_mkl, only : fft_initialize,fft_finalize
   use flowfield, only : flowfield_initialize,wrtd
   use time_integration, only : timeintegration_initialize,rk4n_al
   implicit none
   type(dfti_descriptor),pointer :: des_n_r2c, des_n_c2r ! normal    (length=N)
   type(dfti_descriptor),pointer :: des_p_r2c, des_p_c2r ! padding   (length=N*3/2)
   integer(int32) :: n
   real(real64) :: time
!$   real(real64) :: t1,t2,t3

! pre-process
   print *,'main) initialization start'
   print *,'main) --fft--'
   print *,'main) for N to N transform'
   call fft_initialize(IM,JM,des_n_r2c,des_n_c2r)
   print *,'main) for N*3/2 to N*3/2 transform'
   call fft_initialize(IM*3/2,JM*3/2,des_p_r2c,des_p_c2r)
!  initial flow field setting. zeta(vorticity) -> aij(spectral coefficients)
   print *,'main) --flow field--'
   call flowfield_initialize(IM,JM,aij,des_n_r2c)
!   call flowfield_initialize(IM,JM,aij,des_n_r2c,des_n_c2r)
!  linear term (viscous term) calculation
   print *,'main) --linear term--'
   call timeintegration_initialize(IM,JM,DT)
   print *,'main) initialization end'
!
   print *,'main) output flow field at t=0'
   call wrtd(IM,JM,0,aij,des_n_c2r)

! main sequence
!$   t3 = 0.0d0
   do n=1,NMAX
      time = dble(n-1)*DT
#ifdef debug
      print *,'loop',n
#endif
!      call rk4n_al(IM,JM,DT,aij,des_p_r2c,des_p_c2r,des_n_r2c,des_n_c2r)
!$    t1 = omp_get_wtime()      
      call rk4n_al(IM,JM,DT,aij,des_p_r2c,des_p_c2r)
!$    t2 = omp_get_wtime()      
!$    t3 = t3+t2-t1
      if(mod(n,NOUT).eq.0) then
         call wrtd(IM,JM,n,aij,des_n_c2r)
      end if
   end do

! post-process
!$   print *,'time',t3
   call fft_finalize(des_n_r2c,des_n_c2r)
   call fft_finalize(des_p_r2c,des_p_c2r)

   stop
end program main
