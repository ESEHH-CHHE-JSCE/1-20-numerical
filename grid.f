      program gridconst
c
      implicit none
      include 'common.h'
c
      integer i,k
      integer iscf ! source flug 0:off,1:on
      real x(0:im,0:km),z(0:im,0:km)
      integer ifl(0:im,0:km) ! 0:open,1:block
      integer isu_u(0:im,0:km),isw_w(0:im,0:km),isr_c(0:im,0:km) ! source flag
c
      real u_l(0:im,0:km),w_l(0:im,0:km)
      real yp_l(0:im,0:km)
      real yr_l(0:im,0:km)
      real bc1_u(0:km),bc2_u(0:km) ! b.c.flow at U part on left & right
      real bc3_w(0:im),bc4_w(0:im) ! b.c.flow at W part on top & bottom
      real bc1_r(0:km),bc2_r(0:km) ! b.c.rho at center on left & right
      real bc3_r(0:im),bc4_r(0:im) ! b.c.rho at center on top & bottom
c
      real scu_u(0:im,0:km),scw_w(0:im,0:km),scr_c(0:im,0:km) ! source value
c
      real chleng,height
      real dx,dz
c
      real g,rho,rnu,diffyr
      real tuk,etime
      real lsor,soralp
      real pretime
c
c     /* reading basic condition data */
      open(11,file='condition.d',status='old')
      read(11,*) nx,nz
      read(11,*) iscf
      read(11,*) chleng,height
      read(11,*) g,rho,rnu,diffyr
      read(11,*) tuk,etime,dt
      read(11,*) lsor,soralp
      read(11,*) pretime
      close(11)
c
      dx=chleng/float(nx)
      dz=height/float(nz)
c
c     initialization
      do i=0,nx
        do k=0,nz
            x(i,k)=0.0+dx*float(i)
            z(i,k)=0.0+dz*float(k)
            yr_l(i,k)=rho
            yp_l(i,k)=0.0
            ifl(i,k)=0
        end do
      end do
      do i=0,nx+1
        do k=0,nz+1
            scu_u(i,k)=0.0
            scw_w(i,k)=0.0
            scr_c(i,k)=rho
            isu_u(i,k)=0
            isw_w(i,k)=0
            isr_c(i,k)=0
        end do
      end do
c
      do i=1,nx
          bc3_w(i)=0.0
          bc4_w(i)=0.0
      end do
c
      do k=1,nz
          bc1_u(k)=0.0
          bc2_u(k)=0.0
      end do
c
      do i=1,nx
          bc3_r(i)=rho
          bc4_r(i)=rho
      end do
      do k=1,nz
          bc1_r(k)=rho
          bc2_r(k)=rho
      end do
c
c     source settings
      do i=nx/2,nx/2
        do k=2,4
            isw_w(i,k)=1
            scw_w(i,k)=0.7
        end do
      end do
      do i=nx/2,nx/2
        do k=3,5
            isr_c(i,k)=1
            scr_c(i,k)=rho+5.0
        end do
      end do
c
      open(16,file='geometry.dat',status='unknown',form='formatted')
      write(16,*) nx,nz
      do i=0,nx
        do k=0,nz
c         write(16,'(2f12.6,I8)')
c     & x(i,k),z(i,k)
c     &,ifl(i,k)
         write(16,*) x(i,k),z(i,k),ifl(i,k)
        end do
      end do
      close(16)
c
      if(iscf.eq.1) then
       open(17,file='source.dat',status='unknown',form='formatted')
       write(17,*) nx,nz
       do i=0,nx
        do k=0,nz
         write(17,'(3f12.6,3I8)')
     &  scu_u(i,k),scw_w(i,k),scr_c(i,k)
     & ,isu_u(i,k),isw_w(i,k),isr_c(i,k)
        end do
       end do
       close(17)
      end if
c
c
      open(22,file='initial.dat',status='unknown',form='formatted')
      write(22,*) nx,nz
      do i=0,nx
        do k=0,nz
         write(22,*) u_l(i,k),w_l(i,k),yp_l(i,k),yr_l(i,k) ! at node
        end do
      end do
      close(22)
c
      open(23,file='bc_vert.dat',status='unknown',form='formatted')
      write(23,*) nz
      do k=1,nz
         write(23,*) bc1_u(k),bc2_u(k),bc1_r(k),bc2_r(k) ! at u part
      end do
      close(23)
c 
      open(24,file='bc_horz.dat',status='unknown',form='formatted')
      write(24,*) nx
      do i=1,nx
         write(24,*) bc3_w(i),bc4_w(i),bc3_r(i),bc4_r(i) ! at w part
      end do
      close(24)
c
      end