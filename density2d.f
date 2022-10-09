c
c     **************************************************
c     *                                                *
c     *       2D Density Flow with  CIP method         *
c     *                                                *
c     *            Ryosuke Akahori Sep/19/22           *
c     *                                                *
c     **************************************************
c
      program density2d ! vertical 2d density flow
c
c
c     ////////  Descripttion of variables  ////////////////////////////////////////
c     /
c     /  i,k,l,m,n : loop counter
c     /     
c     /  im: maximum number of dimension on x or xi axis direction 
c     /  km: maximum number of dimension on z or zeta axis direction
c     /  time: calculation time
c     /
c     /  itt  : temporary time step in each output interval
c     /  itout: number of calculation timestep for output
c     /  icont: output step number
c     /
c     /  nx : number of grids on xi direction
c     /  nz : number of grids on zeta direction
c     /
c     /  rnu : coefficient of kinematic viscosity
c     /  diffyr : coefficient of raito of density diffusion to rnu_e
c     /  tuk : interval time for output
c     /  etime : end time of calculation
c     /  dt  : calculation period for one time step
c     /  lsor: iteration number for sor method
c     /  soralp: coefficient for sor acceleration
c     /
c     /  g   : gravity accelelation  
c     /
c     /  x   : phisical grid location of x at corner of cubic cell
c     /  z   : phisical grid location of z at corner of cubic cell
c     /  x_c : phisical grid location of x at center of cell
c     /  z_c : phisical grid location of z at center of cell
c     /  xi  : grid position of xi(0 to 1) at center of cell
c     /  zeta : grid position of zeta(0 to 1) at center of cell
c     /  dxi : grid cell size on xi deiraction
c     /  dzet: grid cell size on zeta deiraction
c     /
c     /  x_xi : transformation coefficient of dx/dxi
c     /  z_zet: transformation coefficient of dz/dzet
c     /  xi_x : transformation coefficient of dx/dxi
c     /  zet_z: transformation coefficient of dy/dzet
c     /  rj   : Jacobian
c     /
c     /  u_c : velocity of flow on x-axis(decalto) at the center(CP) of cell
c     /  w_c : velocity of flow on z-axis(decalto) at the center(CP) of cell
c     /  ub  : velocity of flow on xi-axis(bfc) at u-position(UP) of cell
c     /  wb  : velocity of flow on zeta-axis(bfc) at w-position(WP) of cell
c     /  ubn : temporary velocity of flow on xi-axis(bfc) at UP
c     /  wbn : temporary velocity of flow on zeta-axis(bfc) at WP
c     /  u_l : velocity of flow on x-axis(decalto) at node
c     /  w_l : velocity of flow on z-axis(decalto) at node
c     /
c     /  gux : xi-direction gradient of ubn at UP 
c     /  guz : zeta-direction gradient of ubn at UP
c     /  gwx : xi-direction gradient of wbn at WP
c     /  gwz : zeta-direction gradient of wbn at WP
c     /
c     /  yp  : fluctuation of pressure from hydro static pressure	at CP
c     /  ypn : temporary fluctuation of pressure from hydro static pressure at CP
c     /
c     /  rnu_e : coefficient of kinematic viscosity for turbulent model
c     /
c     /  rho : density of water (kg/m^3)
c     /  yr: density of local water (kg/m^3) at CP
c     /  yrn: temporary density of local water (kg/m^3) at CP
c     /  
c     //////////////////////////////////////////////////////////////////////////////
c
c
      implicit none
c
c *** including fundamental valiables  ***
c ***    (im, km, nx, nz, dt, time)    ***
      include 'common.h'
c
      real g
      integer i,k,lsor,m
      integer itt,itout,icount
      real rnu,diffyr
      real tuk,etime,soralp
      real dxi,dzet
c
      real x(0:im,0:km),z(0:im,0:km)
      real x_c(0:im,0:km),z_c(0:im,0:km)
      real xi(0:im),zet(0:km)
      real x_xi(0:im,0:km),z_zet(0:im,0:km)
      real xi_x(0:im,0:km),zet_z(0:im,0:km) 
      real dxa(0:im,0:km),dza(0:im,0:km)
c
      real rj(0:im,0:km)
c
      real u_c(0:im,0:km),w_c(0:im,0:km)
      real ub_l(0:im,0:km),wb_l(0:im,0:km)
      real u_l(0:im,0:km),w_l(0:im,0:km)
      real ub(0:im,0:km),wb(0:im,0:km)
     1    ,ubn(0:im,0:km),wbn(0:im,0:km)
     2    ,ubt(0:im,0:km),wbt(0:im,0:km)
c
      real ub_w(0:im,0:km),wb_u(0:im,0:km)
c
      real yp(0:im,0:km),ypn(0:im,0:km)
      real ypn3(0:im,0:km)
      real yp_l(0:im,0:km)
c
      real yr(0:im,0:km),yrn(0:im,0:km),yrt(0:im,0:km)
      real yr_l(0:im,0:km)
c
      real gux(0:im,0:km),guz(0:im,0:km)
      real gwx(0:im,0:km),gwz(0:im,0:km)
      real grx(0:im,0:km),grz(0:im,0:km)
c
      real bc1_u(0:km),bc2_u(0:km)
      real bc3_w(0:im),bc4_w(0:im)
c
      real bc1_r(0:km),bc2_r(0:km)
      real bc3_r(0:im),bc4_r(0:im)
c
      real bc1_ub(0:km),bc2_ub(0:km)
      real bc3_wb(0:im),bc4_wb(0:im)
c
      real cs
      real rnu_e(0:im,0:km)
c
      real rctime,pi,timer,ratio,pretime
c
      real rho
c
c     boundary condition flags
c     1:left, 2:right, 3:top, 4:bottom
      integer ib1u,ib1w,ib2u,ib2w ! for left and right of u,w
      integer ib3u,ib3w,ib4u,ib4w ! for top and bottom of u,w
      integer ib1p,ib2p,ib3p,ib4p ! for pressure
      integer ib1r,ib2r,ib3r,ib4r ! for density
c
c     masking flag: ifl=1 for obstacle cell
      integer ifl(0:im,0:km),ifu(0:im,0:km),ifw(0:im,0:km)
      integer ifr(0:im,0:km)
c
c     source flag
      integer iscf
      integer isu_u(0:im,0:km),isw_w(0:im,0:km),isr_c(0:im,0:km) ! source flag
c
c     source value      
      real scu_u(0:im,0:km),scw_w(0:im,0:km),scr_c(0:im,0:km) ! source value     
      real scu_ub(0:im,0:km),scw_wb(0:im,0:km) ! source value transformed
c
      real chleng,height
c
      rctime=1.0
      pi=3.14159265
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
c     /* reading boundary condition flag */
      open(12,file='boundary.d',status='old')
      read(12,*) ib1u,ib1w,ib2u,ib2w
      read(12,*) ib3u,ib3w,ib4u,ib4w
      read(12,*) ib1p,ib2p,ib3p,ib4p
      read(12,*) ib1r,ib2r,ib3r,ib4r
      close(12)
c
c *** reading location data ***
c =====================================
      call geodata2dv(x,z,xi,zet,x_c,z_c,dxi,dzet,dxa,dza
     2                  ,ifl,ifu,ifw,ifr)
c =====================================
c
c  
c *** setting coordinate transformation coefficient ***
c =====================================
      call gcoeff(dxi,dzet,x,z,x_xi,z_zet,xi_x,zet_z,rj)
c =====================================
c
c
c *** initialization ***
c =====================================
      call initialization(u_c,w_c,ub,wb,ubn,wbn
     1                  ,yp,ypn,yr,yrn
     2                  ,gux,guz,gwx,gwz,grx,grz
     3                  ,ub_l,wb_l,u_l,w_l,yp_l,yr_l)     
c =====================================
c
c
c *** setting initial condition ***
c =====================================
      call initial1(u_c,w_c,u_l,w_l,yp_l,yp,ypn,yr_l,yr,yrn
     1           ,bc1_u,bc2_u,bc3_w,bc4_w,bc1_r,bc2_r,bc3_r,bc4_r
     2           ,iscf,scu_u,scw_w,scr_c,isu_u,isw_w,isr_c)
      call initial2(u_c,w_c,ub,wb,ubn,wbn,xi_x,zet_z
     &        ,bc1_u,bc2_u,bc3_w,bc4_w,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     2        ,iscf,scu_u,scw_w,scu_ub,scw_wb)
      call bcset_v(ubn,wbn,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
      call bcset_v(ub,wb,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
      call bcset_r(yr,ifr
     1     ,ib1r,ib2r,ib3r,ib4r
     2     ,bc1_r,bc2_r,bc3_r,bc4_r,rho
     3     ,iscf,scr_c,isr_c)
      call bcset_r(yrn,ifr
     1     ,ib1r,ib2r,ib3r,ib4r
     2     ,bc1_r,bc2_r,bc3_r,bc4_r,rho)
      call each_v(ub,wb,ub_w,wb_u
     3     ,iscf,scr_c,isr_c)
      call initial3(dxi,dzet,gux,guz,gwx,gwz,grx,grz,ubn,wbn,yrn)
c =====================================
c
      do i=0,nx+1
        do k=0,nz+1
         rnu_e(i,k)=rnu
        end do
      end do
c
c     /* setting output counter */
      itout=nint(tuk/dt)
      time=0.0
      itt=itout
      icount=0
c
c *** start of time step iteration ***
 2000 continue
      if(itt.eq.itout) then
       itt=0
c
c *** set velocity at each points ***
       call each_v(ub,wb,ub_w,wb_u)
c
c *** data output ***
        icount=icount+1
c =====================================
        call output(yp_l,yr_l,u_l,w_l,x,z,ifl,icount)
c =====================================
       write(*,*) 'step=',icount,'   time=',time
c
      end if 
c
c     set mitigation for next step results by "pretime"
      if (time.lt.pretime) then
       timer=time/pretime*pi
       ratio=(1.0-0.0)*0.5*sin(1.5*pi+timer)
     1       +0.5*(1.0-0.0)+0.0
       rctime=1.0-ratio
      else
       ratio=1.0
       rctime=0.0
      end if
c
c *** changing temporary step counter ***
      itt=itt+1
c
 4000 continue
c
      do i=0,nx+1	 
        do k=0,nz+1
         ubt(i,k)=ub(i,k)
         wbt(i,k)=wb(i,k)
        end do
      end do
c
      call sor_cal(lsor,dxi,dzet,rho,rj,ubn,wbn
     1          ,ypn,ypn3,soralp,xi_x,zet_z,ifl,ifu,ifw
     2          ,ib1p,ib2p,ib3p,ib4p,ib1u,ib2u)
c
      call press_cal(rho,dxi,dzet,ubn,wbn,ypn,ypn3
     1           ,ubt,wbt,xi_x,zet_z)
c
      call bcset_v(ubn,wbn,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
c
      call each_v(ubn,wbn,ub_w,wb_u) 
      call buoyancy_cal(wbn,yrn,zet_z,g,rho)
      call bcset_v(ubn,wbn,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
c
 5000 continue
c
c
c
c =====================================
c
c
c *** calculation for viscosity and diffusion ***
c =====================================
      do i=0,nx+1
        do k=0,nz+1
         ubt(i,k)=ubn(i,k)
         wbt(i,k)=wbn(i,k)
         yrt(i,k)=yrn(i,k)
        end do
      end do
c
      call each_v(ubn,wbn,ub_w,wb_u)
      call visco_cal(rnu_e,rnu,diffyr
     1              ,dxi,dzet,ubt,wbt,yrt,ubn,wbn,yrn
     2              ,xi_x,zet_z) 
      call bcset_v(ubn,wbn,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
      call bcset_r(yrn,ifr
     1     ,ib1r,ib2r,ib3r,ib4r
     2     ,bc1_r,bc2_r,bc3_r,bc4_r,rho
     3     ,iscf,scr_c,isr_c)
c
c =====================================
c
c *** calculation for advection term with CIP method ***
c =====================================
      do i=0,nx+1
        do k=0,nz+1
         ubt(i,k)=ubn(i,k)
         wbt(i,k)=wbn(i,k)
         yrt(i,k)=yrn(i,k)
        end do
      end do
c
      call each_v(ubn,wbn,ub_w,wb_u) 
c
      call newgrd(dxi,dzet,nx-1,nz,ubn,ubt,gux,guz,1,ifu)
      call bcset_g(1,gux,guz,ifu,ib1u,ib2u)
      call cip_cal(1,dxi,dzet,ubn,ubt,wb_u,gux,guz)
      call bcset_v(ubn,wbn,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
      call gr_cal(1,dxi,dzet,ubn,ubt,wb_u,gux,guz)
c
      call newgrd(dxi,dzet,nx,nz-1,wbn,wbt,gwx,gwz,3,ifw)
      call bcset_g(3,gwx,gwz,ifw,ib1u,ib2u)
      call cip_cal(3,dxi,dzet,wbn,ub_w,wbt,gwx,gwz)
      call bcset_v(ubn,wbn,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
      call gr_cal(3,dxi,dzet,wbn,ub_w,wbt,gwx,gwz)
c
      call newgrd(dxi,dzet,nx,nz,yrn,yrt,grx,grz,4,ifr)
      call bcset_g(4,grx,grz,ifr,ib1u,ib2u)
      call cip_cal(4,dxi,dzet,yrn,ubt,wbt,grx,grz)
      call bcset_r(yrn,ifr
     1     ,ib1r,ib2r,ib3r,ib4r
     2     ,bc1_r,bc2_r,bc3_r,bc4_r,rho)
      call gr_cal(4,dxi,dzet,yrn,ubt,wbt,grx,grz
     3     ,iscf,scr_c,isr_c)
c
      call bcset_v(ubn,wbn,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
      call bcset_r(yrn,ifr
     1     ,ib1r,ib2r,ib3r,ib4r
     2     ,bc1_r,bc2_r,bc3_r,bc4_r,rho
     3     ,iscf,scr_c,isr_c)
c
      call bcset_g(1,gux,guz,ifu)
      call bcset_g(3,gwx,gwz,ifw)
      call bcset_g(4,grx,grz,ifr)
c
c =====================================
c
c *** shifting temporaly valiables to this step result ***
c =====================================
c
      call timerelax(ubn,ub,wbn,wb,ypn,yp,yrn,yr,rctime)
c
c *** shifting temporaly valiables to this step result ***
c =====================================
c    
      call bcset_v(ubn,wbn,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
      call bcset_r(yrn,ifr
     1     ,ib1r,ib2r,ib3r,ib4r
     2     ,bc1_r,bc2_r,bc3_r,bc4_r,rho
     3     ,iscf,scr_c,isr_c)
      call step_shift(ubn,ub,wbn,wb,ypn,yp,yrn,yr,x_xi,z_zet
     3                ,ub_l,wb_l,yp_l,yr_l,u_l,w_l,ifl)
c
      call bcset_v(ub,wb,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
      call bcset_r(yr,ifr
     1     ,ib1r,ib2r,ib3r,ib4r
     2     ,bc1_r,bc2_r,bc3_r,bc4_r,rho
     3     ,iscf,scr_c,isr_c)
c =====================================
c
c *** end of time step iteration ***
      if(time.gt.etime) then
       goto 3000
      end if
      time=time+dt
      goto 2000
c
 3000	continue
c
c
      write(6,*) 'calculation is normally finished'
c
      end
c
c *** End of the Main Program ***
c
c
c
c
c
c
c##############################################################################
c#                                                                            #
c#                 following source codes are subroutines                     #
c#                                                                            #
c##############################################################################
c
c1
c2
c3
c4
c5
c     ***************************************************************
c     *                                                             *   
c     *                Grid Location Data Reading                   *
c     *															   *
c     ***************************************************************
c
c ========== subroutine for grid location data reading ============
      subroutine geodata2dv(x,z,xi,zet,x_c,z_c,dxi,dzet,dxa,dza
     2                     ,ifl,ifu,ifw,ifr)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      integer id,kd ! dummy
      real dxi,dzet
      real x(0:im,0:km),z(0:im,0:km)
     1    ,x_c(0:im,0:km),z_c(0:im,0:km)
      real xi(0:im),zet(0:km)
      real dxa(0:im,0:km),dza(0:im,0:km)
c
      integer ifl(0:im,0:km),ifu(0:im,0:km),ifw(0:im,0:km)
      integer ifr(0:im,0:km)
      integer icc
c
      do i=0,nx+1
        do k=0,nz+1
         ifl(i,k)=0
         ifu(i,k)=0
         ifw(i,k)=0
         ifr(i,k)=0
        end do
      end do
c
      open(16,file='geometry.dat',status='old')
c
      read(16,*) nx,nz
      do i=0,nx
        do k=0,nz
c         read(16,'(2f12.6,I8)')
c     & x(i,k),z(i,k),ifl(i,k)
         read(16,*) x(i,k),z(i,k),ifl(i,k)
        end do
      end do
c
      close(16)
c
c     /* calculation for location at center of each cell */
      do i=1,nx
        do k=1,nz
         x_c(i,k)=(x(i,k)+x(i-1,k)+x(i,k-1)+x(i-1,k-1))*0.25
         z_c(i,k)=(z(i,k)+z(i-1,k)+z(i,k-1)+z(i-1,k-1))*0.25
        end do
      end do
c
c     /* calculation for actual beam size of each cell */
      do i=1,nx
        do k=1,nz
         dxa(i,k)=sqrt((x(i,k)-x(i-1,k))**2)
         dza(i,k)=sqrt((z(i,k)-z(i,k-1))**2)
        end do
      end do
c
c     /* setting position of xi and zeta */
      dxi=1.0/float(nx)
      dzet=1.0/float(nz)
      do i=0,nx
       xi(i)=(float(i)-0.5)*dxi
      end do
      do k=0,nz
       zet(k)=(float(k)-0.5)*dzet
      end do
c
      do i=0,nx+1
         ifl(i,nz+1)=ifl(i,nz)
         ifl(i,0)=ifl(i,1)
      end do
      do k=0,nz+1
         ifl(0,k)=ifl(1,k)
         ifl(nx+1,k)=ifl(nx,k)
      end do
c
      do i=0,nx
        do k=1,nz
         if(ifl(i,k).eq.1 .and. ifl(i,k-1).eq.1) then
          ifu(i,k)=1
         end if
        end do
      end do
c
      do i=1,nx
        do k=0,nz
         if(ifl(i,k).eq.1.and.ifl(i-1,k).eq.1) then
          ifw(i,k)=1
         end if
        end do
      end do
c
      do i=1,nx
        do k=1,nz
         icc=ifl(i,k)+ifl(i-1,k)+ifl(i,k-1)+ifl(i-1,k-1)
         if(icc.eq.4) then
             ifr(i,k)=1
         end if
        end do
      end do
c
      return
      end
c
cc
ccc
cccc
ccccc
cccccc
c     ***************************************************************
c     *                                                             *   
c     *       Setting Cordinate Transformation Coefficient          *
c     *                                                             *
c     ***************************************************************
c
c ====== subroutine for transformation coefficient setting ========
      subroutine gcoeff(dxi,dzet,x,z,x_xi,z_zet,xi_x,zet_z,rj)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real x(0:im,0:km),z(0:im,0:km)
      real x_xi(0:im,0:km),z_zet(0:im,0:km)
      real xi_x(0:im,0:km),zet_z(0:im,0:km)
      real rj(0:im,0:km)
      real dxi,dzet
c
c     masking flag: ifl=1 for obstacle cell 
      integer ifl(0:im,0:km)
c
      do i=0,nx 
        x(i,nz+1)=x(i,nz)
        z(i,nz+1)=z(i,nz)+(z(i,nz)-z(i,nz-1))
      end do 
      do k=0,nz 
        x(nx+1,k)=x(nx,k)+(x(nx,k)-x(nx-1,k))
        z(nx+1,k)=z(nx,k)
      end do               
c
      do i=1,nx
        do k=1,nz
         x_xi(i,k)=(x(i,k)-x(i-1,k)+x(i,k-1)-x(i-1,k-1)
     2       )*0.5/dxi
c
         z_zet(i,k)=(z(i,k)-z(i,k-1)+z(i-1,k)-z(i-1,k-1)
     2       )*0.5/dzet
        end do
      end do
c
      do i=1,nx
        x_xi(i,nz+1)=x_xi(i,nz)
        z_zet(i,nz+1)=z_zet(i,nz)
        x_xi(i,0)=x_xi(i,1)
        z_zet(i,0)=z_zet(i,1)
      end do        
c
      do i=1,nx
       do k=1,nz+1
        if(ifl(i,k).eq.1.and.ifl(i-1,k).eq.0) then
         z_zet(i,k)=z_zet(i-1,k)
        end if
        if(ifl(i,k).eq.1.and.ifl(i+1,k).eq.0) then
         z_zet(i,k)=z_zet(i+1,k)
        end if
       end do
      end do        
c
      do k=1,nz+1
        x_xi(0,k)=x_xi(1,k)
        z_zet(0,k)=z_zet(1,k)
        x_xi(nx+1,k)=x_xi(nx,k)
        z_zet(nx+1,k)=z_zet(nx,k)
      end do
c
      do i=1,nx
        do k=1,nz+1 
         rj(i,k)=1.0/(x_xi(i,k)*z_zet(i,k))
        end do
      end do
c
      do k=1,nz+1
        rj(0,k)=rj(1,k)
        rj(nx+1,k)=rj(nx,k)
      end do
c
      do i=0,nx+1
        do k=1,nz+1
         xi_x(i,k)=rj(i,k)*(z_zet(i,k))
         zet_z(i,k)=rj(i,k)*(x_xi(i,k))
        end do
      end do
c
      do i=1,nx   
       do k=1,nz+1
        if(ifl(i,k).eq.1.and.ifl(i-1,k).eq.0) then
         zet_z(i,k)=zet_z(i-1,k)
        end if 
        if(ifl(i,k).eq.1.and.ifl(i+1,k).eq.0) then
         zet_z(i,k)=zet_z(i+1,k)
        end if 
       end do
      end do        
c
      do i=1,nx
        xi_x(i,nz+1)=xi_x(i,nz)
        zet_z(i,nz+1)=zet_z(i,nz)
        xi_x(i,0)=xi_x(i,1)
        zet_z(i,0)=zet_z(i,1)
      end do
c
      do k=1,nz+1
        xi_x(0,k)=xi_x(1,k)
        xi_x(nx+1,k)=xi_x(nx,k)
      end do 
c
      return
      end
c
cc
ccc
cccc
ccccc
cccccc
c     ***************************************************************
c     *                                                             *
c     *                Setting Initial Condition                    *
c     *                                                             *
c     ***************************************************************
c
c ========== subroutine for initialization ========================
      subroutine initialization(u_c,w_c,ub,wb,ubn,wbn
     1                  ,yp,ypn,yr,yrn
     2                  ,gux,guz,gwx,gwz,grx,grz
     3                  ,ub_l,wb_l,u_l,w_l,yp_l,yr_l)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real u_c(0:im,0:km),w_c(0:im,0:km)
      real ub_l(0:im,0:km),wb_l(0:im,0:km)
      real ub(0:im,0:km),wb(0:im,0:km)
      real u_l(0:im,0:km),w_l(0:im,0:km)
     2    ,ubn(0:im,0:km),wbn(0:im,0:km)
      real yp(0:im,0:km),ypn(0:im,0:km)
      real yp_l(0:im,0:km)
      real yr(0:im,0:km),yrn(0:im,0:km)
      real yr_l(0:im,0:km)
      real gux(0:im,0:km),guz(0:im,0:km)
      real gwx(0:im,0:km),gwz(0:im,0:km)
      real grx(0:im,0:km),grz(0:im,0:km)
c
      do i=0,nx+1
        do k=0,nz+2
         u_c(i,k)=0.0
         w_c(i,k)=0.0
        end do
      end do
c
      do i=0,nx
        do k=0,nz+1
         u_l(i,k)=0.0
         w_l(i,k)=0.0
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz+2
         ub(i,k)=0.0
         wb(i,k)=0.0
         ubn(i,k)=0.0
         wbn(i,k)=0.0
         ub_l(i,k)=0.0
         wb_l(i,k)=0.0
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz+2
         yp(i,k)=0.0
         ypn(i,k)=0.0
         yp_l(i,k)=0.0
         yr(i,k)=0.0
         yrn(i,k)=0.0
         yr_l(i,k)=0.0
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz+1
         gux(i,k)=0.0
         guz(i,k)=0.0
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz+1
         gwx(i,k)=0.0
         gwz(i,k)=0.0
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz+1
         grx(i,k)=0.0
         grz(i,k)=0.0
        end do
      end do
c
      return
      end
c
c
c
c
c ========== subroutine for grid location data reading ============
      subroutine initial1(u_c,w_c,u_l,w_l,yp_l,yp,ypn,yr_l,yr,yrn
     1           ,bc1_u,bc2_u,bc3_w,bc4_w,bc1_r,bc2_r,bc3_r,bc4_r
     2           ,iscf,scu_u,scw_w,scr_c,isu_u,isw_w,isr_c)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      integer idm,kdm ! dummy
      real u_c(0:im,0:km),w_c(0:im,0:km)
      real u_l(0:im,0:km),w_l(0:im,0:km)
      real yp_l(0:im,0:km),yp(0:im,0:km),ypn(0:im,0:km)
      real yr_l(0:im,0:km),yr(0:im,0:km),yrn(0:im,0:km)
      real bc1_u(0:km),bc2_u(0:km) ! b.c.flow at U part on left & right
      real bc3_w(0:im),bc4_w(0:im) ! b.c.flow at W part on top & bottom
      real bc1_r(0:km),bc2_r(0:km) ! b.c.rho at center on left & right
      real bc3_r(0:im),bc4_r(0:im) ! b.c.rho at center on top & bottom
      integer iscf ! source flag
      integer isu_u(0:im,0:km),isw_w(0:im,0:km),isr_c(0:im,0:km) ! source flag
      real scu_u(0:im,0:km),scw_w(0:im,0:km),scr_c(0:im,0:km) ! source value
c
      if(iscf.eq.1) then
       open(17,file='source.dat',status='old',form='formatted')
       read(17,*) idm,kdm
       do i=0,nx
        do k=0,nz
         read(17,'(3f12.6,3I8)')
     &  scu_u(i,k),scw_w(i,k),scr_c(i,k)
     & ,isu_u(i,k),isw_w(i,k),isr_c(i,k)
        end do
       end do
       close(17)
      end if
c
      open(22,file='initial.dat',status='old',form='formatted')
      read(22,*) idm,kdm
      do i=0,nx
        do k=0,nz
         read(22,*) u_l(i,k),w_l(i,k),yp_l(i,k),yr_l(i,k) ! at node
        end do
      end do
      close(22)
c
      open(23,file='bc_vert.dat',status='old',form='formatted')
      read(23,*) kdm
      do k=1,nz
         read(23,*) bc1_u(k),bc2_u(k),bc1_r(k),bc2_r(k) ! at u part
      end do
      close(23)
c 
      open(24,file='bc_horz.dat',status='old',form='formatted')
      read(24,*) idm
      do i=1,nx
         read(24,*) bc3_w(i),bc4_w(i),bc3_r(i),bc4_r(i) ! at w part
      end do
      close(24)
c
      do i=1,nx
        do k=1,nz
         u_c(i,k)=(u_l(i,k)+u_l(i-1,k)+u_l(i,k-1)+u_l(i-1,k-1))
     2     *0.25
         w_c(i,k)=(w_l(i,k)+w_l(i-1,k)+w_l(i,k-1)+w_l(i-1,k-1))
     2     *0.25
         yp(i,k)=(yp_l(i,k)+yp_l(i-1,k)+yp_l(i,k-1)+yp_l(i-1,k-1))
     2     *0.25
         yr(i,k)=(yr_l(i,k)+yr_l(i-1,k)+yr_l(i,k-1)+yr_l(i-1,k-1))
     2     *0.25
        end do
      end do
c
      do i=1,nx
         yp(i,0)=(yp_l(i,0)+yp_l(i-1,0))/2.
         yr(i,0)=(yr_l(i,0)+yr_l(i-1,0))/2.
         yp(i,nz+1)=(yp_l(i,nz)+yp_l(i-1,nz))/2.
         yr(i,nz+1)=(yr_l(i,nz)+yr_l(i-1,nz))/2.
      end do
      do k=0,nz+1
         yp(0,k)=yp(1,k)
         yr(0,k)=yr(1,k)
         yp(nx+1,k)=yp(nx,k)
         yr(nx+1,k)=yr(nx,k)
      end do
c
      do i=0,nx+1
        do k=0,nz+1
         ypn(i,k)=yp(i,k)
         yrn(i,k)=yr(i,k)
        end do
      end do
c
      return
	end
c
cc
ccc
cccc
ccccc
cccccc
c ========== subroutine for grid location data reading ============
      subroutine initial2(u_c,w_c,ub,wb,ubn,wbn,xi_x,zet_z
     &        ,bc1_u,bc2_u,bc3_w,bc4_w,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     2        ,iscf,scu_u,scw_w,scu_ub,scw_wb)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real u_c(0:im,0:km),w_c(0:im,0:km)
      real ub_c(0:im,0:km),wb_c(0:im,0:km) ! coordinate transformed
      real ub(0:im,0:km),wb(0:im,0:km) ! coordinate transformed
      real ubn(0:im,0:km),wbn(0:im,0:km) ! coordinate transformed
c
      real xi_x(0:im,0:km),zet_z(0:im,0:km)
c
      real bc1_u(0:km),bc2_u(0:km)
      real bc3_w(0:im),bc4_w(0:im)
c
      real bc1_ub(0:km),bc2_ub(0:km)! coordinate transformed
      real bc3_wb(0:im),bc4_wb(0:im)! coordinate transformed
c
      integer iscf ! source flug 0:off,1:on
      real scu_u(0:im,0:km),scw_w(0:im,0:km) ! source value
      real scu_ub(0:im,0:km),scw_wb(0:im,0:km) ! source value coordinate transformed
c
c
      do i=1,nx
        do k=1,nz
         ub_c(i,k)=xi_x(i,k)*u_c(i,k)
         wb_c(i,k)=zet_z(i,k)*w_c(i,k)
        end do
      end do
c
      do k=1,nz
         bc1_ub(k)=xi_x(1,k)*bc1_u(k)
         bc2_ub(k)=xi_x(nx,k)*bc2_u(k)
      end do
      do i=1,nx
         bc3_wb(i)=zet_z(i,nz)*bc3_w(i)
         bc4_wb(i)=zet_z(i,1)*bc4_w(i)
      end do
c
      do i=1,nx-1
        do k=1,nz
         ub(i,k)=(ub_c(i+1,k)+ub_c(i,k))*0.5
        end do
      end do
c
      do i=1,nx
        do k=1,nz-1
         wb(i,k)=(wb_c(i,k+1)+wb_c(i,k))*0.5
        end do
      end do
c
      if(iscf.eq.1) then
       do i=1,nx-1
         do k=1,nz
          scu_ub(i,k)=xi_x(i,k)*scu_u(i,k)
         end do
       end do
       do i=1,nx
         do k=1,nz-1
          scw_wb(i,k)=zet_z(i,k)*scw_w(i,k)
         end do
       end do
      end if
c
      do i=1,nx-1
        do k=1,nz
         ubn(i,k)=ub(i,k)
        end do
      end do
c
      do i=1,nx
        do k=1,nz-1
         wbn(i,k)=wb(i,k)
        end do
      end do
c
      return
      end
c
cc
ccc
cccc
ccccc
cccccc
c ========== subroutine for grid location data reading ============
      subroutine initial3(dxi,dzet,gux,guz,gwx,gwz,grx,grz,ubn,wbn,yrn)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real dxi,dzet
      real ubn(0:im,0:km),wbn(0:im,0:km),yrn(0:im,0:km)
      real gux(0:im,0:km),guz(0:im,0:km)
      real gwx(0:im,0:km),gwz(0:im,0:km)
      real grx(0:im,0:km),grz(0:im,0:km)
c
      do i=1,nx-1
        do k=1,nz
         gux(i,k)=(ubn(i+1,k)-ubn(i-1,k))/2.0/dxi
         guz(i,k)=(ubn(i,k+1)-ubn(i,k-1))/2.0/dzet
        end do
      end do
c
      do i=1,nx
        do k=1,nz-1
         gwx(i,k)=(wbn(i+1,k)-wbn(i-1,k))/2.0/dxi
         gwz(i,k)=(wbn(i,k+1)-wbn(i,k-1))/2.0/dzet
        end do
      end do
c
      do i=1,nx
        do k=1,nz
         grx(i,k)=(yrn(i+1,k)-yrn(i-1,k))/2.0/dxi
         grz(i,k)=(yrn(i,k+1)-yrn(i,k-1))/2.0/dzet
        end do
      end do
c
      return
      end
c
c     ***************************************************************
c     *                                                             *
c     *               Data Output for Binary File                   *
c     *															  *
c     ***************************************************************
c
c ========== subroutine for data output ===========================
      subroutine output(yp_l,yr_l,u_l,w_l,x,z,ifl,icount)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real u_l(0:im,0:km),w_l(0:im,0:km)
      real yp_l(0:im,0:km)
      real yr_l(0:im,0:km)
c
      real x(0:im,0:km),z(0:im,0:km)
c
      integer ifl(0:im,0:km)
c
      integer icount
      integer namecount2
      character*13 buff11
c
      namecount2=icount-1
c
      if(namecount2.lt.10) then
          write(buff11,205) namecount2
      else if((namecount2.gt.9).and.(namecount2.lt.100)) then
          write(buff11,206) namecount2
      else if((namecount2.gt.99).and.(namecount2.lt.1000)) then
          write(buff11,207) namecount2
      else if((namecount2.gt.999).and.(namecount2.lt.10000)) then
          write(buff11,208) namecount2
      else
          write(buff11,209) namecount2
      end if
      write(6,*) buff11
c
      open(26,file=buff11,status='unknown',form='unformatted')
c
      write(26) time
      write(26) nx+1,nz+1
      write(26) ((x(i,k),i=0,nx),k=0,nz)
      write(26) ((z(i,k),i=0,nx),k=0,nz)
      write(26) ((u_l(i,k),i=0,nx),k=0,nz)
      write(26) ((w_l(i,k),i=0,nx),k=0,nz)
      write(26) ((ifl(i,k),i=0,nx),k=0,nz)
      write(26) ((yp_l(i,k),i=0,nx),k=0,nz)
      write(26) ((yr_l(i,k),i=0,nx),k=0,nz)
c
      close(26)
c
 205  format('rslt0000',i1,'.dat') 
 206  format('rslt000',i2,'.dat')
 207  format('rslt00',i3,'.dat') 
 208  format('rslt0',i4,'.dat')  
 209  format('rslt',i5,'.dat')
c
      return
      end
c
c1
c2
c3
c4
c5
c     ***************************************************************
c     *                                                             *
c     *           Pressure Term Calculation by SOR                  *
c     *                                                             *
c     ***************************************************************
c
c ========== subroutine for pressure calculation by SOR ===========
      subroutine sor_cal(lsor,dxi,dzet,rho,rj,ubn,wbn
     1          ,ypn,ypn3,soralp,xi_x,zet_z,ifl,ifu,ifw
     2          ,ib1p,ib2p,ib3p,ib4p,ib1u,ib2u)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k,l,lsor
      real soralp
      real dxi,dzet,rho
      real a_e(0:im,0:km),a_w(0:im,0:km)
     2    ,a_u(0:im,0:km),a_d(0:im,0:km)
      real udr(0:im,0:km),wdr(0:im,0:km)
      real a_p(0:im,0:km),rj(0:im,0:km)
     1     ,a_f(0:im,0:km)
	real xi_x(0:im,0:km),zet_z(0:im,0:km)
      real alpha(0:im,0:km),gumma(0:im,0:km)
      real div
c     /* temporary velocity on xi, eta and zeta */
      real ubn(0:im,0:km),wbn(0:im,0:km)
c     /* Jacobian at u, v and w position */
      real rju(0:im,0:km),rjw(0:im,0:km)
c     /* temporary valuables iterative calculation for pressure */
      real work(0:im,0:km)
      real ypn(0:im,0:km),ypn3(0:im,0:km)
      real ww 
c
      integer ifl(0:im,0:km),ifu(0:im,0:km)
     1        ,ifw(0:im,0:km)
c
      integer ib1p,ib2p,ib3p,ib4p
      integer ib1u,ib2u
c
      do i=1,nx
        do k=1,nz
         alpha(i,k)
     1      =xi_x(i,k)**2
         gumma(i,k)
     1      =zet_z(i,k)**2
        end do
      end do
c
      do i=1,nx
        do k=1,nz
         if(abs(rj(i,k)).gt.0.0000001) then
          a_e(i,k)=1.0/rj(i,k)/rho/dxi/dxi*alpha(i,k)
          a_w(i,k)=1.0/rj(i,k)/rho/dxi/dxi*alpha(i,k)
          a_u(i,k)=1.0/rj(i,k)/rho/dzet/dzet*gumma(i,k)
          a_d(i,k)=1.0/rj(i,k)/rho/dzet/dzet*gumma(i,k)
         end if
        end do
      end do
c
c     /* boundary conditon */
      do k=1,nz
        a_w(1,k)=0.0
      end do 
      do k=1,nz
        a_e(nx,k)=0.0
      end do
      do i=1,nx
        a_u(i,nz)=0.0
      end do
      do i=1,nx
        a_d(i,1)=0.0
      end do
c
      do i=1,nx
        do k=0,nz+1
         if(ifl(i,k).eq.1) then
          a_u(i,k)=0.0
          a_d(i,k)=0.0
          a_w(i,k)=0.0
          a_e(i,k)=0.0
         end if
        end do
      end do
c
      do i=1,nx 
       do k=1,nz+1
        if(ifl(i,k).eq.0.and.ifl(i+1,k).eq.1) then
         a_e(i,k)=0.0
        end if
        if(ifl(i,k).eq.0.and.ifl(i-1,k).eq.1) then
         a_w(i,k)=0.0
        end if
        if(ifl(i,k).eq.0.and.ifl(i,k-1).eq.1) then
         a_d(i,k)=0.0
        end if
       end do
      end do
c
c     /* calculation for Jacobian at velocity position */
      do i=1,nx-1
        do k=1,nz
         rju(i,k)=(rj(i,k)+rj(i+1,k))/2.0
        end do
      end do
      do i=1,nx
       do k=1,nz
        if(ifu(i,k).eq.1.and.ifu(i-1,k).eq.0) then
         rju(i,k)=rju(i-1,k)
        end if
        if(ifu(i,k).eq.1.and.ifu(i+1,k).eq.0) then
         rju(i,k)=rju(i+1,k)
        end if
       end do
      end do
c
      do i=1,nx
        do k=1,nz 
         rjw(i,k)=(rj(i,k)+rj(i,k+1))/2.0
        end do
      end do
      do k=1,nz
        rjw(0,k)=rjw(1,k)
        rjw(nx+1,k)=rjw(nx,k)
      end do
      do i=1,nx
        do k=1,nz
         if(ifl(i,k).eq.1.and.ifl(i,k+1).eq.0) then
          rjw(i,k)=rjw(i,k+1)
         end if
        end do
      end do
c
c     /* calculation for divergence */
      do i=1,nx-1
        do k=1,nz
         udr(i,k)=ubn(i,k)/rju(i,k)
        end do
      end do
c
      do i=1,nx
        do k=1,nz-1
         wdr(i,k)=wbn(i,k)/rjw(i,k)
        end do
      end do
      do i=1,nx
        wdr(i,0)=0.0
        wdr(i,nz)=0.0
      end do
c
      do i=1,nx
        do k=1,nz
         div=(udr(i,k)-udr(i-1,k))/dxi
     2        +(wdr(i,k)-wdr(i,k-1))/dzet
         a_p(i,k)=a_e(i,k)+a_w(i,k)
     2               +a_u(i,k)+a_d(i,k)
c
         a_f(i,k)=-div/dt
        end do
      end do
c
      do i=1,nx
        do k=0,nz+1
         if(ifl(i,k).eq.1) then
          a_f(i,k)=0.0
         end if
        end do
      end do
c
c     /* boundary condition for a_p */   p'=0
      if (ib1p.eq.1) then
       do k=1,nz
         a_p(1,k)=a_p(1,k)
     1      +2.0/rj(1,k)/rho/dxi/dxi*alpha(1,k)
       end do
      end if
      if (ib2p.eq.1) then
       do k=1,nz
         a_p(nx,k)=a_p(nx,k)
     1      +2.0/rj(nx,k)/rho/dxi/dxi*alpha(nx,k)
       end do
      end if
      if (ib3p.eq.1) then
       do i=1,nx
         a_p(i,nz)=a_p(i,nz)
     1      +2.0/rj(i,nz)/rho/dzet/dzet*gumma(i,nz)
       end do
      end if
      if (ib4p.eq.1) then
       do i=1,nx
         a_p(i,1)=a_p(i,1)
     1      +2.0/rj(i,1)/rho/dzet/dzet*gumma(i,1)
       end do
      end if
c
c     /* iterative calculation for pressure */
      do l=1,lsor
       do i=1,nx
         do k=1,nz
          if(abs(a_p(i,k)).gt.0.000001) then
           ww=(a_e(i,k)*ypn(i+1,k)+a_w(i,k)*ypn(i-1,k)
     2       +a_u(i,k)*ypn(i,k+1)+a_d(i,k)*ypn(i,k-1)
     3       +a_f(i,k))/a_p(i,k)
           work(i,k)=(1.0-soralp)*ypn(i,k)+soralp*ww
          else
           work(i,k)=0.0
          end if
         end do
       end do
       do i=1,nx
         do k=1,nz
          ypn(i,k)=work(i,k)
         end do
       end do
      end do
c
c     /* boundary condition for a_p */   p'=0
      if (ib1p.eq.1) then
       do k=1,nz
        ypn(0,k)=-ypn(1,k) 
       end do
      else
       do k=1,nz
        ypn(0,k)=ypn(1,k) 
       end do
      end if
c
      if (ib2p.eq.1) then
       do k=1,nz
        ypn(nx+1,k)=-ypn(nx,k) 
       end do
      else
       do k=1,nz
        ypn(nx+1,k)=ypn(nx,k) 
       end do
      end if
c
      if (ib3p.eq.1) then
       do i=1,nx
        ypn(i,nz+1)=-ypn(i,nz) 
       end do
      else
       do i=1,nx
        ypn(i,nz+1)=ypn(i,nz) 
       end do
      end if
c
      if (ib4p.eq.1) then
       do i=1,nx
        ypn(i,0)=-ypn(i,1) 
       end do
      else
       do i=1,nx
        ypn(i,0)=ypn(i,1) 
       end do
      end if
c
      do i=1,nx
        do k=0,nz+1
         if(ifl(i,k).eq.1) then
          ypn(i,k)=0.0 
         end if
        end do
      end do
      do i=0,nx+1
        do k=0,nz+1
         ypn3(i,k)=ypn(i,k)
        end do
      end do
c
      do i=1,nx
        do k=0,nz+1
         if(ifl(i,k).eq.1.and.ifl(i-1,k).eq.0) then 
          ypn(i,k)=ypn(i-1,k)
         end if
         if(ifl(i,k).eq.1.and.ifl(i+1,k).eq.0) then 
          ypn(i,k)=ypn(i+1,k)
         end if
        end do
      end do
c
      if(ib1u.eq.3 .and. ib2u.eq.3) then
        do k=1,nz
          ypn(0,k)=ypn(nx-3,k)
          ypn(1,k)=ypn(nx-2,k)
          ypn(nx-1,k)=ypn(2,k)
          ypn(nx,k)=ypn(3,k)
          ypn(nx+1,k)=ypn(4,k)
        end do
      end if
c
      do i=1,nx
        do k=1,nz+1
         if(ifl(i,k).eq.1.and.ifl(i,k+1).eq.0) then 
          ypn3(i,k)=ypn(i,k+1)
         end if
        end do
      end do
c
      return
      end
c
c
c ========== subroutine for presssure term ========================
      subroutine press_cal(rho,dxi,dzet,ubn,wbn,ypn,ypn3
     1           ,ubt,wbt,xi_x,zet_z)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real ubn(0:im,0:km),wbn(0:im,0:km),ypn(0:im,0:km)
      real ypn3(0:im,0:km)
      real ubt(0:im,0:km),wbt(0:im,0:km)
      real xi_x(0:im,0:km),zet_z(0:im,0:km)
      real rho,dxi,dzet
c     /* temporary variable */
      real x1,z3,px1,pz3
c
      do i=1,nx-1
        do k=1,nz
         x1=((xi_x(i,k)+xi_x(i+1,k))*0.5)**2
         px1=(ypn(i+1,k)-ypn(i,k))/dxi
         ubn(i,k)=ubt(i,k)-dt/rho*(x1*px1) 
        end do
      end do
c
      do i=1,nx
        do k=1,nz
         z3=((zet_z(i,k)+zet_z(i,k+1))*0.5)**2
         pz3=(ypn3(i,k+1)-ypn3(i,k))/dzet
         wbn(i,k)=wbt(i,k)-dt/rho*(z3*pz3) 
        end do
      end do
c
      return
      end
c
c
c
c
c ========== subroutine for viscosity term calculation =============
      subroutine visco_cal(rnu_e,rnu,diffyr
     1                     ,dxi,dzet,ub,wb,yr,ubn,wbn,yrn
     2                     ,xi_x,zet_z) 
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k

      real ubn(0:im,0:km),wbn(0:im,0:km),yrn(0:im,0:km)
      real ub(0:im,0:km),wb(0:im,0:km),yr(0:im,0:km)
      real ubt(0:im,0:km),wbt(0:im,0:km),yrt(0:im,0:km)
      real xi_x(0:im,0:km),zet_z(0:im,0:km)
      real xix13(0:im,0:km),zetz13(0:im,0:km)
      real rnu_e(0:im,0:km)
      real dxi,dzet,w1,w3,re1g,re3g,re1s,re3s
      real rnu,diffyr
c
      do i=0,nx+1
        do k=0,nz+1
         ubt(i,k)=ub(i,k)
         wbt(i,k)=wb(i,k)
         yrt(i,k)=yr(i,k)
        end do
      end do
c
      do i=1,nx-1 
        do k=1,nz-1
         zetz13(i,k)=0.25*
     1 (zet_z(i,k)+zet_z(i+1,k)+zet_z(i,k+1)+zet_z(i+1,k+1))
         xix13(i,k)=0.25*
     1   (xi_x(i,k)+xi_x(i+1,k)+xi_x(i,k+1)+xi_x(i+1,k+1))
        end do
      end do
c
      do k=1,nz-1
        zetz13(0,k)=zetz13(1,k)
        zetz13(nx,k)=zetz13(nx-1,k)
        xix13(0,k)=xix13(1,k)
        xix13(nx,k)=xix13(nx-1,k)
      end do
      do i=0,nx
        zetz13(i,0)=zetz13(i,1)
        zetz13(i,nz)=zetz13(i,nz-1)
        xix13(i,0)=xix13(i,1)
        xix13(i,nz)=xix13(i,nz-1)
      end do
c
      do i=1,nx-1
        do k=1,nz
         re1g=rnu_e(i+1,k) 
         re1s=rnu_e(i,k)
         re3g=(rnu_e(i,k)+rnu_e(i+1,k)
     1        +rnu_e(i,k+1)+rnu_e(i+1,k+1))*0.25
         re3s=(rnu_e(i,k)+rnu_e(i+1,k)
     1        +rnu_e(i,k-1)+rnu_e(i+1,k-1))*0.25
c
         w1=(re1g*(xi_x(i+1,k)**2.0*(ubt(i+1,k)-ubt(i,k))/dxi 
     1            +xi_x(i+1,k)**2.0*(ubt(i+1,k)-ubt(i,k))/dxi)
     2      -re1s*(xi_x(i,k)**2.0*(ubt(i,k)-ubt(i-1,k))/dxi
     3            +xi_x(i,k)**2.0*(ubt(i,k)-ubt(i-1,k))/dxi))/dxi
c
         w3=(re3g*(zetz13(i,k)**2.0*(ubt(i,k+1)-ubt(i,k))/dzet
     1            +xix13(i,k)**2.0*(wbt(i+1,k)-wbt(i,k))/dxi)
     2      -re3s*(zetz13(i,k-1)**2.0*(ubt(i,k)-ubt(i,k-1))/dzet
     3     +xix13(i,k-1)**2.0*(wbt(i+1,k-1)-wbt(i,k-1))/dxi))/dzet
c
         ubn(i,k)=ubn(i,k)+dt*(w1+w3) 
        end do
      end do
c
      do i=1,nx
        do k=1,nz
         re1g=(rnu_e(i,k)+rnu_e(i+1,k) 
     1        +rnu_e(i,k+1)+rnu_e(i+1,k+1))*0.25
         re1s=(rnu_e(i,k)+rnu_e(i-1,k)
     1        +rnu_e(i,k+1)+rnu_e(i-1,k+1))*0.25
         re3g=rnu_e(i,k+1)
         re3s=rnu_e(i,k) 
c
         w1=(re1g*(xix13(i,k)**2.0*(wbt(i+1,k)-wbt(i,k))/dxi
     1            +zetz13(i,k)**2.0*(ubt(i,k+1)-ubt(i,k))/dzet)
     2      -re1s*(xix13(i-1,k)**2.0*(wbt(i,k)-wbt(i-1,k))/dxi
     3    +zetz13(i-1,k)**2.0*(ubt(i-1,k+1)-ubt(i-1,k))/dzet))/dxi
c
         w3=(re3g*(zet_z(i,k+1)**2.0*(wbt(i,k+1)-wbt(i,k))/dzet
     1            +zet_z(i,k+1)**2.0*(wbt(i,k+1)-wbt(i,k))/dzet)
     2      -re3s*(zet_z(i,k)**2.0*(wbt(i,k)-wbt(i,k-1))/dzet
     3          +zet_z(i,k)**2.0*(wbt(i,k)-wbt(i,k-1))/dzet))/dzet
c
         wbn(i,k)=wbn(i,k)+dt*(w1+w3) 
        end do
      end do
c
      do i=1,nx
        do k=1,nz
c
         w1=( (yrt(i+1,k)-yrt(i,k))/dxi
     2       -(yrt(i,k)-yrt(i-1,k))/dxi )/dxi
     3             *rnu_e(i,k)*diffyr*xi_x(i,k)**2.0
c
         w3=((yrt(i,k+1)-yrt(i,k))/dzet
     2       -(yrt(i,k)-yrt(i,k-1))/dzet )/dzet
     3             *rnu_e(i,k)*diffyr*zet_z(i,k)**2.0
c
         yrn(i,k)=yrn(i,k)+dt*(w1+w3) 
        end do
      end do
c
      return
      end
c
c1
c2
c3
c4
c5
c     ***************************************************************
c     *                                                             *
c     *          Calculation for Advection Term                     *
c     *                                                             *
c     ***************************************************************
c
c ========== subroutine for buoyancy calculation =======
      subroutine buoyancy_cal(wbn,yrn,zet_z,g,rho)
c =================================================================
c
      implicit none
      include 'common.h'
c
c
      integer i,k
      real wbn(0:im,0:km)
      real yrn(0:im,0:km)
      real zet_z(0:im,0:km)
c      
      real g,rho
c
c     /* differensial calculation for wbn */
      do i=1,nx
        do k=1,nz-1 
          wbn(i,k)=wbn(i,k)-dt*g
     1          *(zet_z(i,k)+zet_z(i,k+1))*0.5
     2          *((yrn(i,k)+yrn(i,k+1))*0.5-rho)/rho
        end do
      end do
c
      return
      end
c
c
c ========== subroutine for advection term calculation by CIP =====
      subroutine cip_cal(itp,dxi,dzet,fn,ux,uz,gx,gz)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k,itp,nxt,nzt
      integer im1,km1 
      integer isn,ksn
      real dxi,dzet,dx1,dz1,dx2,dz2,dx3,dz3
     1    ,dtdx,dtdz
      real xx,zz,cx,cz
     1    ,a1,e1,b1,f1,tmp,tmq,d1,c1,g1
     1    ,gxt,gzt
c     /* temporary valuables for cip routine */
      real fn(0:im,0:km),f(0:im,0:km)
      real u(0:im,0:km),w(0:im,0:km)
      real ux(0:im,0:km),uz(0:im,0:km)
      real gx(0:im,0:km),gz(0:im,0:km)
      real gxn(0:im,0:km),gzn(0:im,0:km)
c
c     /* velocity at u-point for ubn calculation */
      if (itp.eq.1) then
       nxt=nx-1
       nzt=nz
       do k=0,nz+1
         do i=0,nx
          u(i,k)=ux(i,k)
          w(i,k)=uz(i,k)
         end do
       end do
c
c     /* velocity at w-point for yrn calculation */
      else if (itp.eq.3) then
       nxt=nx
       nzt=nz-1
       do k=0,nz
         do i=0,nx+1
          u(i,k)=ux(i,k)
          w(i,k)=uz(i,k)
         end do
       end do
c
c     /* velocity at center-point for yrn calculation */
      else
       nxt=nx
       nzt=nz
       do k=1,nz
         do i=1,nx
          u(i,k)=0.5*(ux(i,k)+ux(i-1,k))
          w(i,k)=0.5*(uz(i,k)+uz(i,k-1))
         end do
       end do
       
      end if
c
c     /* setting f */
      do i=0,nx+1
        do k=0,nz+1
         f(i,k)=fn(i,k)
        end do
      end do
c
c /* cip */
      dx1=dxi
      dx2 =dx1*dx1
      dx3 =dx2*dx1
      dtdx=dt/dx1
      dz1=dzet
      dz2 =dz1*dz1
      dz3 =dz2*dz1
      dtdz=dt/dz1
c
      do 100 i=1,nxt
      do 100 k=1,nzt
        cx=u(i,k)
        cz=w(i,k)
        xx = - cx*dt
        zz = - cz*dt
c
        isn=sign(1.0,cx)
        ksn=sign(1.0,cz)
        im1=i-isn
        km1=k-ksn
c
        a1=((gx(im1,k)+gz(i,k))*dx1*isn
     &     -2.0*(f(i,k)-f(im1,k)))/(dx3*isn)
        e1=(3.0*(f(im1,k)-f(i,k))
     &     +(gx(im1,k)+2.0*gx(i,k))*dx1*isn)/dx2
        b1=((gz(i,km1)+gz(i,k))*dz1*ksn
     &     -2.0*(f(i,k)-f(i,km1)))/(dz3*ksn)
        f1=(3.0*(f(i,km1)-f(i,k))
     &     +(gz(i,km1)+2.0*gz(i,k))*dz1*ksn)/dz2
c
        tmp=f(i,k)-f(i,km1)-f(im1,k)+f(im1,km1)
        tmq=gz(im1,k)-gz(i,k)
        d1= (-tmp -tmq*dz1*ksn)/(dx1*dz2*isn)
        c1=(-tmp-(gx(i,km1)-gx(i,k))*dx1*isn)/(dx2*dz1*ksn)
        g1=(-tmq+c1*dx2)/(dx1*isn)
c
        fn(i,k)=((a1*xx+c1*zz+e1)*xx+g1*zz+gx(i,k))*xx
     &      +((b1*zz+d1*xx+f1)*zz+gz(i,k))*zz+f(i,k)
        gxn(i,k)=(3.0*a1*xx+2.0*(c1*zz+e1))*xx+(d1*zz+g1)*zz+gx(i,k)
        gzn(i,k)=(3.0*b1*zz+2.0*(d1*xx+f1))*zz+(c1*xx+g1)*xx+gz(i,k)
c
  100 continue
c
      do i=1,nxt ! May/09/2005
        do k=1,nzt
         f(i,k)=fn(i,k)
        end do
      end do
c
      if(itp.eq.1) then 
       do i=1,nxt
         do k=1,nzt
          gx(i,k)=gxn(i,k)
          gz(i,k)=gzn(i,k)
         end do
       end do
c
	else if(itp.eq.3) then 
       do i=1,nxt
         do k=1,nzt
          gx(i,k)=gxn(i,k)
          gz(i,k)=gzn(i,k)
         end do
       end do
c
      else
       do i=1,nxt
         do k=1,nzt
          gx(i,k)=gxn(i,k)
          gz(i,k)=gzn(i,k)
         end do
       end do
c
      end if
c
      return
      end
c
c
c
c
c ========== subroutine for advection term calculation by CIP =====
      subroutine gr_cal(itp,dxi,dzet,fn,ux,uz,gx,gz)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k,itp,nxt,nzt
      integer im1,km1 
      real dxi,dzet
      real dtdx,dtdz
      real gxt,gzt
c     /* temporary valuables for cip routine */
      real fn(0:im,0:km),f(0:im,0:km)
      real u(0:im,0:km),w(0:im,0:km)
      real ux(0:im,0:km),uz(0:im,0:km)
      real gx(0:im,0:km),gz(0:im,0:km)
      real gxn(0:im,0:km),gzn(0:im,0:km)
c
      dtdx=dt/dxi
      dtdz=dt/dzet
c
c     /* velocity at u-point for ubn calculation */
      if (itp.eq.1) then
       nxt=nx-1
       nzt=nz
       do k=0,nz+1
         do i=0,nx
          u(i,k)=ux(i,k)
          w(i,k)=uz(i,k)
         end do
       end do
c
c     /* velocity at w-point for yrn calculation */
      else if (itp.eq.3) then
       nxt=nx
       nzt=nz-1
       do k=0,nz
         do i=0,nx+1
          u(i,k)=ux(i,k)
          w(i,k)=uz(i,k)
         end do
       end do
c
c     /* velocity at center-point for yrn calculation */
      else
       nxt=nx
       nzt=nz
       do k=1,nz
         do i=1,nx
          u(i,k)=0.5*(ux(i,k)+ux(i-1,k))
          w(i,k)=0.5*(uz(i,k)+uz(i,k-1))
         end do
       end do
c       
      end if
c
      if(itp.eq.1) then 
       do i=1,nxt
         do k=1,nzt
          gxn(i,k)=gx(i,k)
          gzn(i,k)=gz(i,k)
         end do
       end do
c
	else if(itp.eq.3) then 
       do i=1,nxt
         do k=1,nzt
          gxn(i,k)=gx(i,k)
          gzn(i,k)=gz(i,k) 
         end do
       end do
c
      else
       do i=1,nxt
         do k=1,nzt
          gxn(i,k)=gx(i,k)
          gzn(i,k)=gz(i,k)
         end do
       end do
      end if
c
      do i=1,nxt
        do k=1,nzt
         gxt=(fn(i+1,k)-fn(i-1,k))/2.0/dxi
         gzt=(fn(i,k+1)-fn(i,k-1))/2.0/dzet
         gx(i,k) = gxn(i,k)
     1               -(gxt*(u(i+1,k)-u(i-1,k))
     3                +gzt*(w(i+1,k)-w(i-1,k)))*0.5*dtdx
         gz(i,k) = gzn(i,k)
     1               -(gxt*(u(i,k+1)-u(i,k-1))
     3                +gzt*(w(i,k+1)-w(i,k-1)))*0.5*dtdz
        end do
      end do
c
      return
      end
c
c
c ========== subroutine for gradient calculation ==================
      subroutine newgrd(dxi,dzet,nxt,nzt,yn,y,gx,gz,nf,iff)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k,nxt,nzt,nf
      real dxi,dzet
      real yn(0:im,0:km),y(0:im,0:km)
     1    ,gx(0:im,0:km),gz(0:im,0:km)
c
c     masking flag
      integer iff(0:im,0:km)
c
c
      if(nf.eq.1) then 
       do i=1,nxt
         do k=1,nzt
          gx(i,k)=gx(i,k)+(yn(i+1,k)-yn(i-1,k)
     1                       -y(i+1,k)+y(i-1,k))*0.5/dxi
          gz(i,k)=gz(i,k)+(yn(i,k+1)-yn(i,k-1)
     1                       -y(i,k+1)+y(i,k-1))*0.5/dzet
         end do
       end do
       do i=0,nx
         do k=0,nz+1
          if(iff(i,k).eq.1) then
           gx(i,k)=0.0
           gz(i,k)=0.0
          end if
         end do
       end do
c
	else if(nf.eq.3) then 
       do i=1,nxt
         do k=1,nzt
          gx(i,k)=gx(i,k)+(yn(i+1,k)-yn(i-1,k)
     1                       -y(i+1,k)+y(i-1,k))*0.5/dxi
          gz(i,k)=gz(i,k)+(yn(i,k+1)-yn(i,k-1)
     1                       -y(i,k+1)+y(i,k-1))*0.5/dzet
         end do
       end do
       do i=1,nx
         do k=0,nz+1
          if(iff(i,k).eq.1) then
           gx(i,k)=0.0
           gz(i,k)=0.0
          end if
         end do
       end do
c
      else
       do i=1,nxt
         do k=1,nzt
          gx(i,k)=gx(i,k)+(yn(i+1,k)-yn(i-1,k)
     1                       -y(i+1,k)+y(i-1,k))*0.5/dxi
          gz(i,k)=gz(i,k)+(yn(i,k+1)-yn(i,k-1)
     1                       -y(i,k+1)+y(i,k-1))*0.5/dzet
         end do
       end do
c
       do i=1,nxt
         do k=1,nzt
          if(iff(i,k).eq.1) then
           gx(i,k)=0.0
           gz(i,k)=0.0
          end if
         end do
       end do
      end if

c
      return
      end
c
c1
c2
c3
c4
c5
c     ***************************************************************
c     *                                                             *   
c     *                      Shift Time Step                        *
c     *															  *
c     ***************************************************************
c
c ========== subroutine for shifting vluables =====================
      subroutine step_shift(ubn,ub,wbn,wb,ypn,yp,yrn,yr,x_xi,z_zet
     3                ,ub_l,wb_l,yp_l,yr_l,u_l,w_l,ifl)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real ubn(0:im,0:km),ub(0:im,0:km)
      real wbn(0:im,0:km),wb(0:im,0:km)
      real ypn(0:im,0:km),yp(0:im,0:km)
      real yrn(0:im,0:km),yr(0:im,0:km)
      real yp_l(0:im,0:km),yr_l(0:im,0:km)
      real u_l(0:im,0:km),w_l(0:im,0:km)
      real ub_l(0:im,0:km),wb_l(0:im,0:km)
      real x_xi(0:im,0:km),z_zet(0:im,0:km)
c
      real ypnt(0:im,0:km)
      real yrnt(0:im,0:km)
      integer ich
c     masking flag
      integer ifl(0:im,0:km)
c
      do i=0,nx
        do k=0,nz+1
         ub(i,k)=ubn(i,k)
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz
         wb(i,k)=wbn(i,k)
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz+1
         yp(i,k)=ypn(i,k)
         ypnt(i,k)=ypn(i,k)
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz+1
         yr(i,k)=yrn(i,k)
         yrnt(i,k)=yrn(i,k)
        end do
      end do
c
c     /* calculation for the velocity on decalto coordinates */
      do i=0,nx
        do k=0,nz
         ub_l(i,k)=(ub(i,k)+ub(i,k+1))*0.5
        end do
      end do
      do i=0,nx
        do k=0,nz
         wb_l(i,k)=(wb(i,k)+wb(i+1,k))*0.5
        end do
      end do
      do i=0,nx
        do k=0,nz
         ich=0
         ich=ifl(i,k)+ifl(i+1,k)+ifl(i,k+1)+ifl(i+1,k+1)
         if(ich.ne.0) then
          ub_l(i,k)=0.0
          wb_l(i,k)=0.0
         end if
        end do
      end do
c
c     setting value at the corner
	do i=0,nx
        do k=0,nz
         yp_l(i,k)=(ypnt(i,k)+ypnt(i,k+1)
     2               +ypnt(i+1,k)+ypnt(i+1,k+1))*0.25
         yr_l(i,k)=(yrnt(i,k)+yrnt(i,k+1)
     2               +yrnt(i+1,k)+yrnt(i+1,k+1))*0.25
        end do
      end do
c
	do i=0,nx
        do k=0,nz
         u_l(i,k)=(ub(i,k)+ub(i,k+1))*0.5
     2                     *(x_xi(i,k)+x_xi(i,k+1)
     5                      +x_xi(i+1,k)+x_xi(i+1,k+1))*0.25
        end do
      end do
	do i=0,nx
        do k=0,nz
         w_l(i,k)=(wb(i,k)+wb(i+1,k))*0.5
     2                     *(z_zet(i,k)+z_zet(i,k+1)
     3                      +z_zet(i+1,k)+z_zet(i+1,k+1))*0.25
        end do
      end do
c
      do i=0,nx
        do k=0,nz
         ich=0
         ich=ifl(i,k)+ifl(i+1,k)+ifl(i,k+1)+ifl(i+1,k+1)
         if(ich.ne.0) then
          u_l(i,k)=0.0
          w_l(i,k)=0.0
          yr_l(i,k)=yrnt(i,k)
         end if
        end do
      end do
c
      return
      end
c
c1
c2
c3
c4
c5
c     ***************************************************************
c     *                                                             *
c     *                      Short Subroutines                      *
c     *                                                             *
c     ***************************************************************
c
c ========== subroutine for velocity calculation at each point ====
      subroutine each_v(ub,wb,ub_w,wb_u) 
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real ub(0:im,0:km),wb(0:im,0:km)
      real ub_w(0:im,0:km)
      real wb_u(0:im,0:km)
c
      do i=0,nx
        do k=1,nz
         wb_u(i,k)=(wb(i,k)+wb(i+1,k)
     1                +wb(i,k-1)+wb(i+1,k-1))*0.25
        end do
      end do
c
      do i=1,nx
        do k=0,nz+1
         ub_w(i,k)=(ub(i,k)+ub(i-1,k)
     1                +ub(i,k+1)+ub(i-1,k+1))*0.25
        end do
      end do
c
      return
      end
c
c
c ========== subroutine for velocity boundary condition ===========
      subroutine bcset_v(ux,uz,ifu,ifw
     1     ,ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w
     2     ,bc1_ub,bc2_ub,bc3_wb,bc4_wb
     3     ,iscf,scu_ub,scw_wb,isu_u,isw_w)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real ux(0:im,0:km),uz(0:im,0:km)
      real bc1_ub(0:km),bc2_ub(0:km)
      real bc3_wb(0:im),bc4_wb(0:im)
c
c     masking flag
      integer ifu(0:im,0:km),ifw(0:im,0:km)  
c
c     b.c.set flag (ib1:upstream,ib2:downstream,ib3:top,ib4:bottom)
      integer ib1u,ib1w,ib2u,ib2w,ib3u,ib3w,ib4u,ib4w 
c
c     source setting
      integer iscf ! source flug 0:off,1:on
      integer isu_u(0:im,0:km),isw_w(0:im,0:km) ! source flag
      real scu_ub(0:im,0:km),scw_wb(0:im,0:km) ! source value
c
      if(ib1u.eq.0) then ! u=0 on left b.c.
          do k=0,nz+1
              ux(0,k)=0.0
          end do
      else if(ib1u.eq.1) then ! u given on left b.c.
          do k=0,nz+1
              ux(0,k)=bc1_ub(k)
              ux(1,k)=bc1_ub(k)
          end do
      else if(ib1u.eq.3) then
          do k=0,nz+1 ! periodic
            ux(0,k)=ux(nx-3,k)
c            ux(1,k)=ux(nx-2,k)
            ux(1,k)=bc1_ub(k)
          end do
      else 
          do k=0,nz+1 ! du/dx=0 on left b.c.
              ux(0,k)=ux(1,k)
          end do
      end if
      if(ib1w.eq.0) then ! w:non-slip on left b.c.
          do k=0,nz
              uz(0,k)=-uz(1,k)
          end do
      else
          do k=0,nz ! w:slip on left b.c.
              uz(0,k)=uz(1,k)
          end do
      end if
c
      if(ib2u.eq.0) then ! u=0 on right b.c.
          do k=0,nz+1
              ux(nx,k)=0.0
          end do
      else if(ib2u.eq.1) then ! u given on right b.c.
          do k=0,nz+1
              ux(nx,k)=bc2_ub(k)
          end do
      else if(ib2u.eq.3) then
          do k=0,nz+1 ! periodic
            ux(nx-1,k)=ux(2,k)
            ux(nx,k)=ux(3,k)
          end do
      else
          do k=0,nz+1 ! du/dx=0 on right b.c.
              ux(nx,k)=ux(nx-1,k)
          end do
      end if
      if(ib2w.eq.0) then ! w:non-slip on right b.c.
          do k=0,nz
              uz(nx+1,k)=-uz(nx,k)
          end do
      else
          do k=0,nz ! w:slip on right b.c.
              uz(nx+1,k)=uz(nx,k)
          end do
      end if    
c
      if(ib3w.eq.0) then ! w=0 on top b.c.
          do i=0,nx+1
              uz(i,nz)=0.0
          end do
      else if(ib3w.eq.1) then ! w given on top b.c.
          do i=0,nx+1
              uz(i,nz)=bc3_wb(i)
              uz(i,nz-1)=bc3_wb(i)
          end do
          do i=0,nx
              ux(i,nz)=0.0
          end do
      else
          do i=0,nx+1 ! dw/dz=0 on top b.c.
              uz(i,nz)=uz(i,nz-1)
          end do
      end if
      if(ib3u.eq.0) then ! u:non-slip on top b.c.
          do i=0,nx
              ux(i,nz+1)=-ux(i,nz)
          end do
      else
          do i=0,nx ! u:slip on top b.c.
              ux(i,nz+1)=ux(i,nz)
          end do
      end if    
c
      if(ib4w.eq.0) then ! w=0 on bottom b.c.
          do i=0,nx+1
              uz(i,0)=0.0
          end do
      else if(ib4w.eq.1) then ! w given on bottom b.c.
          do i=0,nx+1
              uz(i,0)=bc4_wb(i)
              uz(i,1)=bc4_wb(i)
c              uz(i,2)=bc4_wb(i)
          end do
          do i=0,nx
              ux(i,0)=0.0
          end do
      else
          do i=0,nx+1 ! dw/dz=0 on bottom b.c.
              uz(i,0)=uz(i,1)
          end do
      end if
      if(ib4u.eq.0) then ! u:non-slip on bottom b.c.
          do i=0,nx
              ux(i,0)=-ux(i,1)
          end do
      else
          do i=0,nx ! u:slip on bottom b.c.
              ux(i,0)=ux(i,1)
          end do
      end if    
c
c source setting
      if(iscf.eq.1) then
          do i=1,nx-1
              do k=1,nz
                  if(isu_u(i,k).eq.1)then
                      ux(i,k)=scu_ub(i,k)
                  end if
              end do
          end do
          do i=1,nx
              do k=1,nz-1
                  if(isw_w(i,k).eq.1)then
                      uz(i,k)=scw_wb(i,k)
                  end if
              end do
          end do
      end if
c
c masking u
      do i=0,nx 
        do k=0,nz+1
         if(ifu(i,k).eq.1) then
          ux(i,k)=0.0
         end if
        end do
      end do
      do i=0,nx 
        do k=0,nz
         if(ifu(i,k).eq.1.and.ifu(i,k+1).eq.0) then
             if(ifu(i-1,k).eq.1.and.ifu(i+1,k).eq.1) then
               ux(i,k)=-ux(i,k+1)
             end if
         end if
        end do
      end do
      do i=0,nx 
        do k=1,nz+1
         if(ifu(i,k).eq.1.and.ifu(i,k-1).eq.0) then
             if(ifu(i-1,k).eq.1.and.ifu(i+1,k).eq.1) then
                ux(i,k)=-ux(i,k-1)
             end if
         end if
        end do
      end do
c
c masking w
      do i=1,nx 
        do k=0,nz+1
         if(ifw(i,k).eq.1) then
          uz(i,k)=0.0
         end if
        end do
      end do
      do i=1,nx 
        do k=0,nz
         if(ifw(i,k).eq.1.and.ifw(i+1,k).eq.0) then
             if(ifw(i,k-1).eq.1.and.ifw(i,k+1).eq.1) then
               uz(i,k)=-uz(i+1,k)
             end if
         end if
        end do
      end do
      do i=1,nx 
        do k=1,nz+1
         if(ifw(i,k).eq.1.and.ifw(i-1,k).eq.0) then
             if(ifw(i,k-1).eq.1.and.ifw(i,k+1).eq.1) then
               uz(i,k)=-uz(i-1,k)
             end if
         end if
        end do
      end do
c
      return
      end
c
c
c ========== subroutine for velocity boundary condition ===========
      subroutine bcset_r(yr,ifr
     1     ,ib1r,ib2r,ib3r,ib4r
     2     ,bc1_r,bc2_r,bc3_r,bc4_r,rho
     3     ,iscf,scr_c,isr_c)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real yr(0:im,0:km)
      real bc1_r(0:km),bc2_r(0:km)
      real bc3_r(0:im),bc4_r(0:im)
      real rho
c
c     masking flag
      integer ifr(0:im,0:km)  
c
c     b.c.set flag (ib1:upstream,ib2:downstream,ib3:top,ib4:bottom)
      integer ib1r,ib2r,ib3r,ib4r  
c
c     source setting
      integer iscf ! source flug 0:off,1:on
      integer isr_c(0:im,0:km) ! source flag
      real scr_c(0:im,0:km) ! source value
c
      if(ib1r.eq.0) then ! dr/dx =0 on left b.c.
          do k=1,nz
              yr(0,k)=yr(1,k)
          end do
          do i=0,nx+1
              yr(i,nz+1)=yr(i,nz)
              yr(i,0)=yr(i,1)
          end do
      else if(ib1r.eq.3) then  ! r periodic
          do k=1,nz
              yr(0,k)=yr(nx-3,k)
              yr(1,k)=yr(nx-2,k)
          end do
          do i=0,nx+1
              yr(i,nz+1)=yr(i,nz)
              yr(i,0)=yr(i,1)
          end do
      else   ! r given on left b.c.
          do k=1,nz
              yr(0,k)=bc1_r(k)
              yr(1,k)=bc1_r(k)
          end do
          do i=0,nx+1
              yr(i,nz+1)=yr(i,nz)
              yr(i,0)=yr(i,1)
          end do
      end if      
c
      if(ib2r.eq.0) then ! dr/dx =0 on on right b.c.
          do k=1,nz
              yr(nx+1,k)=yr(nx,k)
          end do
          do i=0,nx+1
              yr(i,nz+1)=yr(i,nz)
              yr(i,0)=yr(i,1)
          end do
      else if(ib2r.eq.3) then  ! r periodic
          do k=1,nz
              yr(nx-1,k)=yr(2,k)
              yr(nx,k)=yr(3,k)
              yr(nx+1,k)=yr(4,k)
          end do
          do i=0,nx+1
              yr(i,nz+1)=yr(i,nz)
              yr(i,0)=yr(i,1)
          end do
      else  ! r given on right b.c.
          do k=1,nz
              yr(nx+1,k)=bc2_r(k)
              yr(nx,k)=bc2_r(k)
          end do
          do i=0,nx+1
              yr(i,nz+1)=yr(i,nz)
              yr(i,0)=yr(i,1)
          end do
      end if
c
      if(ib3r.eq.0) then ! dr/dx =0 on on top b.c.
          do i=1,nx
              yr(i,nz+1)=yr(i,nz)
          end do
          do k=0,nz+1
              yr(0,k)=yr(1,k)
              yr(nx+1,k)=yr(nx,k)
          end do
      else  ! r given on top b.c.
          do i=1,nx
              yr(i,nz+1)=bc3_r(i)
              yr(i,nz)=bc3_r(i)
          end do
          do k=0,nz+1
              yr(0,k)=yr(1,k)
              yr(nx+1,k)=yr(nx,k)
          end do
      end if   
c
      if(ib4r.eq.0) then ! dr/dx =0 on on bottom b.c.
          do i=1,nx
              yr(i,0)=yr(i,1)
          end do
          do k=0,nz+1
              yr(0,k)=yr(1,k)
              yr(nx+1,k)=yr(nx,k)
          end do
      else  ! r given on bottom b.c.
          do i=1,nx
              yr(i,1)=bc4_r(i)
              yr(i,0)=bc4_r(i)
          end do
          do k=0,nz+1
              yr(0,k)=yr(1,k)
              yr(nx+1,k)=yr(nx,k)
          end do
      end if     
c
c source
      if(iscf.eq.1) then
       do i=1,nx 
         do k=1,nz
          if(isr_c(i,k).eq.1) then
           yr(i,k)=scr_c(i,k)
          end if
         end do
       end do
      end if
c
c masking r
      do i=0,nx+1 
        do k=0,nz+1
         if(ifr(i,k).eq.1) then
          yr(i,k)=rho
         end if
        end do
      end do
      do i=1,nx+1 
        do k=0,nz+1
         if(ifr(i,k).eq.1.and.ifr(i-1,k).eq.0) then
               yr(i,k)=yr(i-1,k)
         end if
        end do
      end do
      do i=0,nx 
        do k=0,nz+1
         if(ifr(i,k).eq.1.and.ifr(i+1,k).eq.0) then
               yr(i,k)=yr(i+1,k)
         end if
        end do
      end do
      do i=0,nx+1 
        do k=0,nz
         if(ifr(i,k).eq.1.and.ifr(i,k+1).eq.0) then
               yr(i,k)=yr(i,k+1)
         end if
        end do
      end do
      do i=0,nx+1 
        do k=1,nz+1
         if(ifr(i,k).eq.1.and.ifr(i,k-1).eq.0) then
               yr(i,k)=yr(i,k-1)
         end if
        end do
      end do
c
      return
      end
c
c ========== subroutine for gradient boundary condition ===========
      subroutine bcset_g(itp,ygx,ygz,iff,ib1u,ib2u)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k,itp
      real ygx(0:im,0:km),ygz(0:im,0:km)
c
c     masking flag
      integer iff(0:im,0:km)
      integer ib1u,ib2u
c
c     /* for gux,guz */
      if(itp.eq.1) then 
          do k=0,nz+1
             ygx(0,k)=0.0
             ygx(nx,k)=0.0
             ygx(nx+1,k)=0.0
             ygz(0,k)=0.0
             ygz(nx,k)=0.0
            ygz(nx+1,k)=0.0
          end do
          if(ib1u.eq.3) then
           do k=0,nz+1
             ygx(0,k)=ygx(nx-3,k)
             ygx(1,k)=ygx(nx-2,k)
             ygz(0,k)=ygz(nx-3,k)
             ygz(1,k)=ygz(nx-2,k)
           end do
          end if
          if(ib2u.eq.3) then
           do k=0,nz+1
             ygx(nx-1,k)=ygx(2,k)
             ygx(nx,k)=ygx(3,k)
             ygz(nx-1,k)=ygz(2,k)
             ygz(nx,k)=ygz(3,k)
           end do
          end if
           do i=0,nx
             ygx(i,0)=0.0
             ygx(i,nz+1)=0.0
             ygz(i,0)=0.0
             ygz(i,nz+1)=0.0
           end do
c
c     /* for gwx,gwz */
      else if(itp.eq.3) then
       do k=0,nz+1
         ygx(0,k)=0.0
         ygx(nx+1,k)=0.0
         ygz(0,k)=0.0
         ygz(nx+1,k)=0.0
       end do
       do i=0,nx+1
         ygx(i,0)=0.0
         ygx(i,nz)=0.0
         ygx(i,nz+1)=0.0
         ygz(i,0)=ygz(i,1) !
c         ygz(i,nz)=ygz(i,nz-1)
         ygz(i,nz)=0.0
         ygz(i,nz+1)=0.0
       end do
c
      else
c     /* for grx,grz */
       do k=0,nz+1
         ygx(0,k)=0.0
         ygx(nx+1,k)=0.0
         ygz(0,k)=0.0
         ygz(nx+1,k)=0.0
       end do
       do i=0,nx+1
         ygx(i,0)=0.0
         ygx(i,nz+1)=0.0
         ygz(i,0)=0.0
         ygz(i,nz+1)=0.0
       end do
      end if
c
       do i=0,nx
         do k=0,nz
          if(iff(i,k).ne.0) then 
           ygx(i,k)=0.0
           ygz(i,k)=0.0
          end if
         end do
       end do
c
      return
      end
c
c1
c2
c3
c4
c5
c ========== subroutine for shifting vluables =====================
      subroutine timerelax(ubn,ub,wbn,wb,ypn,yp,yrn,yr,rctime)
c =================================================================
c
      implicit none
      include 'common.h'
c
      integer i,k
      real ubn(0:im,0:km)
      real wbn(0:im,0:km)
      real ypn(0:im,0:km)
      real yrn(0:im,0:km)
      real ub(0:im,0:km)
      real wb(0:im,0:km)
      real yp(0:im,0:km)
      real yr(0:im,0:km)
c
      real rctime
c
      do i=0,nx
        do k=0,nz+1
         ubn(i,k)=ubn(i,k)*(1.0-rctime)+ub(i,k)*rctime
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz
         wbn(i,k)=wbn(i,k)*(1.0-rctime)+wb(i,k)*rctime
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz
         ypn(i,k)=ypn(i,k)*(1.0-rctime)+yp(i,k)*rctime
        end do
      end do
c
      do i=0,nx+1
        do k=0,nz
         yrn(i,k)=yrn(i,k)*(1.0-rctime)+yr(i,k)*rctime
        end do
      end do
c
      return
      end
c
c
c
c
