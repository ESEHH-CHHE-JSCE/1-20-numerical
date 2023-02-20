      program csvdata
c
      implicit none
      include 'common.h'
c
      integer i,k,j,ny,istopn
      real x(0:im,0:km),z(0:im,0:km)
      integer ifl(0:im,0:km)
      integer iscf
c
      real u_l(0:im,0:km),w_l(0:im,0:km)
      real yp_l(0:im,0:km)
      real yr_l(0:im,0:km)
      real bc1_u(0:km),bc2_u(0:km) ! b.c.flow at U part on left & right
      real bc3_w(0:im),bc4_w(0:im) ! b.c.flow at W part on top & bottom
      real bc1_r(0:km),bc2_r(0:km) ! b.c.rho at center on left & right
      real bc3_r(0:im),bc4_r(0:im) ! b.c.rho at center on top & bottom
c
      real chleng,height
      real dx,dz
c
      real g,rho,rnu,diffyr
      real tuk,etime
      real lsor,soralp
      real pretime
c
      integer m,ixs,ixe,namecount
      integer id,kd
      integer icount,istep
      character*13 buff1
      character*13 buff2
      character*13 buff5
c
      ixs=100
      ixe=300
c
      write(*,*) 'Now reading'
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
      icount=1
      istep=1
c
      istopn=int(etime/tuk)
c      
c      icontrol=1 !output timing control
      do m=1,istopn
        namecount=icount-1
c
c 
        if(namecount.lt.10) then
          write(buff5,505) namecount
        else if((namecount.gt.9).and.(namecount.lt.100)) then
          write(buff5,506) namecount
        else if((namecount.gt.99).and.(namecount.lt.1000)) then
          write(buff5,507) namecount
        else if((namecount.gt.999).and.(namecount.lt.10000)) then
          write(buff5,508) namecount
        else
          write(buff5,509) namecount
        end if
        write(6,*) buff5
c
       open(56,file=buff5,status='old',form='unformatted')
c
       read(56,end=888) time
       read(56,end=888) id,kd
       read(56,end=888) ((x(i,k),i=0,nx),k=0,nz)
       read(56,end=888) ((z(i,k),i=0,nx),k=0,nz)
       read(56,end=888) ((u_l(i,k),i=0,nx),k=0,nz)
       read(56,end=888) ((w_l(i,k),i=0,nx),k=0,nz)
       read(56,end=888) ((ifl(i,k),i=0,nx),k=0,nz)
       read(56,end=888) ((yp_l(i,k),i=0,nx),k=0,nz)
       read(56,end=888) ((yr_l(i,k),i=0,nx),k=0,nz)
c
       close(56)
c
        if(namecount.lt.10) then
          write(buff1,205) namecount
          write(buff2,305) namecount
        else if((namecount.gt.9).and.(namecount.lt.100)) then
          write(buff1,206) namecount
          write(buff2,306) namecount
        else if((namecount.gt.99).and.(namecount.lt.1000)) then
          write(buff1,207) namecount
          write(buff2,307) namecount
        else if((namecount.gt.999).and.(namecount.lt.10000)) then
          write(buff1,208) namecount
          write(buff2,308) namecount
        else
          write(buff1,209) namecount
          write(buff2,309) namecount
        end if
        write(6,*) buff1
c
        open(26,file=buff1,status='unknown',form='formatted')
        write(26,*)'i',',','k',',','x',',','y',',','concentration'
        do i=ixs,ixe
         do k=0,nz
           write(26,601) i,k,x(i,k),z(i,k),yr_l(i,k)
         end do
        end do 
        close(26)
c
        open(27,file=buff2,status='unknown',form='formatted')
        write(27,*)'i',',','k',',','x',',','y',',','u',',','w'
        do i=ixs,ixe
         do k=0,nz
           write(27,701) i,k,x(i,k),z(i,k),u_l(i,k),w_l(i,k)
         end do
        end do 
        close(27)
c
 888  continue
c
	 icount=icount+1
c
c
      end do
 999  continue
c
 205  format('cntr0000',i1,'.csv') 
 206  format('cntr000',i2,'.csv')
 207  format('cntr00',i3,'.csv') 
 208  format('cntr0',i4,'.csv')  
 209  format('cntr',i5,'.csv')
 305  format('vctr0000',i1,'.csv') 
 306  format('vctr000',i2,'.csv')
 307  format('vctr00',i3,'.csv') 
 308  format('vctr0',i4,'.csv')  
 309  format('vctr',i5,'.csv')
 505  format('rslt0000',i1,'.dat') 
 506  format('rslt000',i2,'.dat')
 507  format('rslt00',i3,'.dat') 
 508  format('rslt0',i4,'.dat')  
 509  format('rslt',i5,'.dat')
 601  format(i4,',',i4,','f8.4,','f8.4,','f8.4)
 701  format(i4,',',i4,','f8.4,','f8.4,','f8.4,','f8.4)
c
      end