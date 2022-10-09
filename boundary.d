0,0,0,0 ! nonslip,u=0 
1,2,0,0 ! top slip&open, bottom given,other nonslip,w=0
0,0,1,0 ! top open, dp/dx=0, dp/dz=0
0,0,0,0 ! bottom given, other dr/dx=0, dr/dz=0



      read(12,*) ib1u,ib1w,ib2u,ib2w
      read(12,*) ib3u,ib3w,ib4u,ib4w
      read(12,*) ib1p,ib2p,ib3p,ib4p
      read(12,*) ib1r,ib2r,ib3r,ib4r
