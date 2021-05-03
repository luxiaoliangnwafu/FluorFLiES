

      implicit none

      real s,x,y,z,ux,uy,uz,x1,y1,z1
      real int(3)
      real th,ph,pi,tobj(5)
      integer io, face
      
      pi=3.14159265358979

      int(1)=30.0
      int(2)=30.0
      int(3)=10.0

      x1=0.0
      y1=0.0
      z1=0.0
d     comment out


      x=8.
      y=15.
      z=8.0

      th=35.0
      ph=0.0
      th=th*pi/180.0
      ph=ph*pi/180.0
      
      ux=sin(th)*cos(ph)
      uy=sin(th)*sin(ph)
      uz=cos(th)

      tobj(1)=15.
      tobj(2)=15.
      tobj(3)=10.
      tobj(4)=5.
c     tobj(4)=tobj(4)/2.
      tobj(5)=5.

      face=-1

c      call planes(s,x,y,z,ux,uy,uz,x1,y1,z1,face,int)
c     call  elpss(s ,x ,y ,z ,ux ,uy ,uz ,tobj, io)
c     call  cyls(s, x, y, z, ux, uy, uz, tobj, face, io)
c      call  cones(s,x,y,z,ux,uy,uz,tobj,face,io)

      call  helpss(s ,x ,y ,z ,ux ,uy ,uz ,tobj, face,io)

      write(*,*) s,io,face
      
      stop
      end
