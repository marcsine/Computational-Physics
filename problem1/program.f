      program metodetir
      implicit none
      integer i,i3,nmax
      parameter(nmax=100000)
      double precision xxg,xxm,xxb,pi,v0,x0,y0,theta,cond(4),dt,t
      double precision yyin(4),yyout(4),x3(nmax),y3(nmax),v01,v02
      common/dades/xxg,xxm,xxb

      xxg=9.8d0
      xxm=0.625d0
      xxb=0.2d0
      pi=4.d0*datan(1.d0)

      x0=-6.75d0
      y0=2.d0
      theta=(2d0/9d0)*pi
      

      dt=0.005d0

      open(11,file="res11.dat")

      v0=5.d0

      call condinicials(v0,theta,x0,y0,cond)

      yyin=cond
      yyout=yyin

      i=0
      do while((yyout(1).lt.0).and.(yyout(3).gt.0))
      	i=i+1
      	t=0.d0+dble(i)*dt
      	call mirungekutta4(4,t,dt,yyin,yyout)
      	write(11,*) t,yyout(1),yyout(3)
      	yyin=yyout
      enddo
      write(11,*) ""
      write(11,*) ""

      v0=7.5d0

      call condinicials(v0,theta,x0,y0,cond)

      yyin=cond
      yyout=yyin

      i=0
      do while((yyout(1).lt.0).and.(yyout(3).gt.0))
        i=i+1
        t=0.d0+dble(i)*dt
        call mirungekutta4(4,t,dt,yyin,yyout)
        write(11,*) t,yyout(1),yyout(3)
        yyin=yyout
      enddo
      write(11,*) ""
      write(11,*) ""

      v0=10.d0
      
      call condinicials(v0,theta,x0,y0,cond)

      yyin=cond
      yyout=yyin

      i=0
      do while((yyout(1).lt.0).and.(yyout(3).gt.0))
        i=i+1
        t=0.d0+dble(i)*dt
        call mirungekutta4(4,t,dt,yyin,yyout)
        write(11,*) t,yyout(1),yyout(3)
        yyin=yyout
      enddo
      write(11,*) ""
      write(11,*) ""

      v0=12.5d0

      call condinicials(v0,theta,x0,y0,cond)

      yyin=cond
      yyout=yyin
      
      i=0
      do while((yyout(1).lt.0).and.(yyout(3).gt.0))
        i=i+1
        t=0.d0+dble(i)*dt
        call mirungekutta4(4,t,dt,yyin,yyout)
        write(11,*) t,yyout(1),yyout(3)
        yyin=yyout
      enddo
      write(11,*) ""
      write(11,*) ""

      v0=15.d0

      call condinicials(v0,theta,x0,y0,cond)

      yyin=cond
      yyout=yyin

      i=0
      do while((yyout(1).lt.0).and.(yyout(3).gt.0))
        i=i+1
        t=0.d0+dble(i)*dt
        call mirungekutta4(4,t,dt,yyin,yyout)
        write(11,*) t,yyout(1),yyout(3)
        yyin=yyout
      enddo
      write(11,*) ""
      write(11,*) ""

      close(11)

      

c-----------------------------------------------------------------------
      open(22,file="res22.dat")

      theta=(2d0/9d0)*pi
      v01=10.d0
      v02=12.5d0

      call disparo(theta,v01,v02,x3,y3,i3)
 
      do i=1,i3
        write(22,*) x3(i),y3(i)
      enddo


      end program

c.......................................................................
      subroutine condinicials(v0,theta,x0,y0,cond)
      implicit none
      double precision v0,x0,y0,theta,cond(4)
      cond(1)=x0
      cond(2)=v0*dcos(theta)
      cond(3)=y0
      cond(4)=v0*dsin(theta)
      end subroutine


      subroutine derivades(nequ,x,yin,dyout)
c yin(4)=(x,\cdot{x},y,\cdot{y})
      implicit none
      integer nequ
      double precision x,yin(nequ),dyout(nequ),xxg,xxm,xxb
      common/dades/xxg,xxm,xxb
      dyout(1)=yin(2) !\cdot{x}
      dyout(2)=(-xxb/xxm)*yin(2) !\cdot\cdot{x}
      dyout(3)=yin(4) !\cdot{y}
      dyout(4)=(-xxb/xxm)*yin(4)-xxg !\cdot\cdot{y}
      return
      end subroutine

     
      subroutine mirungekutta4(nequs,x,dx,yyin,yyout)
      implicit none
      integer nequs
      double precision yyin(nequs),yyout(nequs),x,x2,x3,x4,dx
      double precision k1(nequs),k2(nequs),k3(nequs),k4(nequs)

      call derivades(nequs,x,yyin,k1) !k1

      x2=x+dx/2.d0 
      call derivades(nequs,x2,yyin+dx*k1/2.d0,k2) !k2

      x3=x+dx/2.d0 
      call derivades(nequs,x3,yyin+dx*k2/2.d0,k3) !k3

      x4=x+dx
      call derivades(nequs,x4,yyin+dx*k3,k4) !k4

      yyout=yyin+dx/6d0*(k1+2.d0*k2+2.d0*k3+k4) !y1
      
      return
      end subroutine


      subroutine disparo(theta,v01,v02,x3,y3,i3)
      implicit none
      integer i1,i2,i3,nmax
      parameter(nmax=100000)
      double precision x0,y0,dt,t,d1,d2,d3,v01,v02,v03,theta,cond(4)
      double precision yyin(4),yyout(4),x1(nmax),y1(nmax),x2(nmax)
      double precision y2(nmax),x3(nmax),y3(nmax)

      x0=-6.75d0
      y0=2.d0
      dt=0.005d0

      d3=1000d0

      do while (d3.gt.0.01d0)

        call condinicials(v01,theta,x0,y0,cond)
        yyin=cond
        yyout=cond
        i1=0
        do while (yyout(1).lt.0.d0)
          i1=i1+1
          t=0.d0+dble(i1)*dt
          call mirungekutta4(4,t,dt,yyin,yyout)
          x1(i1)=yyout(1)
          y1(i1)=yyout(3)
          yyin=yyout
        enddo
        d1=dabs(3.05d0-y1(i1))

        call condinicials(v02,theta,x0,y0,cond)
        yyin=cond
        yyout=cond
        i2=0
        do while (yyout(1).lt.0.d0)
          i2=i2+1
          t=0.d0+dble(i2)*dt
          call mirungekutta4(4,t,dt,yyin,yyout)
          x2(i2)=yyout(1)
          y2(i2)=yyout(3)
          yyin=yyout
        enddo
        d2=dabs(3.05d0-y2(i2))

        v03=dabs((v01*d2-v02*d1)/(d2-d1))
        call condinicials(v03,theta,x0,y0,cond)
        yyin=cond
        yyout=cond
        i3=0
        do while (yyout(1).lt.0.d0)
          i3=i3+1
          t=0.d0+dble(i3)*dt
          call mirungekutta4(4,t,dt,yyin,yyout)
          x3(i3)=yyout(1)
          y3(i3)=yyout(3)
          yyin=yyout
        enddo
        d3=dabs(3.05d0-y3(i3))

        v01=v02
        v02=v03

      enddo
      print*,'la resposta a apartat 2 es: ',v03
      return
      end subroutine






      