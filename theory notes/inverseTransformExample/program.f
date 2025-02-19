      implicit none
      integer ndat,nbox
      parameter(ndat=1000000)
      parameter(nbox=int(sqrt(dble(ndat))))
      double precision canviVar,x,r,xnums(ndat),xhis(nbox),vhis(nbox),
     + ERRHIS(nbox),boxsize,a,b
      integer i
      external canviVar
      call srand(123)

      open(1,file='res.dat')
      do i=1,ndat
      	r=rand()

      	xnums(i)=canviVar(r)
      	
      enddo
      
      
      a=0d0
      b=dacos(-1d0)


      call HISTOGRAMA(ndat,xnums,a,b,NBOX,XHIS,VHIS,ERRHIS
     + ,BOXSIZE)

      do i=1,nbox
      	write(1,*) xhis(i),vhis(i),ERRHIS(i)
      enddo
      
      end


      double precision function canviVar(x)
      implicit none
      double precision x
      canviVar=dacos(1d0-2d0*x)
      end

      SUBROUTINE HISTOGRAMA(NDAT,XDATA,XA,XB,NBOX,XHIS,VHIS,ERRHIS
     + ,BOXSIZE)
       IMPLICIT NONE

C INPUT/OUTPUT VARIABLES

       INTEGER NDAT,NBOX
       DOUBLE PRECISION XDATA(NDAT),XA,XB
       DOUBLE PRECISION XHIS(NBOX),VHIS(NBOX),ERRHIS(NBOX)
       INTEGER IERR
C 
       INTEGER I,IBOX,ICOUNT
       DOUBLE PRECISION BOXSIZE

       IF (XA.GE.XB) THEN 
          IERR=1
          RETURN
       ENDIF
C BOX SIZE
         BOXSIZE=(XB-XA)/NBOX

C COUNTS NUMBER OF POINTS WITHIN THE INTERVAL XA,XB
       ICOUNT=0

C SETS ALL TO ZERO
       DO I=1,NBOX
          VHIS(I)=0
          ERRHIS(I)=0
       ENDDO

C WE RUN THROUGH THE DATASET
       DO I=1,NDAT
C CHECKS IF DATA LIES WITHIN XA,XB
         IF (XDATA(I).GE.XA.AND.XDATA(I).LE.XB) THEN 
            IBOX=INT((XDATA(I)-XA)/BOXSIZE)+1
C PUTS XB INTO THE LAST BOX, IF NEEDED
            IF (IBOX.EQ.NBOX+1) IBOX=NBOX 

            VHIS(IBOX)=VHIS(IBOX)+1
            ICOUNT=ICOUNT+1
         ENDIF
        ENDDO

        IF (ICOUNT.EQ.0) THEN 
           IERR=2
           RETURN
        ENDIF

        IERR=0
        PRINT*,"ACCEPTED:",ICOUNT,
     1  " OUT OF:",NDAT

        DO I=1,NBOX
c CENTRAL VALUE OF THE BAR
           XHIS(I)=XA+BOXSIZE/2.D0+(I-1)*BOXSIZE
C  ERROBAR, STANDARD DEVIATION OF CORRESPONDING BINOMIAL
           ERRHIS(I)=SQRT(VHIS(I)/ICOUNT*(1.D0-VHIS(I)/ICOUNT))
     1      /BOXSIZE / SQRT(DBLE(ICOUNT))
C NORMALIZED VALUE OF THE BAR
           VHIS(I)=VHIS(I)/ICOUNT/BOXSIZE
        ENDDO
        END
