C*******************************************************************C
C*******************************************************************C
C Monte-Carlo Simulation In Honeycomb Lattice With Heisenberg Spin  C
C                Developing By Nepal Banerjee                       C
C                University of Seoul,South-Korea                    C
C                Email:nb.uos1989@gmail.com                         C
C                     Date:28/11/2025                               C
C*******************************************************************C
C*******************************************************************C
          IMPLICIT NONE 
          INTEGER L,LSQ,itmx,itequi,itav,LMAX
          REAL*8 pi,th,phi,S
          PARAMETER(L=30,LSQ=L*L,itmx=900000,itequi=500000,
     &    itav=itmx-itequi,pi=4*atan(1.d0),LMAX=L) !LMAX=>Segment Fault
          INTEGER i,j,temp,tmax,tmin,dt 
          INTEGER ic,id,ip,in,jp,jn,ipp,inn,jl,j3NNNN
          REAL*8 X1J1,X2J1,X3J1,X1J2,X2J2,X3J2,X1J3,X2J3,X3J3
          REAL*8 KX,KY,KZ,GamaX,GamaY,GamaZ
          REAL*8 sx(LMAX,LMAX),sy(LMAX,LMAX),sz(LMAX,LMAX)
          REAL*8 flipx(LMAX,LMAX),flipy(LMAX,LMAX),flipz(LMAX,LMAX)
          REAL*8 t,E,E_flip,dE,prob,am1,am2,am3,energy1,energy2
          REAL*8 xm,xm1,xm2,xm3,en,flucM,flucE,UL
          REAL*8 r1,r2,r
C===================================================C
C         DATA CARD                                 C
C===================================================C
          OPEN(1,file='KGH.ini',status='old')
            read(1,*) S                  !SPIN MAG
            read(1,*) X1J1,X2J1,X3J1     !J1
            read(1,*) X1J2,X2J2,X3J2     !J2
            read(1,*) X1J3,X2J3,X3J3     !J3
            read(1,*) KX,KY,KZ           !KX,KY,KZ
            read(1,*) GamaX,GamaY,GamaZ  !GammaX,GammaY,GammaZ 
          CLOSE(1)

C*****************************************
         CALL SYSTEM('rm HSO')
         call SYSTEM('rm spin_config_ini')

         OPEN(4,file='HSO',status='new')
           WRITE(4,*)'#### T #### M #### Chi #### E #### Cv #### U ##'
         CLOSE(4)         

         OPEN(5,file='spin_config_ini',status='new')
           WRITE(5,*) ' ## i ##  j ## Sx ## Sy ## Sz ## '
         CLOSE(5)


C**********************************C
C         Screen Output for User   C
C**********************************C
           WRITE(*,29)
           WRITE(*,31) S,X1J1,X2J1,X3J1,X1J2,X2J2,X3J2,X1J3,X2J3,
     &X3J3,KX,KY,KZ,GamaX,GamaY,GamaZ,LSQ,itmx,itequi
           WRITE(*,22)

C----------------------------------C
C        Format                    C
C----------------------------------C
   22    Format(/,99('*'),/,10X,'Temp',15X,'Magnetization',15X,'Energy',
     &   20X,'Heat capacity',/,99('-'))

   27     Format(5X,'#Temp',23X,'Magnetization',20X,'Susceptibility',20X
     &   'Energy',20X,'Cv')


   29     Format(/,5X,50('*'),/,5X,50('*'),/,30X,'Hi !!',/,
     &25X,'Welcome to',/,
     &20X,'Classical Monte-Carlo Simulation With',/,
     &20X,'Heisenberg  Spin in a Honeycomb Lattice',/,
     &20X,'In Presence Of J1(NN) ,J2(NNN) ,J3(NNNN) ',/,
     &20X,'KX,KY,KZ,GammaX,GammaY,GammaZ',/,5X,50('*'),/,
     &20X,'Developed by: Nepal Banerjee',/,
     &20X,'Email Id:nb.uos1989@gmail.com',/,      
     &20X,'University of Seoul,Seoul,South-Korea ',/,
     &25X,'Date:25/11/2025',/,5X,50('*'),
     &/,5X,50('*'))
     

  31     Format(/,5X,50('-'),/,
     &20X,'Spin Mag:',F5.2,/,
     &20X,'JX_NN,JY_NN,JZ_NN:',F5.2,',',F5.2,',',F5.2/,
     &20X,'JX_NNN,JY_NNN,JZ_NNN:',F5.2,',',F5.2,',',F5.2/,
     &20X,'JX_NNNN,JY_NNNN,JZ_NNNN:',F5.2,',',F5.2,',',F5.2/,
     &20X,'KX,KY,KZ:',F5.2,',',F5.2,',',F5.2/,
     &20X,'GammaX,GammaY,GammaZ:',F5.2,',',F5.2,',',F5.2/,
     &20X,'Lattice Size:',I10,/,
     &20X,'Maximum MC steps:',I10,/,
     &20x,'Minimum MC steps for Equillibrium:',I10,/,
     &5X,50('-'))

  32     Format(/,10X,60('*'),/,
     &30X,'Congratulation ! ',/,
     &20X,'Simulation has completed Successfully. ',/,
     &20X,'Thanks for your attention,interest and patience.',/,
     &20X,'Please collect all the results from data files.',/,
     &20X,'Your valuable feed-back and suggestions are welcome.',/,
     &30X,'Wish you  a great Day !',/,10X,60('*'))

C******************************************C
C*** Simulating Initial Spin Config  ******C
C******************************************C
         DO i=1,L
           DO j=1,L
               CALL RANDOM_NUMBER(r1)
               CALL RANDOM_NUMBER(r2)
               th=acos(2*r1-1)
               phi= 2*pi*r2
               sx(i,j)=S*sin(th)*cos(phi)  ! Sx
               sy(i,j)=S*sin(th)*sin(phi)  ! Sy
               sz(i,j)=S*cos(th)           ! Sz
           ENDDO
         ENDDO  
C******************************************
C******************************************

      OPEN(5, file='spin_config_ini', status='unknown')
        DO i = 1, L
           DO j = 1, L
              WRITE(5,*) i, j, sx(i,j), sy(i,j), sz(i,j)
           ENDDO
        ENDDO
      CLOSE(5)

C******** Main Simulation Loop Start Here *******
        tmax=50         ;          ! Maximum Temperature
        tmin=1          ;          ! Minimum Temperature
        dt=1            ;          ! Temperature Step

        DO temp=tmax,tmin,-dt      ! Temperature Loop 
           t=temp/100.0 ;

           am1=0.d0     ;
           am2=0.d0     ;
           am3=0.d0     ;
           energy1=0.d0 ;
           energy2=0.d0 ;

            DO ic=1,itmx           ! Iteration  Loop
      
            DO id=1,LSQ
            CALL RANDOM_NUMBER(r1)
            CALL RANDOM_NUMBER(r2)
      
            i=1+(r1*L)
            j=1+(r2*L)

            call Periodic_Boundary_Condition(i,j,ip,in,jp,jn,ipp,inn,
     &jl,j3NNNN,L)

C**** Calculate the energy of the Chosen spin *************C

        E=-sx(i,j)*X1J1*(sx(ip,j)+sx(in,j)+sx(i,jl))  !J1(NN)
     &  -sy(i,j)*X2J1*(sy(ip,j)+sy(in,j)+sy(i,jl))    !J1(NN)
     &  -sz(i,j)*X3J1*(sz(ip,j)+sz(in,j)+sz(i,jl))    !J1(NN)

     &  -sx(i,j)*X1J2*(sx(ip,jp)+sx(ip,jn)+sx(in,jn)+sx(in,jp))!J2(NNN)
     &  -sy(i,j)*X2J2*(sy(ip,jp)+sy(ip,jn)+sy(in,jn)+sy(in,jp))!J2(NNN)
     &  -sz(i,j)*X3J2*(sz(ip,jp)+sz(ip,jn)+sz(in,jn)+sz(in,jp))!J2(NNN)

     &  -sx(i,j)*X1J2*(sx(ipp,j)+sx(inn,j)) !J2(NNN)
     &  -sy(i,j)*X2J2*(sy(ipp,j)+sy(inn,j)) !J2(NNN)
     &  -sz(i,j)*X3J2*(sz(ipp,j)+sz(inn,j)) !J2(NNN)

     &  -sx(i,j)*X1J3*(sx(ipp,jl)+sx(inn,jl)+sx(i,j3NNNN)) !J3(NNNN)
     &  -sy(i,j)*X2J3*(sy(ipp,jl)+sy(inn,jl)+sy(i,j3NNNN)) !J3(NNNN)  
     &  -sz(i,j)*X3J3*(sz(ipp,jl)+sz(inn,jl)+sz(i,j3NNNN)) !J3(NNNN)

     &  -sx(i,j)*KX*sx(ip,j)  !KX-X-bond
     &  -sy(i,j)*KY*sy(in,j)  !KY-Y-bond
     &  -sz(i,j)*KZ*sz(i,jl)  !KZ-Z-bond  

     &  -GamaX*(sy(i,j)*sz(ip,j)+sy(ip,j)*sz(i,j)) !GammaX
     &  -GamaY*(sx(i,j)*sz(in,j)+sx(in,j)*sz(i,j)) !GammaY
     &  -GamaZ*(sx(i,j)*sy(i,jl)+sx(i,jl)*sy(i,j)) !GammaZ

C***************************************************C
C**Creating a new orientation for the selected spin C       
C***************************************************C
           CALL RANDOM_NUMBER(r1)
           CALL RANDOM_NUMBER(r2)
 
           th=acos(2*r1-1)
           phi= 2*pi*r2
 
           flipx(i,j)=S*sin(th)*cos(phi)
           flipy(i,j)=S*sin(th)*sin(phi)
           flipz(i,j)=S*cos(th)

C**** Calculate the spin energy after flipping *********C
      E_flip=-flipx(i,j)*X1J1*(sx(ip,j)+sx(in,j)+sx(i,jl))  !J1(NN)
     &-flipy(i,j)*X2J1*(sy(ip,j)+sy(in,j)+sy(i,jl))         !J1(NN)
     &-flipz(i,j)*X3J1*(sz(ip,j)+sz(in,j)+sz(i,jl))         !J1(NN)

     &-flipx(i,j)*X1J2*(sx(ip,jp)+sx(ip,j) +sx(in,jn)+sx(in,jp))!J2(NNN)
     &-flipy(i,j)*X2J2*(sy(ip,jp)+sy(ip,jn)+sy(in,jn)+sy(in,jp))!J2(NNN)
     &-flipz(i,j)*X3J2*(sz(ip,jp)+sz(ip,jn)+sz(in,jn)+sz(in,jp))!J2(NNN) 

     &-flipx(i,j)*X1J2*(sx(ipp,j)+sx(inn,j))                !J2(NNN)
     &-flipy(i,j)*X2J2*(sy(ipp,j)+sy(inn,j))                !J2(NNN)
     &-flipz(i,j)*X3J2*(sz(ipp,j)+sz(inn,j))                !J2(NNN)

     &-flipx(i,j)*X1J3*(sx(ipp,jl)+sx(inn,jl)+sx(i,j3NNNN)) !J3(NNNN)
     &-flipy(i,j)*X2J3*(sy(ipp,jl)+sy(inn,jl)+sy(i,j3NNNN)) !J3(NNNN)
     &-flipz(i,j)*X3J3*(sz(ipp,jl)+sz(inn,jl)+sz(i,j3NNNN)) !J3(NNNN)

     &-flipx(i,j)*KX*sx(ip,j)                               !KX-X-BOND
     &-flipy(i,j)*KY*sy(in,j)                               !KY-Y-BOND
     &-flipz(i,j)*KZ*sz(i,jl)                               !KZ-Z-BOND

     &-GamaX*(flipy(i,j)*sz(ip,j)+sy(ip,j)*flipz(i,j))      !GammaX
     &-GamaY*(flipx(i,j)*sz(in,j)+sx(in,j)*flipz(i,j))      !GammaY
     &-GamaZ*(flipx(i,j)*sy(i,jl)+sx(i,jl)*flipy(i,j))      !GammaZ
          
C*********************************C
C        Flipping Probablity      C
C*********************************C
          dE=(E_flip-E)
          prob=DEXP(-dE/t) 
          CALL RANDOM_NUMBER(r)
          IF(r.le.prob) then
          sx(i,j)= flipx(i,j)
          sy(i,j)= flipy(i,j)
          sz(i,j)= flipz(i,j)
          ELSE
          sx(i,j)=sx(i,j)
          sy(i,j)=sy(i,j)
          sz(i,j)=sz(i,j)
          ENDIF
          ENDDO
C**********************************************C
C******* Averaging after equillibrium *********C
C**********************************************C          
       if(ic.ge.itequi)then

       xm1=0.d0
       xm2=0.d0
       xm3=0.d0
       en=0.d0

       DO i=1,L
       DO j=1,L

       call Periodic_Boundary_Condition(i,j,ip,in,jp,jn,ipp,inn,jl,
     & j3NNNN,L)

       xm1=xm1+sx(i,j)
       xm2=xm2+sy(i,j)
       xm3=xm3+sz(i,j)

       en=en-0.5d0*(sx(i,j)*X1J1*(sx(ip,j)+sx(in,j)+sx(i,jl))   !J1(NN)
     &  +sy(i,j)*X2J1*(sy(ip,j)+sy(in,j)+sy(i,jl))              !J1(NN)
     &  +sz(i,j)*X3J1*(sz(ip,j)+sz(in,j)+sz(i,jl))              !J1(NN)

     &  +sx(i,j)*X1J2*(sx(ip,jp)+sx(ip,jn)+sx(in,jn)+sx(in,jp)) !J2(NNN)
     &  +sy(i,j)*X2J2*(sy(ip,jp)+sy(ip,jn)+sy(in,jn)+sy(in,jp)) !J2(NNN)
     &  +sz(i,j)*X3J2*(sz(ip,jp)+sz(ip,jn)+sz(in,jn)+sz(in,jp)) !J2(NNN)

     &  +sx(i,j)*X1J2*(sx(ipp,j)+sx(inn,j))                     !J2(NNN)
     &  +sy(i,j)*X2J2*(sy(ipp,j)+sy(inn,j))                     !J2(NNN)
     &  +sz(i,j)*X3J2*(sz(ipp,j)+sz(inn,j))                     !J2(NNN)
     
     &  +sx(i,j)*X1J3*(sx(ipp,jl)+sx(inn,jl)+sx(i,j3NNNN))      !J3(NNNN)
     &  +sy(i,j)*X2J3*(sy(ipp,jl)+sy(inn,jl)+sy(i,j3NNNN))      !J3(NNNN)
     &  +sz(i,j)*X3J3*(sz(ipp,jl)+sz(inn,jl)+sz(i,j3NNNN))      !J3(NNNN)

     &  +KX*sx(i,j)*sx(ip,j)                                    !KX-X bond
     &  +KY*sy(i,j)*sy(in,j)                                    !KY-Y bond
     &  +KZ*sz(i,j)*sz(i,jl)                                    !KZ-Z bond  

     &  +GamaX*(flipy(i,j)*sz(ip,j)+sy(ip,j)*flipz(i,j))        !GammaX
     &  +GamaY*(flipx(i,j)*sz(in,j)+sx(in,j)*flipz(i,j))        !GammaY
     &  +GamaZ*(flipx(i,j)*sy(i,jl)+sx(i,jl)*flipy(i,j)))       !GammaZ

      ENDDO
      ENDDO

      xm1=xm1/LSQ
      xm2=xm2/LSQ
      xm3=xm3/LSQ
      en=en/LSQ
      xm1=dabs(xm1)
      xm2=dabs(xm2)
      xm3=dabs(xm3)
      xm=sqrt(xm1**2+ xm2**2+xm3**2)
      am1=am1+xm
      am2=am2+xm**2
      am3=am3+xm**4
      energy1=energy1+en
      energy2=energy2+en*en

      ENDIF
      ENDDO

      am1=am1/itav                 ! <m>
      am2=am2/itav                 ! <m**2>
      am3=am3/itav                 ! <m**4>
      energy1=energy1/itav         ! <E>
      energy2=energy2/itav         ! <E**2>
      UL=1-(am3)/(3.d0*am2**2)     ! Binder-Cumulant
      flucE=LSQ*(energy2-energy1**2)/(t*t)      !Cv
      flucM=LSQ*(am2-am1**2)/t                  !Chi

C***** Writing the Data File ****************C

      OPEN(4,file='HSO',status='old',POSITION='APPEND')
        WRITE(4,*)t,am1,flucM,energy1,flucE,UL
      CLOSE(4)     
C**************************************************************C
C******* This file is for spin dynamics during transition *****C
C**************************************************************C

      CALL Spin_Configuration_During_Transition(i,j,L,
     &sx,sy,sz,t)

      WRITE(*,*)t,am1,energy1,flucE

      IF(temp.eq.tmin)then 
        WRITE(*,32)
      ENDIF

      ENDDO

      CALL SYSTEM('date')
      END  
      
C****** Main Program End Here **********C

C*****CONTAINS SUBROUTINE **************C

      SUBROUTINE Periodic_Boundary_Condition(i,j,ip,in,jp,jn,ipp,inn,jl,
     &j3NNNN,L)
      IMPLICIT NONE
      INTEGER i,j,ip,in,jp,jn,L
      INTEGER ipp,inn,jl,j3NNNN

       ip=i+1
       in=i-1
       jp=j+1
       jn=j-1
       ipp=i+2
       inn=i-2
 
        if(i.eq.L)ip=1
        if(i.eq.1)in=L
        if(j.eq.L)jp=1
        if(j.eq.1)jn=L

        if(i.eq.L-1)ipp=1
        if(i.eq.2)inn=L
        if(i.eq.L)ipp=2
        if(i.eq.1)inn=L-1
 
        IF (mod(i+j,2)==0) THEN
          jl=j+1
          j3NNNN=j-1
          if(j.eq.L)jl=1
          if(j.eq.1)j3NNNN=L
        ELSE
          jl=j-1
          j3NNNN=j+1
          if(j.eq.1)jl=L
          if(j.eq.L)j3NNNN=1
        ENDIF
 
       RETURN
       END
      
      SUBROUTINE Spin_Configuration_During_Transition(i,j,L,
     &sx,sy,sz,t)
      IMPLICIT NONE          
      INTEGER i,j,L
      REAL*8 sx(L,L),sy(L,L),sz(L,L),t

      CALL SYSTEM ('rm spin_dynamics')
      OPEN(9,file='spin_dynamics',status='new')
      WRITE(9,*) '## t ## i ## j ## S_x ## S_y ## S_z ##'

      DO i =1,L
        DO j =1,L 
          WRITE(9,*) t,i,j,sx(i,j),sy(i,j),sz(i,j)
        ENDDO
      ENDDO
      CLOSE(9)

      RETURN
      END



