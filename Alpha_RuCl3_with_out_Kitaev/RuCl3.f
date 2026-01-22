
C*******************************************************************C
C*******************************************************************C
C*    Classical Monte-Carlo Simulation In Honeycomb  Lattice       *C
C*    In Presence of NN(J1) ,NNN(J2), NNNN(J3) Interaction         *C
C*                   Alpha-RuCl3                                   *C
C*    Here T  is in meV unit.Here 1meV=11.6045 K                   *C
C*                                                                 *C      
C*    -------------------------------------------------------      *C
C*    Deleloped by: Nepal Banerjee ,University of Seoul,Seoul      *C
C*    Contact:-nepalbanerjee36@gmail.com/nb.uos1989@gail.com       *C
C*    Mob:+91-9474172495/+91-9816075495/+91-8250773046             *C
C*    -------------------------------------------------------      *C
C*                                                                 *C
C*******************************************************************C
C*******************************************************************C
         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
         PARAMETER(L=4,LSQ=L*L,itmx=900000,
     &itset=500000,idiv=itmx-itset,pi=4*atan(1.d0))
         DIMENSION hs1(L,L),hs2(L,L),hs3(L,L)
         DIMENSION flip1(L,L),flip2(L,L),flip3(L,L)  

         S=0.5d0   !SPIN MAG
         
         CALL SYSTEM('rm HSO')
   
         OPEN(4,file='HSO',status='new')
            WRITE(4,*)'## T ## M ## Chi ## E ## Cv ##'
         CLOSE(4)

C****************************************************C
C***Initial Spin Configuration                 ******C
C****************************************************C
        DO i=1,L
           DO j=1,L
              CALL RANDOM_NUMBER(r1)
              CALL RANDOM_NUMBER(r2)
              th=acos(2*r1-1)   
              phi=2*pi*r2
              hs1(i,j)=s*sin(th)*cos(phi) 
              hs2(i,j)=s*sin(th)*sin(phi)
              hs3(i,j)=s*cos(th)   
           ENDDO
        ENDDO

C*** First-Nearest Neighbour  *****C
          X1J1=10.69d0     !J1(FM)
          X2J1=10.69d0     !J1(FM)
          X3J1=10.69d0     !J1(FM)

C****Second-Nearest Neighbour****C 
          X1J2=-1.69d0     !J2(AFM)
          X2J2=-1.69d0     !J2(AFM)
          X3J2=-1.69d0     !J2(AFM)

C****Third-Nearest Neighbour ****C
          X1J3=2.54d0      !J3(FM)
          X2J3=2.54d0      !J3(FM)
          X3J3=2.54d0      !J3(FM)

C****Temperature loop start here ****C

         it1=500  ;             !Maximum temperature
         it2=10   ;             !Minimum temperature
         idt=10   ;             !Temperature steps 


         PRINT*, 'Simulation is running ....'
         
         DO it=it1,it2,-idt     !Temperature Loop 

            t=it/100.0
            am1=0.d0
            am2=0.d0
            energy1=0.d0
            energy2=0.d0


            DO ic=1,itmx        !Iteration for Equillibrium and Average

              DO id=1,LSQ
                
                 CALL RANDOM_NUMBER(r1)
                 CALL RANDOM_NUMBER(r2)

                 i=1+r1*L
                 j=1+r2*L

                 CALL Periodic_Boundary_Condition(i,j,ip,in,jp,jn,jl,
     &j3NNNN,ipp,inn,L)

C****** First Nearest neighbour ************************C
       E=-hs1(i,j)*X1J1*(hs1(ip,j)+hs1(in,j)+hs1(i,jl))              !J1
     &  -hs2(i,j)*X2J1*(hs2(ip,j)+hs2(in,j)+hs2(i,jl))               !J1
     &  -hs3(i,j)*X3J1*(hs3(ip,j)+hs3(in,j)+hs3(i,jl))               !J1

C******* Second Nearest Neihbour ***********************C
     &  -hs1(i,j)*X1J2*(hs1(ip,jp)+hs1(ip,jn)+hs1(in,jn)+hs1(in,jp)) !J2
     &  -hs2(i,j)*X2J2*(hs2(ip,jp)+hs2(ip,jn)+hs2(in,jn)+hs2(in,jp)) !J2
     &  -hs3(i,j)*X3J2*(hs3(ip,jp)+hs3(ip,jn)+hs3(in,jn)+hs3(in,jp)) !J2
     &  -hs1(i,j)*X1J2*(hs1(ipp,j)+hs1(inn,j))                       !J2
     &  -hs2(i,j)*X2J2*(hs2(ipp,j)+hs2(inn,j))                       !J2
     &  -hs3(i,j)*X3J2*(hs3(ipp,j)+hs3(inn,j))                       !J2

C******** Third Nearest Neighbour **********************C
     &  -hs1(i,j)*X1J3*(hs1(ipp,jl)+hs1(inn,jl)+hs1(i,j3NNNN))       !J3
     &  -hs2(i,j)*X2J3*(hs2(ipp,jl)+hs2(inn,jl)+hs2(i,j3NNNN))       !J3
     &  -hs3(i,j)*X3J3*(hs3(ipp,jl)+hs3(inn,jl)+hs3(i,j3NNNN))       !J3

C***************************************************C
C**New Preferable Orientation of Heisenberg spin ***C
C***************************************************C
        CALL RANDOM_NUMBER(r1)
        CALL RANDOM_NUMBER(r2)

        th=acos(2*r1-1)
        phi=2*pi*r2
        
        flip1(i,j)=s*sin(th)*cos(phi)
        flip2(i,j)=s*sin(th)*sin(phi)
        flip3(i,j)=s*cos(th)

C******First Nearest neighbour (NN)************************C 
      E1=-flip1(i,j)*X1J1*(hs1(ip,j)+hs1(in,j)+hs1(i,jl))            !J1
     &-flip2(i,j)*X2J1*(hs2(ip,j)+hs2(in,j)+hs2(i,jl))               !J1
     &-flip3(i,j)*X3J1*(hs3(ip,j)+hs3(in,j)+hs3(i,jl))               !J1

C******Second Nearest Neighbour (NNN)**********************C
     &-flip1(i,j)*X1J2*(hs1(ip,jp)+hs1(ip,jn)+hs1(in,jn)+hs1(in,jp)) !J2
     &-flip2(i,j)*X2J2*(hs2(ip,jp)+hs2(ip,jn)+hs2(in,jn)+hs2(in,jp)) !J2
     &-flip3(i,j)*X3J2*(hs3(ip,jp)+hs3(ip,jn)+hs3(in,jn)+hs3(in,jp)) !J2
     &-flip1(i,j)*X1J2*(hs1(ipp,j)+hs1(inn,j))                       !J2
     &-flip2(i,j)*X2J2*(hs2(ipp,j)+hs2(inn,j))                       !J2
     &-flip3(i,j)*X3J2*(hs3(ipp,j)+hs3(inn,j))                       !J2

C******Third Nearest Neighbour (NNNN)***********************C
     &-flip1(i,j)*X1J3*(hs1(ipp,jl)+hs1(inn,jl)+hs1(i,j3NNNN))       !J3
     &-flip2(i,j)*X2J3*(hs2(ipp,jl)+hs2(inn,jl)+hs2(i,j3NNNN))       !J3
     &-flip3(i,j)*X3J3*(hs3(ipp,jl)+hs3(inn,jl)+hs3(i,j3NNNN))       !J3

C******Calculating Flipping Probablity******C
        dE=(E1-E)
        prob=dexp(-dE/t) 

        if(rand().le.prob) then
        hs1(i,j)= flip1(i,j)
        hs2(i,j)= flip2(i,j)
        hs3(i,j)= flip3(i,j)
        else
        hs1(i,j)=hs1(i,j)
        hs2(i,j)=hs2(i,j)
        hs3(i,j)=hs3(i,j)
        endif

        ENDDO

C********* Averaging After Equillibrium ********C

       if(ic.ge.itset)then      

       xm1=0.d0
       xm2=0.d0
       xm3=0.d0
       en=0.d0

       DO i=1,L
       DO j=1,L
       

       CALL Periodic_Boundary_Condition(i,j,ip,in,jp,jn,jl,j3NNNN,
     & ipp,inn,L)
C***********************************************
       xm1=xm1+hs1(i,j)
       xm2=xm2+hs2(i,j)
       xm3=xm3+hs3(i,j)       

C****** First Nearest Neighbour ************************C
       en=en-0.5d0*(hs1(i,j)*X1J1*(hs1(ip,j)+hs1(in,j)+hs1(i,jl))    !J1
     &  +hs2(i,j)*X2J1*(hs2(ip,j)+hs2(in,j)+hs2(i,jl))               !J1
     &  +hs3(i,j)*X3J1*(hs3(ip,j)+hs3(in,j)+hs3(i,jl))               !J1

C******* Second Nearest Neighbour ***********************C
     &  +hs1(i,j)*X1J2*(hs1(ip,jp)+hs1(ip,jn)+hs1(in,jn)+hs1(in,jp)) !J2
     &  +hs2(i,j)*X2J2*(hs2(ip,jp)+hs2(ip,jn)+hs2(in,jn)+hs2(in,jp)) !J2
     &  +hs3(i,j)*X3J2*(hs3(ip,jp)+hs3(ip,jn)+hs3(in,jn)+hs3(in,jp)) !J2
     &  +hs1(i,j)*X1J2*(hs1(ipp,j)+hs1(inn,j))                       !J2
     &  +hs2(i,j)*X2J2*(hs2(ipp,j)+hs2(inn,j))                       !J2
     &  +hs3(i,j)*X3J2*(hs3(ipp,j)+hs3(inn,j))                       !J2

C******** Third Nearest Neighbour *****************C
     &  +hs1(i,j)*X1J3*(hs1(ipp,jl)+hs1(inn,jl)+hs1(i,j3NNNN))       !J3
     &  +hs2(i,j)*X2J3*(hs2(ipp,jl)+hs2(inn,jl)+hs2(i,j3NNNN))       !J3
     &  +hs3(i,j)*X3J3*(hs3(ipp,jl)+hs3(inn,jl)+hs3(i,j3NNNN)))      !J3
 
      ENDDO
      ENDDO

C*************************************************C
      xm1=xm1/LSQ
      xm2=xm2/LSQ
      xm3=xm3/LSQ
      en=en/LSQ
      xm1=dabs(xm1)
      xm2=dabs(xm2)
      xm3=dabs(xm3)
      xm=sqrt(xm1**2+ xm2**2+xm3**2)
      am1=am1+xm
      am2=am2+xm*xm
      energy1=energy1+en
      energy2=energy2+en*en
      endif
      ENDDO

C**************************************************C      
      am1=am1/idiv
      am2=am2/idiv
      energy1=energy1/idiv
      energy2=energy2/idiv
      flucE=LSQ*(energy2-energy1*energy1)/(t*t)
      flucM=LSQ*(am2-(am1)**2)/t

C**************************************************C
      OPEN(4,file='HSO',status='old')
      CALL fseek(4,0,2)
      WRITE(4,*)t,am1,flucM,energy1,flucE
      CLOSE(4)
      ENDDO       
      PRINT*,'Congratulation!'
      PRINT*,'Simulation is completed sucessfully' 
      END  

C**********MAIN PROGRAMM END HERE ******************C

      SUBROUTINE Periodic_Boundary_Condition(i,j,ip,in,jp,jn,jl,j3NNNN,
     & ipp,inn,L)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

       ip=i+1
       in=i-1
       jp=j+1
       jn=j-1
       ipp=i+2
       inn=i-2

       if(i.eq.L)ip=1
       if(i.eq.1)in=L
       if(i.eq.L-1)ipp=1
       if(i.eq.L)ipp=2
       if(i.eq.2)inn=L
       if(i.eq.1)inn=L-1
       if(j.eq.L)jp=1
       if(j.eq.1)jn=L
       
       if (mod(j,2).gt.0)then 
          if (mod(i,2).gt.0) then
            jl=j+1
            j3NNNN=j-1
            if(j.eq.L)jl=1
            if(j.eq.1)j3NNNN=L
          else
            jl=j-1
            j3NNNN=j+1
            if(j.eq.1)jl=L
            if(j.eq.L)j3NNNN=1
          endif
       else
          if (mod(i,2).gt.0)then
            jl=j-1
            j3NNNN=j+1
            if(j.eq.1)jl=L
            if(j.eq.L)j3NNNN=1
          else
            jl=j+1
            j3NNNN=j-1
            if(j.eq.L)jl=1
            if(j.eq.1)j3NNNN=L
          endif
       endif 

       RETURN
       END
