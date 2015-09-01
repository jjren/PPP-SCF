    ! THIS PROGRAM GENERATES THE POSITION CORDINATES OF ALL 60
! CARBON ATOMS OF THE FULLERENE STRUCTURE EXPRESSED AS
! X(NP,NV,I) WHERE NP=PENTAGON INDEX,NV IDENTIFIES THE VERTEX  
! WITHIN THE PENTAGON &I IS THE COORDINATE (=1forX, etc.)
! INPUT ARE TWO BOND LENGTHS (LINES 4,5 AFTER COMMENTS) :
! S=LENGTH OF SINGLE BOND; D=LENGTH OF DOUBLE BOND.
! UNITS:OUT UNIT WILL BE SAME AS INPUT UNIT
! THIS IS NOT THE "STANDARD" ORIENTATION, BUT ONE WITH PENTAGONS
! AT TOP & BOTTOM	(IN "ST." D.BONDS INTERSECT COORDINATE AXES)
! Written by Prof. Keya Dharamvir
!
    subroutine c60_gen(S,DB,rr_atom,rc)
    
    
 !   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 !   DIMENSION X(12,5,3),EM1(3,3),EM2(3,3)
 !   dimension rr_atom(3,60),rc(3)
 
 IMPLICIT NONE
 
    integer::ok,IJ,I,J,K,NP,NMP,NV,icount,KK
    real(kind=8),intent(in)::S,DB,rc(3)
    real(kind=8)::BLIMIT,RATIO,STANDRD,PI,RHO,F,A,ALPHA,SS,Z0,G,GAMMA,E,&
                  ETA,Q,SQ,BUCKYR,XX,ZE,D,D1,C,C1,RADF,DIS,XNOW,X(12,5,3),&
		  EM1(3,3),EM2(3,3),M,rr_atom(3,60)
!    real(kind=8),allocatable::X,EM1,EM2
!    allocate(X(12,5,3),EM1(3,3),EM2(3,3),stat=ok)
!      if(ok/=0)then
!        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING X,EM1 or EM2'
!	stop
!      end if   
    
    OPEN(UNIT=3,FILE='c60cordandists.DAT',STATUS='UNKNOWN')
!	S=1.45D0
!	D=1.40D0
!	S=1.45D0
! DB=1.270D0
! print*,'Give the Single and Double Bond Lengths in that order'
! read*,s,db
        BLIMIT=2.1d0/0.6d0
! this will give coords of only those atoms which r linked to 1st
! ,when max linear compression=0.6
!	DO 56 NC=1,31
    !	DB= DB + 0.005
    RATIO=S/DB
    STANDRD=3.5485013897d0
        PI=4.d0*atan(1.d0)
! WRITE(3,220)S,DB,RATIO
        RHO=2.d0+ DB/S
    F=PI/180.d0
    A=54.d0*F
    ALPHA=DCOS(A)
    SS=S/(2.d0*ALPHA)  
    Z0=(SS/2.d0)*(1.d0+(RHO-1.d0)*DSQRT(5.0D0))
    G=72.0D0*F
    GAMMA=DCOS(G)
    E=144.d0*F
    ETA=DCOS(E)
    X(1,1,1)=0.d0
    X(1,1,2)=SS
    X(1,2,1)=SS*DSQRT(1.d0-GAMMA*GAMMA)
    X(1,2,2)=SS*GAMMA
    X(1,3,1)=SS*DSQRT(1.d0-ETA*ETA)   
    X(1,3,2)=SS*ETA
    DO  M=1,5
        X(1,M,3)=Z0
    END DO
        X(1,4,1)= -X(1,3,1)
        X(1,4,2)=  X(1,3,2)
        X(1,5,1)= -X(1,2,1)
        X(1,5,2)=  X(1,2,2)
                SQ=0.
        DO IJ=1,3
            Q =X(1,1,IJ)
            SQ=SQ+ Q*Q
        END DO
            BUCKYR=SQRT(SQ)
        !	WRITE(3,110)BUCKYR
        !	WRITE(3,11)
            11 FORMAT(/)
                        C=1.0-(8.d0/5.d0)*ALPHA*ALPHA
            C1=DSQRT(1.-C*C)
            EM1(1,1)=-1.d0
            EM1(1,2)=0.d0
            EM1(1,2)=0.d0
            EM1(1,3)=0.d0
            EM1(2,1)=0.d0	
            EM1(2,2)=-C
            EM1(2,3)=C1
            EM1(3,1)=0.
            EM1(3,2)=C1
            EM1(3,3)=C
            DO I=1,5
                DO J=1,3
                    XX=0.
                    DO K=1,3
                        XX=EM1(J,K)*X(1,I,K)+XX
                    END DO
                        X(2,I,J)=XX
                END DO
            END DO
                        ZE=72.d0*F
                        D=DCOS(ZE)
                        D1=DSIN(ZE)
                        EM2(3,3)=1.d0
                        EM2(1,3)=0.d0
                        EM2(3,1)=0.d0
                        EM2(2,3)=0.d0
                        EM2(3,2)=0.d0
                        EM2(1,1)=D
                        EM2(2,2)=D
                        EM2(1,2)=D1
                        EM2(2,1)=-D1
                        DO NP = 3,6
                            NMP=NP-1
                            DO NV=1,5
                                DO J=1,3
                                    XX=0.
                                    DO K=1,3
                                        XX=EM2(J,K)*X(NMP,NV,K)+XX
                                    END DO        
                                        X(NP,NV,J)=XX
                                 END DO
		            END DO  
			 END DO      
                                    ! Generating the remaining coordinates
                                        DO NP=7,12
                                            NMP=NP-6
                                            DO NV=1,5
                                                DO J=1,3
                                                    X(NP,NV,J)= -X(NMP,NV,J)
                                                END DO
				             END DO
					 END DO
                                                    RADF=STANDRD/BUCKYR
                                                !	x(1,1,1)=x(1,1,1)*RADF
                                                !	x(1,1,2)=x(1,1,2)*RADF
                                                !	x(1,1,3)=x(1,1,3)*RADF
                                                ! only upto I=6 (leaving 4,5) becos rest have v large distances
                                                    icount=0
                                                    DO I=1,12
                                                    !	if (I.eq.4)go to 100
                                                    !	if (I.eq.5)go to 100
                                                        do j=1,5
                                                            icount=icount+1
                                                            DIS=0.d0
                                                            DO KK=1,3
                                                                XNOW=X(I,J,KK)*RADF
                                                                DIS=DIS+(XNOW-X(1,1,KK))**2
                                                            END DO
                                                            DIS=SQRT(DIS)
                                                            ! IF(DIS.GT.BLIMIT)GO TO 100
                                                                do k=1,3
                                                                    rr_atom(k,icount)=x(i,j,k)+rc(k)
                                                                end do
                                                            !
                                                            ! WRITE(3,450)I,J,(X(I,J,K),K=1,3),DIS
                                                            ! write(*,33)(x(i,j,k),k=1,3)
                                                          END DO
                                                      END DO
                                                            ! 6	CONTINUE
                                                            !    450	FORMAT (2I3,4F10.6)
                                                            !    220	FORMAT('BONDLENGTHS:  S=',F11.7,3X,'D=',F11.7,'Ratio=',F6.4) 
                                                            !    110	FORMAT('RADIUS OF THAT BUCKYBALL=',F14.10)
                                                            !    33 format(f11.7, 3x, f11.7, 3x, f11.7)
                                                                return
                                                                END
