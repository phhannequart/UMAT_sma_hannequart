!		UMAT Subroutine for ABAQUS
!		For modelling a polycrystalline shape memory alloy wire
!		Authors : Philippe Hannequart (1,2), Michael Peigney(1), Jean-François Caron(1), Emmanuel Viglino(2)
!		(1) Universite Paris-Est, Laboratoire Navier (UMR 8205), CNRS, Ecole des Ponts ParisTech, IFSTTAR, 77455 Marne la vallee, France
!		(2) Arcora, Groupe Ingerop, Rueil-Malmaison, France
!
!		It refers to the model published in paper 
!
!		The number of state variables NSTATV in Abaqus must be twice the number of crystal orientations
!		NSTATV = 2 is a SMA monocrystal
!
!		The material properties in Abaqus must be set as follows :
!		PROPS(1) : Young's Modulus
!		PROPS(2) : Latent heat parameter
!       PROPS(3) : Reference temperature
!		PROPS(4) : Dissipation parameter for self-accomodated martensite 
!		PROPS(5) : Dissipation parameter for oriented martensite
!		PROPS(6) to PROPS(5 + NSTATV/2) are the values of the maximum phase transformation strains
!
!
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      IMPLICIT NONE
      CHARACTER*80 CMNAME
      INTEGER I,J,NPROPS
!
      DOUBLE PRECISION STRESS, STRAN, DSTRAN, TEMP, DDSDDE, PREDEF
      DOUBLE PRECISION TIME(2), COORDS(3),DDSDDT(1),DFGRD0(3,3)
      DOUBLE PRECISION DFGRD1(3,3), DPRED(1),DROT(3,3),DRPLDE
      INTEGER NTENS
      INTEGER NSTATV ! Number of state variables, must be a multiple of 2 !
      INTEGER CELENT,DRPLDT,DTEMP,DTIME,KSPT,KSTEP,LAYER,NDI
      INTEGER NOEL,NPT,NSHR,PNEWDT,RPL,SCD,SPD,SSE,KINC
      DOUBLE PRECISION PROPS(NPROPS), STATEV(NSTATV), DNORI
      DOUBLE PRECISION EPSTR(NSTATV), Q(7*NSTATV/2),S(7*NSTATV/2)
      DOUBLE PRECISION X(7*NSTATV/2),MULTI(NSTATV),CR(NSTATV/2)
      DOUBLE PRECISION MATRIX(7*NSTATV/2,7*NSTATV/2)
      DOUBLE PRECISION EV, ETA, ZERO, ONE
      INTEGER NORI ! Number of crystalline orientations
      INTEGER NMAT ! Size of matrix M in the LCP
!
      ZERO=0.D0;  ONE=1.D0; 
      NORI=NSTATV/2
      NMAT=7*NORI
      DNORI=DBLE(NORI)
      EV = STRAN + DSTRAN   ! Update strain
!
!     Texture Functions
!
      DO I=1,NORI
         EPSTR(I) = ZERO     ! Self-accomodated martensite
         EPSTR(NORI+I) = PROPS(5+I) ! Maximum phase transformation strain (Texture function)
         CR(I) = 1/DNORI ! Volume fraction of each crystalline orientation, can be adapted to the chosen model
      END DO	 
!
!     Initialize M
!
      DO I=1,NMAT
          DO J=1,NMAT
              MATRIX(I,J)=ZERO
          END DO
      END DO
!
!     Fill M - Upper left corner
!
      DO I=2*NORI+1,3*NORI
          DO J=2*NORI+1,3*NORI
              MATRIX(I,J)= EPSTR(I-NORI)/EPSTR(J-NORI)
              MATRIX(I,J+NORI)= -EPSTR(I-NORI)/EPSTR(J-NORI)
              MATRIX(I+NORI,J)= -EPSTR(I-NORI)/EPSTR(J-NORI)
              MATRIX(I+NORI,J+NORI)= EPSTR(I-NORI)/EPSTR(J-NORI)
          END DO
      END DO
!	  
!     Fill M - Upper right corner
!
      DO I=1,NORI
          MATRIX(I,I+4*NORI)=-ONE
          MATRIX(I,I+6*NORI)=+ONE
          MATRIX(I+NORI,I+4*NORI)=+ONE
          MATRIX(I+NORI,I+6*NORI)=-ONE
          MATRIX(I+2*NORI,I+5*NORI)=-ONE
          MATRIX(I+2*NORI,I+6*NORI)=+ONE
          MATRIX(I+3*NORI,I+5*NORI)=+ONE
          MATRIX(I+3*NORI,I+6*NORI)=-ONE
      END DO
!
!     Fill M - Lower left corner  
!
      DO I=1,NORI
          MATRIX(I+4*NORI,I)=+ONE
          MATRIX(I+4*NORI,I+NORI)=-ONE
          MATRIX(I+5*NORI,I+2*NORI)=+ONE
          MATRIX(I+5*NORI,I+3*NORI)=-ONE
          MATRIX(I+6*NORI,I)=-ONE
          MATRIX(I+6*NORI,I+NORI)=+ONE
          MATRIX(I+6*NORI,I+2*NORI)=-ONE
          MATRIX(I+6*NORI,I+3*NORI)=+ONE
      END DO
!
!     Fill Q      
!
      DO I=1,NSTATV
          MULTI(I)=EPSTR(I)*STATEV(I)
      END DO
      DO I=1,NORI
          Q(I)=PROPS(4)+PROPS(2)*((TEMP-PROPS(3))/PROPS(3))
          Q(I+NORI)=PROPS(4)-PROPS(2)*((TEMP-PROPS(3))/PROPS(3))
          Q(I+2*NORI)=PROPS(5)
     1 +PROPS(2)*((TEMP-PROPS(3))/PROPS(3)) 
     2 -PROPS(1)*EPSTR(I+NORI)*(EV-SUM(MULTI))
          Q(I+3*NORI)=PROPS(5)
     1 -PROPS(2)*((TEMP-PROPS(3))/PROPS(3)) 
     2 +PROPS(1)*EPSTR(I+NORI)*(EV-SUM(MULTI)) 
         Q(I+4*NORI)=STATEV(I)*PROPS(1)*(EPSTR(I+NORI))**2
         Q(I+5*NORI)=STATEV(I+NORI)*PROPS(1)*(EPSTR(I+NORI))**2
         Q(I+6*NORI)=(CR(I)-STATEV(I)-STATEV(I+NORI))
     1  *PROPS(1)*(EPSTR(I+NORI))**2
      END DO
!
!     Compute the updated volume fractions of martensite for each crystal orientation
!     through the solving of a Linear Complementarity Problem (LCP)
!
      CALL LEMKE (MATRIX,Q,NMAT,X)
      S = MATMUL(MATRIX,X)   ! S = M X + Q
      DO I=1,NMAT
         S(I) = S(I) + Q(I)
      END DO
!
      ETA = 1.0D-6   ! Tolerance parameter
      DO I = 1,NORI
          STATEV(I) = S(4*NORI+I)/(PROPS(1)*(EPSTR(NORI+I))**2)   
          STATEV(NORI+I) = S(5*NORI+I)/(PROPS(1)*(EPSTR(NORI+I))**2) 
		  ! Update state variables
      END DO
!
      DO I = 1,NSTATV
         MULTI(I) =  STATEV(I)*EPSTR(I)      
      END DO
      STRESS = PROPS(1)*(EV-SUM(MULTI))   ! Compute the new stress
!
      DDSDDE = PROPS(1) ! Default value of the tangent modulus
!	  
      DO I=1,NORI
         IF ((ABS(X(2*NORI+I)-X(3*NORI+I)) .GT. ETA) .AND. 
     1 (STATEV(NORI+I) .GT. ETA) .AND.
     2 (STATEV(NORI+I) .LT. CR(I)-ETA))	 THEN
             DDSDDE = ZERO ! Zero tangent modulus
         END IF
      END DO	  
!	  
      RETURN
      END SUBROUTINE
!
!     --------- Solving the LCP by Lemke's method -----
!
!     Adapted to Fortran from the Matlab code num4lcp 
!     Published online at code.google.com/p/num4lcp/
!     Kenny Erleben (2012)
!     The algorithm uses a modified Lemke's algorithm (complementary pivoting)
!     with a covering ray of ones.  

      SUBROUTINE LEMKE (MATRXX,QU,NMATT,XI)
!
      IMPLICIT NONE
!
      INTEGER I,J,ITER,MAXITR,T,ENTER,LEAVE,LVIND,NMATT,INFO
      DOUBLE PRECISION ZERTOL, PIVTOL, TVAL, RATIO,THETA
      DOUBLE PRECISION ZERO, ONE
      DOUBLE PRECISION B(NMATT,NMATT),BSTOR(NMATT,NMATT)
      DOUBLE PRECISION MATRXX(NMATT,NMATT)
      DOUBLE PRECISION XSUB(NMATT),XI(NMATT),MX(NMATT),BE(NMATT)
      DOUBLE PRECISION D(NMATT),THETAV(NMATT)
      DOUBLE PRECISION ONES(NMATT),BONES(NMATT),QU(NMATT)
      DOUBLE PRECISION Z(2*NMATT)
      INTEGER BAS(NMATT),IPVT(NMATT)
      LOGICAL INDICE(NMATT)
!
      ZERO=0.D0; ONE=1.D0;  
      ZERTOL = 1.0D-5; PIVTOL = 1.0D-8; MAXITR = 10000;
!	  
! 	Check if trivial solution exists	  	  
!
      IF (MINVAL(QU) .GE. ZERO) THEN
          DO I=1,NMATT
              XI(I) = ZERO
          END DO
          RETURN
      END IF	  
      DO I=1,2*NMATT
         Z(I) = ZERO
      END DO
!
! 	Determine initial basis	
!	  
      DO I=1,NMATT
          INDICE(I) = .FALSE. 
          BAS(I)=NMATT+I 
          DO J=1,NMATT
              B(I,J) = ZERO
          END DO
          B(I,I)=-ONE
          XI(I) = ZERO
      END DO	
!	  
! 	Determine initial values
!	  	  
      DO I=1,NMATT  ! XSUB=-(B\QU)
	     XSUB(I) = QU(I)
           DO J=1,NMATT
               BSTOR(I,J)=B(I,J)
           END DO
      END DO  ! Initialisation de XSUB
      CALL DGESV(NMATT,1,BSTOR,NMATT,IPVT,XSUB,NMATT,INFO)  
	  ! Now XSUB (copy of -QU) has been overwritten by the solution XSUB	 
      DO I=1,NMATT
          XSUB(I) = -XSUB(I)
      END DO
!
!	Check if initial basis provides solution 
!
      IF (MINVAL(XSUB) .GE. ZERO) THEN
          DO I=1,NMATT
              XI(I) = ZERO
          END DO
          RETURN
      END IF		  
!
      T = 2*NMATT+1    ! Artificial variable
      ENTER = T        ! is the first entering variable
!
!	Determine initial leaving variable	  
!
      DO I=1,NMATT
	     MX(I) = -XSUB(I)
      END DO
      TVAL = MAXVAL(MX) 
      LVIND = MAXLOC(MX,1)
      LEAVE = BAS(LVIND)
      BAS(LVIND) = T   ! Pivot in the artificial variable
      DO I=1,NMATT
         XSUB(I) = XSUB(I)+TVAL
         ONES(I) = ONE
         BE(I) = ZERO
      END DO		 
      XSUB(LVIND) = TVAL	 
      BONES = MATMUL(B,ONES)
      DO I=1,NMATT
         B(I,LVIND) = -BONES(I)
      END DO		 
!
!   Main iterations begin here
!	 
      DO ITER=1,MAXITR	 
         IF (LEAVE .EQ. T) THEN
               B(1,1)=0
		     EXIT  !  Check if done; if not, get new entering variable
         ELSE IF (LEAVE .LE. NMATT) THEN
             ENTER = NMATT + LEAVE
             DO I=1,NMATT
                 BE(I) = ZERO
             END DO
             BE(LEAVE) = -ONE 
         ELSE
             ENTER = LEAVE - NMATT
             DO I = 1,NMATT
                 BE(I) = MATRXX(I,ENTER)
             END DO
         END IF
         DO I=1,NMATT  ! D = B\BE
             D(I) = BE(I)
             DO J=1,NMATT
                 BSTOR(I,J) = B(I,J)
             END DO
         END DO  ! Initialisation de D
         CALL DGESV(NMATT,1,BSTOR,NMATT,IPVT,D,NMATT,INFO)  
		 ! Now D (copy of BE) has been overwritten by the solution D	
!
!		 Find new leaving variable
!
         DO I=1,NMATT   ! indices of d>0
             IF (D(I) .GT. PIVTOL) THEN
                 INDICE(I) = .TRUE.
             ELSE
                 INDICE(I) = .FALSE.
             END IF
         END DO
         DO I=1,NMATT
             IF (INDICE(I)) THEN
                 THETAV(I) = (XSUB(I)+ZERTOL)/D(I)
             ELSE
                 THETAV(I) = HUGE(ONE)
             END IF
         END DO 
         THETA = MINVAL(THETAV) ! minimal ratios, d>0
         DO I=1,NMATT
             IF ((XSUB(I)/D(I)) .GT. THETA) THEN
                 INDICE(I) = .FALSE.
             END IF
         END DO	!  indices of minimal ratios, d>0, in INDICE
         J = 0
         DO I=1,NMATT  ! check if artificial among these
             IF (INDICE(I) .AND. (BAS(I).EQ.T)) THEN
			     J = I
			 END IF
         END DO
         IF (J .GT. 0) THEN  ! Always use artifical if possible
             LVIND = J   
         ELSE   ! otherwise pick among set of max d 
            DO I=1,NMATT
               IF (INDICE(I)) THEN
                   THETAV(I) = D(I)
               ELSE
                   THETAV(I) = -HUGE(ONE)
               END IF
            END DO 
            LVIND = MAXLOC(THETAV,1)
         END IF
         LEAVE = BAS(LVIND) 
!
!		Perform pivot
!        
         RATIO = XSUB(LVIND)/D(LVIND)
         DO I=1,NMATT
             XSUB(I) = XSUB(I) - RATIO*D(I)
         END DO
         XSUB(LVIND) = RATIO 
         DO I=1,NMATT
             B(I,LVIND) = BE(I)  
         END DO
         BAS(LVIND) = ENTER 
      END DO    ! end of iterations	
      IF ((ITER .GE. MAXITR) .AND. (LEAVE .NE. T)) THEN
         RETURN
      END IF
      DO I=1,NMATT
         Z(BAS(I)) = XSUB(I);
      END DO
      DO I=1,NMATT
         XI(I) = Z(I);
      END DO
      RETURN
!
      END SUBROUTINE