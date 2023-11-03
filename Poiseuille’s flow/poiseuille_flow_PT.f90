! ************************************************************************************
!                     Particle Tracking Code For Colloid Transport                   *                   
!                               in Poiseuille Flow                                   *
!										     *	
!										     *
!										     *
!                            Written By: Bashar Al Zghoul		             *
!			     Last Updated by: May-31-2023			     *
!                                    						     *
! ************************************************************************************


MODULE MPI_Module
	USE MPI
	IMPLICIT NONE

	INTEGER :: myProc, nprocs
	INTEGER :: Status(MPI_STATUS_SIZE), MPI_err

	INTEGER :: BOT, TOP

	INTEGER, PARAMETER :: root = 0

END MODULE MPI_Module



MODULE SHARE 

        IMPLICIT NONE

	INTEGER :: Ly

        INTEGER,PARAMETER :: Np = 384000
	INTEGER,PARAMETER :: IS = 6 
        INTEGER,PARAMETER :: tf = 100d+6

	REAL(8),PARAMETER :: channel_width = 100.0d-6
	REAL(8),PARAMETER :: U0 = 3.7037d-4/4.0
	REAL(8),PARAMETER :: stop_at = 0.2	
	REAL(8), PARAMETER :: ns_zone = 200.0d-9
	REAL(8), PARAMETER :: dp = 1.1d-6
	REAL(8), PARAMETER :: ap = dp/2.0
	
	! Fluid PARAMETERS
	REAL(8), PARAMETER :: miou = 9.98d-4
	REAL(8), PARAMETER :: g = 9.981
	REAL(8), PARAMETER :: rhof = 998.0
	REAL(8), PARAMETER :: rhop = 1055.0
	REAL(8), PARAMETER :: T = 293.15
	REAL(8), PARAMETER :: pi = 3.14159265358979
	REAL(8), PARAMETER :: mp = (4/3)*pi*ap**3*rhop
	REAL(8), PARAMETER :: masses = mp + (2/3) * ap**3*pi*rhof
	REAL(8), PARAMETER :: FG = -(4/3)*pi*ap**3*(rhop-rhof)*g
	REAL(8), PARAMETER :: Boltz = 1.380649d-23
	REAL(8), PARAMETER :: con = pi*miou*ap
	REAL(8), PARAMETER :: dtmrt = mp/(6*pi*miou*ap)

	! DLVO PARAMETERS
	REAL(8), PARAMETER :: A123 = 3.83d-21
	REAL(8), PARAMETER :: eoer = 7.083d-10
	REAL(8), PARAMETER :: e_charge = 1.602176621d-19
	REAL(8), PARAMETER :: Na = 6.02214086d+23
	REAL(8), PARAMETER :: z = 1.0
	REAL(8), PARAMETER :: Debye_length = (eoer*Boltz*T/2/Na/z**2/e_charge**2/IS)**0.5
	REAL(8), PARAMETER :: k = 1.0/Debye_length;
	REAL(8), PARAMETER :: RZOI = 2.0*(Debye_length*ap)**0.5
	REAL(8), PARAMETER :: AAp = 75.0
	REAL(8), PARAMETER :: AEFF = 2.0*(ap*Ap)/(ap+Ap)
	REAL(8), PARAMETER :: LAMBDAAB = 6.0d-10	
	REAL(8), PARAMETER :: GAMMA0AB = -0.0270305	
	REAL(8), PARAMETER :: SIGMAC = 5.0d-10	
	REAL(8), PARAMETER :: lambda = 100.0d-9	
	REAL(8), PARAMETER :: aSte = 5.0d-8	
	REAL(8), PARAMETER :: gammaoSte = 1.7d-2	
	REAL(8), PARAMETER :: lambdaSte = 4.1d-10	

        REAL(8), ALLOCATABLE, DIMENSION(:) ::  yp, xp, yp0, dt, time, ns_counter, h_prev, FBn, FBt
	REAL(8), ALLOCATABLE, DIMENSION(:) ::  h, Vx, Vy, Ux, Uy, M, f1, f2, f3, f4, FVDW, FBORN, FEDL
	LOGICAL, ALLOCATABLE, DIMENSION(:) ::  pass_log, attached_log

END MODULE SHARE



! *************************************************************************************

PROGRAM PARTICLE_TRACKING

	USE SHARE, ONLY: yp, xp, yp0, ap, miou, pi, k, ns_counter, Boltz, T, dt, time, FBn, masses, con, & 
					FBt, h, Vx, Vy, Uy, Ux, M, dtmrt, pass_log, attached_log, Ly,    &
					ns_zone, f1, f2, f3, f4, RZOI, FVDW, FBORN, FEDL, Debye_length,  &
					Np, channel_width, stop_at, U0, h_prev, tf
	USE MPI_Module
        IMPLICIT NONE

	REAL(8)    :: Cc, Cp, yp_tot(1:Np), output(1:Np,1:8), random
	REAL(8)	   :: my_RC_max, RC_max, my_RC_min, RC_min, tc, run_time  
  	LOGICAL    :: logical_total(1:Np,1:2)
	INTEGER(8) :: t1, t2, Rate
	INTEGER(2) :: hour, min, sec
	INTEGER    :: i, j, n, c, step


	CALL SYSTEM_CLOCK( t1, Rate )
	CALL MPI_INIT( MPI_err )
	CALL MPI_COMM_RANK( MPI_Comm_World, myProc, MPI_err )
	CALL MPI_COMM_SIZE( MPI_Comm_World, nprocs, MPI_err )



	step = INT(tf/100)
	tc = tf

	CALL Allocate_Data_MPI_2D(Np)
	CALL initialize_other_variables
	
	IF (myProc == root) THEN
	CALL Read_inlet_colloids_locations(yp_tot)
	END IF
	
	CALL MPI_SCATTER( yp_tot(:), Ly, MPI_DOUBLE_PRECISION, &
					  yp(1:Ly), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


	yp0 = yp

	
	! Shuffle the Seed
	CALL random_seed()
	
	
	CALL EDL_param(Cp, Cc)



	DO i = 1, Ly


		IF (yp(i) <= channel_width/2.0) THEN
			h(i) = yp(i) - ap 
		ELSE
			h(i) = channel_width - yp(i) - ap
		END IF

		Vx(i) = 1.5*U0*(1 - ((yp(i)-channel_width/2.0)/channel_width*2.0)**2.0)
		Vy(i) = 0.0

		Ux(i)  = Vx(i)
		Uy(i)  = Vy(i)
		M(i)   = DSQRT(Ux(i)**2.0 + Uy(i)**2.0)
		xp(i)  = xp(i) + Ux(i)*dt(i)
		yp(i)  = yp(i) + Uy(i)*dt(i)
		time(i) = time(i) + dt(i)


	END DO


Print*, 'Start'

	! The Main loops

	DO n = 1, tf
	

			IF (n == 1 .OR. MOD(n,step) == 0) THEN
			

			

			CALL MPI_GATHER(  xp(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
								output(:,1), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( yp(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,2), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( time(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,3), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )
			

			CALL MPI_GATHER( yp0(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,4), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )
			

			CALL MPI_GATHER( ns_counter(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,5), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( h_prev(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,6), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( h(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,7), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( dt(1:Ly)/dtmrt, Ly, MPI_DOUBLE_PRECISION, &
							   output(:,8), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( pass_log(1:Ly), Ly, MPI_LOGICAL, &
							   logical_total(:,1), Ly, MPI_LOGICAL, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( attached_log(1:Ly), Ly, MPI_LOGICAL, &
							   logical_total(:,2), Ly, MPI_LOGICAL, root, MPI_COMM_WORLD, MPI_err )


			my_RC_max = MAXVAL(xp)
			CALL MPI_REDUCE( my_RC_max, RC_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, MPI_COMM_WORLD, MPI_err )

			my_RC_min = MINVAL(xp)
			CALL MPI_REDUCE( my_RC_min, RC_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, root, MPI_COMM_WORLD, MPI_err )



			IF (myProc == root ) THEN
				IF (n == 1) THEN
					OPEN (2, File = 'RUN_PROGRESS.txt', ACCESS = 'append')
					WRITE(2,'(/A6,4A12,4A10/)'), 'ts', 'min x', 'max x'
					WRITE(2,'(A)'),' *******************************'
					CLOSE(2)
					WRITE(*,'(/A6,4A12,4A10/)'), 'ts', 'min x', 'max x'
					WRITE(*,'(A)'),' *******************************'
				END IF

				OPEN (2, File = 'RUN_PROGRESS.txt', ACCESS = 'append')
				
				WRITE(2,'(F8.3,2F12.7,E12.3,2F10.1,F11.3,2F9.3)'), &
						n/tc, RC_min, RC_max
				CLOSE(2)
				WRITE(*,'(F8.3,2F12.7,E12.3,2F10.1,F11.3)'), &
						n/tc, RC_min, RC_max

          		IF ( MOD(n,tf)==0 ) THEN

			
			OPEN (121, file = 'results.dat')
			OPEN (122, file = 'log_att.dat')
			OPEN (123, file = 'log_pass.dat')
	
			DO j = 1, Np

			IF (logical_total(j,2)) THEN
				WRITE(122,*)  1
			ELSE
				WRITE(122,*)  0
			END IF 
			IF (logical_total(j,1)) THEN
				WRITE(123,*)  1
			ELSE
				WRITE(123,*)  0
			END IF 
			
				DO c = 1, 8
				WRITE(121,*)  output(j,c)
				
				END DO

				

				
			!END IF
			END DO
			

			END IF

				




			END IF



			END IF 




		!print*, n
		DO i = 1, Ly

		IF (pass_log(i) .OR. attached_log(i)) CYCLE

		IF (xp(i) >= stop_at) THEN
			pass_log(i) = .TRUE.
			CYCLE
		END IF 

				
		h_prev(i) = h(i)

		

		IF (yp(i) <= channel_width/2.0) THEN
			h(i) = yp(i)
		ELSE
			h(i) = channel_width - yp(i)
		END IF



		IF (	h(i) < 1.0d-9 	) THEN
			attached_log(i) = .TRUE.
			CYCLE
		END IF
			

		Vx(i) = 1.5*U0*(1 - ((yp(i)-channel_width/2.0)/channel_width*2.0)**2.0)
		Vy(i) = 0.0


		IF (	h(i) <= ns_zone	) THEN

			ns_counter(i) = ns_counter(i) + dt(i)

		END IF
		
		
		f1(i) = 1.0 - 0.443 * EXP(-1.299*h(i)/ap) - 0.5568 * EXP(-0.32*(h(i)/ap)**0.75)
		f2(i) = 1.0 + 1.455 * EXP(-1.2596*h(i)/ap) + 0.7951 * EXP(-0.56*(h(i)/ap)**0.5)
		f3(i) = 1.0 - 0.487 * EXP(-5.423*h(i)/ap) - 0.5905 * EXP(-37.83*(h(i)/ap)**0.5)
		f4(i) = 1.0 - 0.35 * EXP(-0.25*h(i)/ap) - 0.4 * EXP(-10.0*h(i)/ap)

		

		CALL DLVO(FVDW(i), FBORN(i), FEDL(i), -Cc, Cp, h(i))


		IF (yp(i) > channel_width/2.0) THEN
			FEDL(i) = FEDL(i)*-1.0
			FVDW(i) = FVDW(i)*-1.0
		END IF

		

		
		IF (h(i) > 30.0*ns_zone) THEN
			dt(i) = ns_zone*1.5/M(i)

		ELSE 

			dt(i) = ns_zone/2.0/M(i)		

		END IF

		
		IF (	h(i) < 20.0*Debye_length	) THEN
			
			dt(i) = Debye_length/5.0/M(i) 

		END IF

		IF (	dt(i) > 1000.0*dtmrt	) THEN
			
			dt(i) = 1000.0*dtmrt

		END IF

		IF (	dt(i) < 10.0*dtmrt	) THEN
			
			dt(i) = 10.0*dtmrt

		END IF

		CALL random_stdnormal(random)
		FBn(i) = random*DSQRT(12.0*pi*ap*miou*Boltz*T/dt(i))
		CALL random_stdnormal(random)
		FBt(i) = random*DSQRT(12.0*pi*ap*miou*Boltz*T/dt(i))



		Uy(i) = (masses*Uy(i) + ((FEDL(i) + FVDW(i) + FBn(i))*dt(i) & 
					+ con * 6.0*dt(i)*Vy(i)*f2(i)))/ (masses + 6.0*con*dt(i)/f1(i))

		Ux(i) = ((masses*Ux(i) + FBt(i)*dt(i)) & 
					+ con * 6.0*dt(i)*Vx(i)*f3(i)/f4(i))/ (masses + 6.0*con*dt(i)/f4(i))


		M(i)   = DSQRT(Ux(i)**2.0 + Uy(i)**2.0)

		
		xp(i)  = xp(i) + Ux(i)*dt(i)
		yp(i)  = yp(i) + Uy(i)*dt(i)
		
		time(i) = time(i) + dt(i)


		END DO

	

	END DO


	IF (myProc == root) THEN
		CALL SYSTEM_CLOCK( t2 )

		OPEN (2, File = 'RUN_PROGRESS.txt', ACCESS = 'append')
		WRITE(2,'(/A)'),' *****************************************'
		WRITE(*,'(/A)'),' *****************************************'
		hour = DBLE(t2-t1)/Rate /3600
		min  = MOD( DBLE(t2-t1)/Rate, 3600. )/60
		sec  = MOD( MOD(DBLE(t2-t1)/Rate,3600.), 60.)
		WRITE(2,'(A,I5.2,A,I2.2,A,I2.2)'), ' Time Elapsed:', hour, ':', min, ':', sec
		WRITE(*,'(A,I5.2,A,I2.2,A,I2.2)'), ' Time Elapsed:', hour, ':', min, ':', sec
		WRITE(2,'(A)'),' *****************************************'
		WRITE(*,'(A)'),' *****************************************'
		CLOSE(2)
	END IF

	CALL MPI_FINALIZE( MPI_err )

END 



! ********************************************************************************************
SUBROUTINE initialize_other_variables
	USE SHARE, ONLY: dt, dtmrt, pass_log, attached_log, time, xp, ns_counter,  & 
				 h_prev, f1, f2, f3, f4, FVDW, FBORN, FEDL, yp0, yp
	IMPLICIT NONE


	dt = 1000.0 *dtmrt ! initialize time step
	pass_log = .FALSE.
	attached_log = .FALSE.
	time = 0.0
	xp = 0.0
	ns_counter = 0.0
	h_prev = 0.0
	f1 = 0.0
	f2 = 0.0
	f3 = 0.0
	f4 = 0.0
	FVDW = 0.0
	FBORN = 0.0
	FEDL = 0.0
	
END

! ***********************************************************************************
SUBROUTINE random_stdnormal(x)
	USE SHARE, ONLY: pi
	IMPLICIT NONE 
	REAL(8), INTENT(OUT) :: x
	REAL(8) :: u1, u2
	
	CALL random_stduniform(u1)
	CALL random_stduniform(u2)
	x = DSQRT(-2.0*LOG(u1))*COS(2*pi*u2)
END


! ****************************************************************************************
SUBROUTINE random_stduniform(u)
	IMPLICIT NONE 
	REAL(8), INTENT(OUT) :: u
	REAL(8) :: r
	
	CALL random_number(r)
	u = 1.0 - r
END

! ***************************************************************************************

SUBROUTINE Allocate_Data_MPI_2D(Np)
	USE SHARE, ONLY: ns_counter, dt, time, FBn, FBt, h, Vx, Vy, Uy, Ux, M,  & 
				   pass_log, attached_log, h_prev, Ly, xp,      &
				   yp, yp0, f1, f2, f3, f4, FVDW, FBORN, FEDL 				
	USE MPI_Module, ONLY: nprocs	
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: Np


	IF( MOD(Np,nprocs) /= 0 ) THEN
		PRINT*, 'Np must be divisible by the number of processors'
		PRINT*, 'Np, nprocs =', Np, nprocs
		PRINT*, 'program terminated'
		STOP
	END IF



	Ly = Np / nprocs

	ALLOCATE(	ns_counter(1:Ly)	)
	ALLOCATE(	dt(1:Ly)		)
	ALLOCATE(	time(1:Ly)		)
	ALLOCATE(	FBn(1:Ly)		)
	ALLOCATE(	FBt(1:Ly)		)
	ALLOCATE(	h(1:Ly)			)
	ALLOCATE(	xp(1:Ly)		)
	ALLOCATE(	yp(1:Ly)		)
	ALLOCATE(	yp0(1:Ly)		)
	ALLOCATE(	Vx(1:Ly)		)
	ALLOCATE(	Vy(1:Ly)		)
	ALLOCATE(	Ux(1:Ly)		)
	ALLOCATE(	Uy(1:Ly)		)
	ALLOCATE(	M(1:Ly)			)
	ALLOCATE(	pass_log(1:Ly)		)
	ALLOCATE(	attached_log(1:Ly)	)
	ALLOCATE(	h_prev(1:Ly)		)
	ALLOCATE(	f1(1:Ly)		)
	ALLOCATE(	f2(1:Ly)		)
	ALLOCATE(	f3(1:Ly)		)
	ALLOCATE(	f4(1:Ly)		)
	ALLOCATE(	FVDW(1:Ly)		)
	ALLOCATE(	FBORN(1:Ly)		)
	ALLOCATE(	FEDL(1:Ly)		)



END

! *********************************************************************************************
SUBROUTINE DLVO(FVDW, FBORN, FEDL, Cc, Cp, sep)

	USE SHARE, ONLY: A123, ap, AAp, lambda, SIGMAC, eoer, k, pi
	IMPLICIT NONE
	REAL(8), INTENT(IN)  :: sep, Cc, Cp
	REAL(8), INTENT(OUT) :: FVDW, FBORN, FEDL
	


	FVDW = - A123*ap*lambda*(lambda +22.232*sep)/((6.0*sep**2.0)*(lambda + 11.116*sep)**2.0)
!	FVDW = -(A123*ap*AAp/(6.0*(ap+AAp)*sep**2.0))*(lambda/(lambda+5.32*sep))
	
	FBORN = (ABS(A123)*SIGMAC**6.0/1260.0)*((7.0*ap-sep)/sep**8.0+(9.0*ap+sep)/(2.0*ap+sep)**8.0)



	FEDL = 4.0*pi*eoer*k*ap*Cp*Cc* & 
		(exp(-k*sep)/(1.0+exp(-k*sep))-(((Cp-Cc)**2.0)/(2.0*Cp*Cc))* & 
					(exp(-2.0*k*sep)/(1.0-exp(-2.0*k*sep))))



END 


! *****************************************************************************************************
SUBROUTINE EDL_param(Cp, Cc)
	USE SHARE, ONLY: ap, IS
	IMPLICIT NONE
	REAL(8), INTENT(OUT)  :: Cp, Cc

	IF (	IS == 6	) THEN
		
		Cc = -70.0d-3

		IF ( ap == 0.05d-6) THEN
		
			Cp = -45.3d-3
		END IF
	
		IF ( ap == 0.15d-6) THEN
		
			Cp = -18.3d-3
		END IF

		IF ( ap == 0.55d-6) THEN
		
			Cp = -65.4d-3
		END IF

		IF ( ap == 1.0d-6) THEN
		
			Cp = -29.9d-3
		END IF

		IF ( ap == 2.2d-6) THEN
		
			Cp = -65.0d-3
		END IF

		IF ( ap == 3.4d-6) THEN
		
			Cp = -10.2d-3
		END IF


	END IF



	IF (	IS == 20	) THEN
		
		Cc = -53.5d-3

		IF ( ap == 0.05d-6) THEN
		
			Cp = -35.9d-3
		END IF
	
		IF ( ap == 0.15d-6) THEN
		
			Cp = -10.5d-3
		END IF

		IF ( ap == 0.55d-6) THEN
		
			Cp = -50.1d-3
		END IF

		IF ( ap == 1.0d-6) THEN
		
			Cp = -8.2d-3
		END IF

		IF ( ap == 2.2d-6) THEN
		
			Cp = -42.8d-3
		END IF

		IF ( ap == 3.4d-6) THEN
		
			Cp = -4.5d-3
		END IF


	END IF

END

! ******************************************************************************************

SUBROUTINE Read_inlet_colloids_locations(yp_tot)
	
	USE SHARE, ONLY: Np
	IMPLICIT NONE
	REAL(8), INTENT(OUT) :: yp_tot(1:Np)
	CHARACTER(LEN=20) :: FileName1
	INTEGER :: X, I
		


	WRITE(FileName1,'(A)') 'yp.txt'

	OPEN(100, file = FileName1, status = 'old', IOSTAT = I)

	IF( I > 0 ) THEN
		PRINT*, 'BASHAR'
		PRINT*, 'Obstacles File "', FileName1, '" does not exist.'
		PRINT*, '... PROGRAM TERMINATED ...'
		STOP
	END IF

	READ (100,*) ( yp_tot(X), X = 1, Np)

        CLOSE(100)

END


! ******************************************************************************************


! THE END








