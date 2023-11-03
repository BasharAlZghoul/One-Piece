! ************************************************************************************
!                     Particle Tracking Code For Colloid Transport                   *                   
!                               in Complex Porous Media                              *
!					(One Piece)				     *	
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

	INTEGER,PARAMETER :: Ny = 501			! Number of Cells in the x direction (consistent with input velocity field)       
	INTEGER,PARAMETER :: Nx = 501                   ! Number of Cells in the y direction 										
        INTEGER,PARAMETER :: Np = 192000		! Number of Particles injected	(consistent with input colloid locations)	
	
	INTEGER,PARAMETER :: IS = 20 			! Ionic Strength in mMolar
	INTEGER,PARAMETER :: case_num = 9		! As this number increase (becomes more unfavorable) see functions
	INTEGER,PARAMETER :: size1 = 50			! Grain smallest size
	INTEGER,PARAMETER :: size2 = 75			! Grain medium size
	INTEGER,PARAMETER :: size3 = 150		! Grain max size (This code can handle three sizes of grains in the media)
        INTEGER,PARAMETER :: tf = 100d+6		! The final time in seconds

	REAL(8), PARAMETER :: v_multiplier = 2.5847/4.0	! velocity nultiplier (to change velocity as needed if differs from input)
	REAL(8), PARAMETER :: D_multiplier = 3.0	! The box size you want each colloid to look at when calculating the nearest grain i.e. 3 times the grain size
	REAL(8), PARAMETER :: num_rounds = 402.0	! number of rounds around the domain
	REAL(8), PARAMETER :: BTC_save = 200.0		! the round number to save the breakthrough curve at
	REAL(8), PARAMETER :: BTC_spatial = 125.0	! the spatial location to save the breakthrough curve at
	REAL(8), PARAMETER :: ns_zone = 200.0d-9	! near surface zone thickness 
	REAL(8), PARAMETER :: dp = 1.1d-6		! colloid diameter
	REAL(8), PARAMETER :: ap = dp/2.0		! colloid radius
	REAL(8), PARAMETER :: length_scale = 1.0d-6	! the numerical cell edge length in meter

	
	! Fluid PARAMETERS
	REAL(8), PARAMETER :: miou = 9.98d-4		! Dynamic viscosity
	REAL(8), PARAMETER :: g = 9.981			! gravitational acceleration
	REAL(8), PARAMETER :: rhof = 998.0		! fluid density
	REAL(8), PARAMETER :: rhop = 1055.0		! colloid density
	REAL(8), PARAMETER :: T = 293.15		! Temperature
	REAL(8), PARAMETER :: pi = 3.14159265358979
	REAL(8), PARAMETER :: mp = (4/3)*pi*ap**3*rhop
	REAL(8), PARAMETER :: masses = mp + (2/3) * ap**3*pi*rhof
	REAL(8), PARAMETER :: FG = -(4/3)*pi*ap**3*(rhop-rhof)*g
	REAL(8), PARAMETER :: Boltz = 1.380649d-23
	REAL(8), PARAMETER :: con = pi*miou*ap
	REAL(8), PARAMETER :: dtmrt = mp/(6*pi*miou*ap)

	! DLVO PARAMETERS
	REAL(8), PARAMETER :: A123 = 3.83d-21		! Hamaker Constant
	REAL(8), PARAMETER :: eoer = 7.083d-10
	REAL(8), PARAMETER :: e_charge = 1.602176621d-19
	REAL(8), PARAMETER :: Na = 6.02214086d+23
	REAL(8), PARAMETER :: z = 1.0
	REAL(8), PARAMETER :: Debye_length = (eoer*Boltz*T/2/Na/z**2/e_charge**2/IS)**0.5
	REAL(8), PARAMETER :: k = 1.0/Debye_length
	REAL(8), PARAMETER :: RZOI = 2.0*(Debye_length*ap)**0.5
	REAL(8), PARAMETER :: AAp = 75.0
	REAL(8), PARAMETER :: AEFF = 2.0*(ap*Ap)/(ap+Ap)
	REAL(8), PARAMETER :: LAMBDAAB = 6.0d-10	
	REAL(8), PARAMETER :: GAMMA0AB = -0.0270305	
	REAL(8), PARAMETER :: SIGMAC = 5.0d-10	
	REAL(8), PARAMETER :: lambda = 100.0d-9		! Van Der Waals length of interaction
	REAL(8), PARAMETER :: aSte = 5.0d-8	
	REAL(8), PARAMETER :: gammaoSte = 1.7d-2	
	REAL(8), PARAMETER :: lambdaSte = 4.1d-10

	
	REAL(8), ALLOCATABLE, DIMENSION(:,:) ::  D
        REAL(8), ALLOCATABLE, DIMENSION(:) ::  yp, xp, xp0, yp0, dt, time, ang, ns_counter_prev, ns_counter
	REAL(8), ALLOCATABLE, DIMENSION(:) ::  ns_counter1, ns_counter2, ns_counter3, ns_counter4, nsg_counter
	REAL(8), ALLOCATABLE, DIMENSION(:) ::  FBn, FBt, x_enter_ns, y_enter_ns, Rp, ratios, ang_enter_ns
	REAL(8), ALLOCATABLE, DIMENSION(:) ::  round_counter, h, nnx, nny, Vx, Vy, Vn, Vtx, Vty, Un, Utx, Uty
	REAL(8), ALLOCATABLE, DIMENSION(:) ::  Ux, Uy, M, prev_grain, h_prev, prev_ang, ns_counter_acc, nearest_g 
	REAL(8), ALLOCATABLE, DIMENSION(:) ::  FVDW, FBORN, FEDL, FEDL_het, R,x_prev, BTC_time, BTC_v, f1, f2, f3, f4
	INTEGER, ALLOCATABLE, DIMENSION(:) ::  aux_ng, aux_ngs
	LOGICAL, ALLOCATABLE, DIMENSION(:) ::  pass_log, attached_log

END MODULE SHARE



! *************************************************************************************

PROGRAM PARTICLE_TRACKING

	USE SHARE, ONLY: Nx, Ny, yp, xp, xp0, yp0, ap, length_scale, miou, size1, size2, size3, pi, Boltz,  &
					k, ns_counter, ns_counter1, ns_counter2, ns_counter3, ns_counter4,  &
					dt, time, ang, ns_counter_prev, nsg_counter, FBn, masses, con, FBt, & 
					x_enter_ns, y_enter_ns, Rp, ratios, ang_enter_ns, T, round_counter, &
					nearest_g, h, aux_ng, nnx, nny, Vx, Vy, Vn, Vtx, Vty, Un, Utx, Uty, &
					Uy, Ux, M, dtmrt, num_rounds, pass_log, attached_log, prev_grain,   &
					prev_ang, ns_zone, ns_counter_acc, f1, f2, f3, f4, RZOI, FVDW, R,   & 
					FBORN, FEDL, FEDL_het, Debye_length, h_prev, aux_ngs, Ly, x_prev,   &
					BTC_spatial, BTC_save, BTC_time, BTC_v, D, Np, tf, v_multiplier,    &
					D_multiplier
	USE MPI_Module
        IMPLICIT NONE

	REAL(8)    :: Cc, Cp, Cc_het, het_radi_small,num_small,het_radi_large,num_large, surface_coverage
	REAL(8)    :: radi1, radi2, radi3, het_angle10_S, het_angle15_S, het_angle20_S, het_radi_ratio
	REAL(8)    :: het_angle10_L, het_angle15_L, het_angle20_L, angle_freq, min_R, Ufx(Nx,Ny), Ufy(Nx,Ny)
	REAL(8)	   :: yp_tot(1:Np), xp_tot(1:Np), output(1:Np,1:23), my_RC_max, RC_max, max_Gsize
	REAL(8)	   :: random, run_time, tc, my_RC_min, RC_min, con_ang, aux_f1, aux_f2
	INTEGER(8) :: t1, t2, Rate
	INTEGER(2) :: hour, min, sec
	INTEGER    :: i, j, n, c, Ng, size_aux_ngs, step
  	LOGICAL    :: logical_total(1:Np,1:2)


	CALL SYSTEM_CLOCK( t1, Rate )
	CALL MPI_INIT( MPI_err )
	CALL MPI_COMM_RANK( MPI_Comm_World, myProc, MPI_err )
	CALL MPI_COMM_SIZE( MPI_Comm_World, nprocs, MPI_err )

	step = INT(tf/100)
	tc = tf
	
	IF (myProc == root) THEN
		CALL Read_inlet_colloids_locations(yp_tot,xp_tot)
		CALL Read_Velocity(Ufx, Ufy)
		Ufx = Ufx * v_multiplier
		Ufy = Ufy * v_multiplier
	END IF

	CALL Read_Geometry(Ng)

	CALL Allocate_Data_MPI_2D(Np,Ng)
	CALL initialize_other_variables

	CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_Err )
	
	CALL MPI_SCATTER( yp_tot(:), Ly, MPI_DOUBLE_PRECISION, &
					  yp(1:Ly), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

	CALL MPI_SCATTER( xp_tot(:), Ly, MPI_DOUBLE_PRECISION, &
					  xp(1:Ly), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )
	
	CALL MPI_BCAST( Ufx, Nx*Ny, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_Err )
	CALL MPI_BCAST( Ufy, Nx*Ny, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_Err )


	yp0 = yp
	xp0 = xp

	max_Gsize = D_multiplier*MAXVAL(D(3,:))

	
	! Shuffle the Seed
	CALL random_seed()
		
	CALL EDL_param(Cp, Cc)
	Cc_het = -Cc
	CALL het_domain(het_radi_small,num_small,het_radi_large,num_large)
	CALL cov(surface_coverage, D, het_radi_small, het_radi_large, length_scale, num_small, num_large)

	radi1 = size1 * length_scale
	radi2 = size2 * length_scale
	radi3 = size3 * length_scale

	IF ( het_radi_small > 0.0	) THEN

		het_radi_ratio = het_radi_small
	END IF

	IF ( het_radi_large > 0.0	) THEN

		het_radi_ratio = het_radi_large
	END IF


	het_angle10_S = (het_radi_small/(size1*length_scale))*(180.0/pi)
	het_angle15_S = (het_radi_small/(size2*length_scale))*(180.0/pi)
	het_angle20_S = (het_radi_small/(size3*length_scale))*(180.0/pi)


	het_angle10_L = (het_radi_large/(size1*length_scale))*(180.0/pi)
	het_angle15_L = (het_radi_large/(size2*length_scale))*(180.0/pi)
	het_angle20_L = (het_radi_large/(size3*length_scale))*(180.0/pi)
	angle_freq = 360.0/num_large	


	! Perform First time Step to initialize variables
	DO i = 1, Ly


		aux_ngs = PACK([(j, j=1,Ng)], D(1,:) >= xp(i)-max_Gsize .AND. D(1,:) <= xp(i) + max_Gsize .AND. D(2,:) >= yp(i)-max_Gsize .AND. D(2,:) <= yp(i) + max_Gsize )
		size_aux_ngs = SIZE(aux_ngs)
		DO c = 1, size_aux_ngs
			R(c) = DSQRT( (xp(i) - D(1,aux_ngs(c)))**2.0 + (yp(i) - D(2,aux_ngs(c)))**2.0 	) - D(3,aux_ngs(c))
		END DO

		min_R = MINVAL(R(1:size_aux_ngs))
		aux_ng = PACK([(j, j=1,size_aux_ngs)], R(1:size_aux_ngs) == min_R)


		IF (min_R == 0.0) THEN
			PRINT*, 'ERROR: problem in finding nearest grain '
		END IF

		h(i) = min_R*length_scale
		nearest_g(i) = aux_ngs(aux_ng(1))
		Rp(i) = D(3,nearest_g(i))*length_scale
		nnx(i) = (xp(i) - D(1,nearest_g(i)))/(min_R+Rp(i)/length_scale)
		nny(i) = (yp(i) - D(2,nearest_g(i)))/(min_R+Rp(i)/length_scale)


		CALL BILINEAR_INTERPOLATION_VEL(Ufx, Ufy, Vx(i), Vy(i), xp(i), yp(i))
	

		Vn(i)  = Vx(i)*nnx(i) + Vy(i)*nny(i)
		Vtx(i) = Vx(i) - Vn(i)*nnx(i)
		Vty(i) = Vy(i) - Vn(i)*nny(i)
		Un(i)  = Vn(i)
		Utx(i) = Vtx(i)
		Uty(i) = Vty(i)
		Ux(i)  = Un(i)*nnx(i) + Utx(i)
		Uy(i)  = Un(i)*nny(i) + Uty(i)
		M(i)   = DSQRT(Ux(i)**2.0 + Uy(i)**2.0)
		x_prev(i) = xp(i)
		xp(i)  = (xp(i)*length_scale + Ux(i)*dt(i))/length_scale
		yp(i)  = (yp(i)*length_scale + Uy(i)*dt(i))/length_scale
		time(i) = time(i) + dt(i)
		

	END DO


Print*, 'Start'


	! The Main loop
	DO n = 1, tf
	

			IF (n == 1 .OR. MOD(n,step) == 0 .OR. n == tf) THEN
			

			

			CALL MPI_GATHER(  xp(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
								output(:,1), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( yp(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,2), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( time(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,3), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )
			

			CALL MPI_GATHER( nsg_counter(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,4), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )
			

			CALL MPI_GATHER( ns_counter(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,5), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( ns_counter1(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,6), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( ns_counter2(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,7), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( ns_counter3(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,8), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( ns_counter4(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,9), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( ns_counter_acc(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,10), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( h_prev(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,11), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( h(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,12), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( ang(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,13), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( nearest_g(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,14), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( dt(1:Ly)/dtmrt, Ly, MPI_DOUBLE_PRECISION, &
							   output(:,15), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( x_enter_ns(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,16), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )


			CALL MPI_GATHER( y_enter_ns(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,17), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( ang_enter_ns(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,18), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( round_counter(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,19), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( xp0(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,20), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( yp0(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,21), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( BTC_time(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,22), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( BTC_v(1:Ly), Ly, MPI_DOUBLE_PRECISION, &
							   output(:,23), Ly, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, MPI_err )




			CALL MPI_GATHER( pass_log(1:Ly), Ly, MPI_LOGICAL, &
							   logical_total(:,1), Ly, MPI_LOGICAL, root, MPI_COMM_WORLD, MPI_err )

			CALL MPI_GATHER( attached_log(1:Ly), Ly, MPI_LOGICAL, &
							   logical_total(:,2), Ly, MPI_LOGICAL, root, MPI_COMM_WORLD, MPI_err )



			my_RC_max = MAXVAL(round_counter)
			CALL MPI_REDUCE( my_RC_max, RC_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, MPI_COMM_WORLD, MPI_err )

			my_RC_min = MINVAL(round_counter)
			CALL MPI_REDUCE( my_RC_min, RC_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, root, MPI_COMM_WORLD, MPI_err )



				IF (myProc == root ) THEN


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
			
							DO c = 1, 23
								WRITE(121,*)  output(j,c)
						
							END DO

						END DO
			

					END IF

				END IF
	
			END IF 




		! Loop over the colloids
		DO i = 1, Ly

		IF (pass_log(i) .OR. attached_log(i)) CYCLE

		
		IF (xp(i) > Nx	.AND.  x_prev(i) <= Nx) THEN
			round_counter(i) = round_counter(i) + 1.0
		END IF

		IF (round_counter(i) > num_rounds) THEN
			pass_log(i) = .TRUE.
			CYCLE
		END IF 

		
		IF (xp(i) > Nx) THEN
			xp(i) = xp(i) - Nx + 1.0 

		END IF

		IF (yp(i) > Ny) THEN
			yp(i) = yp(i) - Ny + 1.0 
		END IF

		IF (yp(i) < 1.0) THEN
			yp(i) = yp(i) + Ny - 1.0 
		END IF

		IF (xp(i) < 1.0) THEN
			xp(i) = xp(i) + Nx - 1.0 
			round_counter(i) = round_counter(i) - 1.0
		END IF


		IF (xp(i) > BTC_spatial .AND. xp(i) < (BTC_spatial+1.0) .AND. round_counter(i) == BTC_save) THEN

			BTC_time(i) = time(i)
			BTC_v(i) = M(i)

		END IF
		
		ns_counter_prev(i) = ns_counter(i)	
		prev_grain(i) = nearest_g(i)
		h_prev(i) = h(i)
		prev_ang(i) = ang(i)
		

		aux_ngs = PACK([(j, j=1,Ng)], D(1,:) >= xp(i)-max_Gsize .AND. D(1,:) <= xp(i) + max_Gsize .AND. D(2,:) >= yp(i)-max_Gsize .AND. D(2,:) <= yp(i) + max_Gsize)

		size_aux_ngs = SIZE(aux_ngs)

		DO c = 1, size_aux_ngs
			R(c) = DSQRT( (xp(i) - D(1,aux_ngs(c)))**2.0 + (yp(i) - D(2,aux_ngs(c)))**2.0 	) - D(3,aux_ngs(c))

		END DO



		min_R = MINVAL(R(1:size_aux_ngs))

		aux_ng = PACK([(j, j=1,size_aux_ngs)], R(1:size_aux_ngs) == min_R)


		IF (min_R == 0.0) THEN
			PRINT*, 'ERROR: problem in finding nearest grain '
		END IF

		h(i) = min_R*length_scale

		IF (h(i) < 1.0d-9) THEN
			attached_log(i) = .TRUE.
			CYCLE
		END IF

		nearest_g(i) = aux_ngs(aux_ng(1))
		Rp(i) = D(3,nearest_g(i))*length_scale
		nnx(i) = (xp(i) - D(1,nearest_g(i)))/(min_R+Rp(i)/length_scale)
		nny(i) = (yp(i) - D(2,nearest_g(i)))/(min_R+Rp(i)/length_scale)
		
		con_ang = ACOS(ABS(nnx(i)))*180.0/pi

		IF (	nnx(i) <= 0.0 .AND. nny(i) > 0.0	) THEN 
			ang(i) = con_ang
		
		ELSEIF (	nnx(i) > 0.0 .AND. nny(i) >= 0.0	) THEN 
			ang(i) = 180.0 - con_ang
		
		ELSEIF (	nnx(i) >= 0.0 .AND. nny(i) < 0.0	) THEN 
			ang(i) = 180.0 + con_ang
		
		ELSEIF (	nnx(i) < 0.0 .AND. nny(i) <= 0.0	) THEN 
			ang(i) = 360.0 - con_ang
		END IF



		IF (	nearest_g(i) /= prev_grain(i) .AND. ns_counter(i) > 0.0	) THEN

			ns_counter(i) = 0.0
		END IF




		CALL BILINEAR_INTERPOLATION_VEL(Ufx, Ufy, Vx(i), Vy(i), xp(i), yp(i))


		Vn(i) = Vx(i)*nnx(i) + Vy(i) * nny(i)
		Vtx(i) = Vx(i) - Vn(i)*nnx(i)
		Vty(i) = Vy(i) - Vn(i)*nny(i)

		IF (	h(i) <= ns_zone	) THEN

			ns_counter(i) = ns_counter(i) + dt(i)
			ns_counter_acc(i) = ns_counter_acc(i) + dt(i)
		END IF
		
		IF (	ns_counter(i) > 0.0 .AND. ns_counter_prev(i) == 0.0	) THEN
			x_enter_ns(i) = xp(i)
			y_enter_ns(i) = yp(i)
			ang_enter_ns(i) = ang(i)
			nsg_counter(i) = nsg_counter(i) + 1.0

		END IF


		IF (	h(i) <= ns_zone .AND. nsg_counter(i) == 1.0	) THEN

			ns_counter1(i) = ns_counter1(i) + dt(i)		

		END IF

		IF (	h(i) <= ns_zone .AND. nsg_counter(i) == 2.0	) THEN

			ns_counter2(i) = ns_counter2(i) + dt(i)		

		END IF

		IF (	h(i) <= ns_zone .AND. nsg_counter(i) == 3.0	) THEN

			ns_counter3(i) = ns_counter3(i) + dt(i)		

		END IF

		IF (	h(i) <= ns_zone .AND. nsg_counter(i) == 4.0	) THEN

			ns_counter4(i) = ns_counter4(i) + dt(i)		

		END IF


		aux_f1 = h(i)/ap
		aux_f2 = (h(i)/ap)**0.5

		f1(i) = 1.0 - 0.443 * EXP(-1.299*aux_f1) - 0.5568 * EXP(-0.32*(aux_f1)**0.75)
		f2(i) = 1.0 + 1.455 * EXP(-1.2596*aux_f1) + 0.7951 * EXP(-0.56*aux_f2)
		f3(i) = 1.0 - 0.487 * EXP(-5.423*aux_f1) - 0.5905 * EXP(-37.83*aux_f2)
		f4(i) = 1.0 - 0.35 * EXP(-0.25*aux_f1) - 0.4 * EXP(-10.0*aux_f1)


		CALL het_angles(ratios(i), ang(i), het_angle10_S, het_angle15_S, het_angle20_S, & 
			het_angle10_L, het_angle15_L, het_angle20_L, angle_freq, Rp(i), radi1, & 
			radi2, radi3, het_radi_ratio)

		

		CALL DLVO(FVDW(i), FBORN(i), FEDL(i), Cc, Cp, h(i))
		CALL DLVO_het(FEDL_het(i), Cc_het, Cp, h(i))

		FEDL(i) = FEDL_het(i) * ratios(i) + (1.0 - ratios(i))* FEDL(i)
		
		IF (h(i) > 30.0*ns_zone			) THEN
			dt(i) = ns_zone*1.5/M(i)
		ELSE 
			dt(i) = ns_zone/2.0/M(i)		
		END IF

		IF (	h(i) < 20.0*Debye_length	) THEN		
			dt(i) = Debye_length/5.0/M(i) 
		END IF

		IF (	dt(i) > 1000.0*dtmrt		) THEN		
			dt(i) = 1000.0*dtmrt
		END IF

		IF (	dt(i) < 10.0*dtmrt		) THEN	
			dt(i) = 10.0*dtmrt
		END IF

		CALL random_stdnormal(random)
		FBn(i) = random*DSQRT(12.0*pi*ap*miou*Boltz*T/dt(i))
		CALL random_stdnormal(random)
		FBt(i) = random*DSQRT(12.0*pi*ap*miou*Boltz*T/dt(i))



		Un(i) = (masses*Un(i) + (FEDL(i) + FVDW(i) + FBn(i))*dt(i) & 
					+ con * 6.0*dt(i)*Vn(i)*f2(i))/ (masses + 6.0*con*dt(i)/f1(i))

		Utx(i) = (f4(i)*(masses*Utx(i) + FBt(i)*dt(i) * nny(i)) & 
					+ con * 6.0*dt(i)*Vtx(i)*f3(i))/ (f4(i)*masses + 6.0*con*dt(i))

		Uty(i) = (f4(i)*(masses*Uty(i) + FBt(i)*dt(i) * nnx(i)) & 
					+ con * 6.0*dt(i)*Vty(i)*f3(i))/ (f4(i)*masses + 6.0*con*dt(i))



		Ux(i)  = Un(i)*nnx(i) + Utx(i)
		Uy(i)  = Un(i)*nny(i) + Uty(i)
		M(i)   = DSQRT(Ux(i)**2.0 + Uy(i)**2.0)

		x_prev(i) = xp(i)
		xp(i)  = (xp(i)*length_scale + Ux(i)*dt(i))/length_scale
		yp(i)  = (yp(i)*length_scale + Uy(i)*dt(i))/length_scale
		
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
! **************************************************************************************
SUBROUTINE het_angles(ratios, ang, het_angle10_S, het_angle15_S, het_angle20_S, & 
			het_angle10_L, het_angle15_L, het_angle20_L, angle_freq, Rp, radi1, & 
			radi2, radi3, het_radi_ratio)
	USE SHARE, ONLY: RZOI
	IMPLICIT NONE
	REAL(8), INTENT(IN)  :: ang, het_angle10_S, het_angle15_S, het_angle20_S
	REAL(8), INTENT(IN)  :: het_angle10_L, het_angle15_L, het_angle20_L, angle_freq, Rp, radi1
	REAL(8), INTENT(IN)  :: radi2, radi3, het_radi_ratio
	REAL(8), INTENT(OUT) :: ratios
	REAL(8) :: base_angle_L, aux_ang_L, base_angle_S, aux_ang_S
	
	base_angle_L = 10.0*INT(ang/10.0)
	aux_ang_L = ang - base_angle_L

	base_angle_S = INT(ang)
	aux_ang_S = ang - base_angle_S
	
	IF ((((	aux_ang_S > 0.1 .AND. aux_ang_S < (0.1 + het_angle10_S)) .OR. &
		(aux_ang_L > 9.0 .AND. aux_ang_L < (9.0 + het_angle10_L) .AND. & 
		MOD(base_angle_L, angle_freq) == 0.0)) .AND. Rp == radi1) .OR. &
		(((aux_ang_S > 0.45 .AND. aux_ang_S < (0.45 + het_angle15_S)) .OR. &
		(aux_ang_L > 6.0 .AND. aux_ang_L < (6.0 + het_angle15_L) .AND. & 
		MOD(base_angle_L, angle_freq) == 0.0)) .AND. Rp == radi2) .OR. &
 		(((aux_ang_S > 0.75 .AND. aux_ang_S < (0.75 + het_angle20_S)) .OR. &
		(aux_ang_L > 3.0 .AND. aux_ang_L < (3.0 + het_angle20_L) .AND. & 
		MOD(base_angle_L, angle_freq) == 0.0)) .AND. Rp == radi3)) THEN
		

		ratios = het_radi_ratio/RZOI
		IF (	ratios >1.0) THEN
			ratios = 1.0
		END IF
	ELSE
		ratios = 0.0
	END IF


END 
! ********************************************************************************************
SUBROUTINE initialize_other_variables
	USE SHARE, ONLY: dt, dtmrt, round_counter, pass_log, attached_log, & 
				time, ns_counter, ns_counter1, ns_counter2, & 
				ns_counter3, ns_counter4, ns_counter_prev, prev_grain, &
				 h_prev, prev_ang, ang, x_enter_ns, y_enter_ns, ang_enter_ns, &
				 nsg_counter, ns_counter_acc, f1, f2, f3, f4, ratios, nsg_counter, &
				FVDW, FBORN, FEDL, FEDL_het, aux_ngs, yp0,yp, x_prev, BTC_time, BTC_v
	IMPLICIT NONE


	dt = 1000.0 *dtmrt ! initialize time step
	round_counter = 1.0
	pass_log = .FALSE.
	attached_log = .FALSE.
	time = 0.0
	ns_counter = 0.0
	ns_counter1 = 0.0
	ns_counter2 = 0.0
	ns_counter3 = 0.0
	ns_counter4 = 0.0
	ns_counter_acc = 0.0
	ns_counter_prev = 0.0
	prev_grain = 0.0
	h_prev = 0.0
	prev_ang = 0.0
	ang = 0.0
	x_enter_ns = 0.0
	y_enter_ns = 0.0
	ang_enter_ns = 0.0
	f1 = 0.0
	f2 = 0.0
	f3 = 0.0
	f4 = 0.0
	ratios = 0.0
	FVDW = 0.0
	FBORN = 0.0
	FEDL = 0.0
	FEDL_het = 0.0
	aux_ngs = 0.0
	nsg_counter = 0.0
	x_prev = 0.0
	BTC_time = 0.0
	BTC_v = 0.0
	
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

! ************************************************************************************

SUBROUTINE BILINEAR_INTERPOLATION_VEL(Ufx, Ufy, Vx, Vy, xp, yp)
	USE SHARE, ONLY:  Nx, Ny
	IMPLICIT NONE 
	
	REAL(8), INTENT(IN)  :: Ufx(1:Nx,1:Ny), Ufy(1:Nx,1:Ny), xp, yp
	REAL(8), INTENT(OUT) :: Vx, Vy
	REAL(8) :: bot_left, bot_right, y_bot, y_top, Rx1, Rx2, Ry1, Ry2, x1, x2, y1, y2	


	x1 = INT(xp)
	y1 = INT(yp)
	x2 = x1 + 1.0
	y2 = y1 + 1.0

	bot_left = (x2 - xp)/(x2 - x1)
	bot_right = (xp - x1)/(x2 - x1)
	y_bot = (y2 - yp)/(y2 - y1)
	y_top = (yp - y1)/(y2 - y1)


	Rx1 = bot_left*Ufx(x1, y1) + bot_right*Ufx(x2, y1)
	Rx2 = bot_left*Ufx(x1, y2) + bot_right*Ufx(x2, y2)

	Ry1 = bot_left*Ufy(x1, y1) + bot_right*Ufy(x2, y1)
	Ry2 = bot_left*Ufy(x1, y2) + bot_right*Ufy(x2, y2)

	Vx = y_bot * Rx1 + y_top * Rx2	
	Vy = y_bot * Ry1 + y_top * Ry2	




END

! ***************************************************************************************

SUBROUTINE Allocate_Data_MPI_2D(Np, Ng)
	USE SHARE, ONLY: ns_counter, ns_counter1, ns_counter2, ns_counter3, ns_counter4, &
					 dt, time, ang, ns_counter_prev, nsg_counter, FBn, & 
					FBt, x_enter_ns, y_enter_ns, Rp, ratios, ang_enter_ns, &
					 round_counter, nearest_g, h, aux_ng, nnx, nny, Vx, Vy, &
					Vn, Vtx, Vty, Un, Utx, Uty, Uy, Ux, M, pass_log, attached_log, &
					prev_grain, h_prev, prev_ang, aux_ng, R, & 
					ns_counter_acc, f1, f2, f3, f4, FVDW, FBORN, FEDL, FEDL_het, aux_ngs, &
					Ly, xp, xp0, yp, yp0, x_prev, BTC_time, BTC_v
	USE MPI_Module, ONLY: nprocs	
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: Np, Ng


	IF( MOD(Np,nprocs) /= 0 ) THEN
		PRINT*, 'Np must be divisible by the number of processors'
		PRINT*, 'Np, nprocs =', Np, nprocs
		PRINT*, 'program terminated'
		STOP
	END IF



	Ly = Np / nprocs


	ALLOCATE(	ns_counter(1:Ly)	)
	ALLOCATE(	ns_counter1(1:Ly)	)
	ALLOCATE(	ns_counter2(1:Ly)	)
	ALLOCATE(	ns_counter3(1:Ly)	)
	ALLOCATE(	ns_counter4(1:Ly)	)
	ALLOCATE(	ns_counter_acc(1:Ly)	)
	ALLOCATE(	dt(1:Ly)		)
	ALLOCATE(	time(1:Ly)		)
	ALLOCATE(	ang(1:Ly)		)
	ALLOCATE(	ns_counter_prev(1:Ly)	)
	ALLOCATE(	nsg_counter(1:Ly)	)
	ALLOCATE(	FBn(1:Ly)		)
	ALLOCATE(	FBt(1:Ly)		)
	ALLOCATE(	x_enter_ns(1:Ly)	)
	ALLOCATE(	y_enter_ns(1:Ly)	)
	ALLOCATE(	Rp(1:Ly)		)
	ALLOCATE(	ratios(1:Ly)		)
	ALLOCATE(	ang_enter_ns(1:Ly)	)
	ALLOCATE(	round_counter(1:Ly)	)
	ALLOCATE(	nearest_g(1:Ly)		)
	ALLOCATE(	h(1:Ly)			)
	ALLOCATE(	xp(1:Ly)		)
	ALLOCATE(	yp(1:Ly)		)
	ALLOCATE(	xp0(1:Ly)		)
	ALLOCATE(	yp0(1:Ly)		)
	ALLOCATE(	x_prev(1:Ly)		)
	ALLOCATE(	aux_ngs(1:Ng)		)
	ALLOCATE(	aux_ng(1:Ng)		)
	ALLOCATE(	R(1:Ng)		)
	ALLOCATE(	nnx(1:Ly)		)
	ALLOCATE(	nny(1:Ly)		)
	ALLOCATE(	Vx(1:Ly)		)
	ALLOCATE(	Vy(1:Ly)		)
	ALLOCATE(	Vn(1:Ly)		)
	ALLOCATE(	Vtx(1:Ly)		)
	ALLOCATE(	Vty(1:Ly)		)
	ALLOCATE(	Un(1:Ly)		)
	ALLOCATE(	Utx(1:Ly)		)
	ALLOCATE(	Uty(1:Ly)		)
	ALLOCATE(	Ux(1:Ly)		)
	ALLOCATE(	Uy(1:Ly)		)
	ALLOCATE(	M(1:Ly)			)
	ALLOCATE(	pass_log(1:Ly)		)
	ALLOCATE(	attached_log(1:Ly)	)
	ALLOCATE(	prev_grain(1:Ly)	)
	ALLOCATE(	h_prev(1:Ly)		)
	ALLOCATE(	prev_ang(1:Ly)		)
	ALLOCATE(	f1(1:Ly)		)
	ALLOCATE(	f2(1:Ly)		)
	ALLOCATE(	f3(1:Ly)		)
	ALLOCATE(	f4(1:Ly)		)
	ALLOCATE(	FVDW(1:Ly)		)
	ALLOCATE(	FBORN(1:Ly)		)
	ALLOCATE(	FEDL(1:Ly)		)
	ALLOCATE(	FEDL_het(1:Ly)		)
	ALLOCATE(	BTC_time(1:Ly)		)
	ALLOCATE(	BTC_v(1:Ly)		)


END

! *********************************************************************************************
SUBROUTINE DLVO(FVDW, FBORN, FEDL, Cc, Cp, sep)

	USE SHARE, ONLY: A123, ap, AAp, lambda, SIGMAC, eoer, k, pi
	IMPLICIT NONE
	REAL(8), INTENT(IN) :: sep, Cc, Cp
	REAL(8), INTENT(OUT) :: FVDW, FBORN, FEDL
	


	FVDW = - A123*ap*lambda*(lambda +22.232*sep)/((6.0*sep**2.0)*(lambda + 11.116*sep)**2.0)
!	FVDW = -(A123*ap*AAp/(6.0*(ap+AAp)*sep**2.0))*(lambda/(lambda+5.32*sep))
	
	FBORN = (ABS(A123)*SIGMAC**6.0/1260.0)*((7.0*ap-sep)/sep**8.0+(9.0*ap+sep)/(2.0*ap+sep)**8.0)



	FEDL = 4.0*pi*eoer*k*ap*Cp*Cc* & 
		(exp(-k*sep)/(1.0+exp(-k*sep))-(((Cp-Cc)**2.0)/(2.0*Cp*Cc))* & 
					(exp(-2.0*k*sep)/(1.0-exp(-2.0*k*sep))))





END 

! *****************************************************************************************

SUBROUTINE DLVO_het(FEDL, Cc, Cp, sep)

	USE SHARE, ONLY: ap, eoer, k, pi
	IMPLICIT NONE
	REAL(8), INTENT(IN)  :: sep, Cc, Cp
	REAL(8), INTENT(OUT) :: FEDL



	FEDL = 4.0*pi*eoer*k*ap*Cp*Cc* & 
		(exp(-k*sep)/(1.0+exp(-k*sep))-(((Cp-Cc)**2.0)/(2.0*Cp*Cc))* & 
					(exp(-2.0*k*sep)/(1.0-exp(-2.0*k*sep))))



END

! *********************************************************************************************
SUBROUTINE cov(	surface_coverage, D, het_radi, het_radi1, length_scale, num_small, num_large)
	USE SHARE, ONLY: pi
	IMPLICIT NONE

	REAL(8), INTENT(OUT) :: surface_coverage
	REAL(8), INTENT(IN)  :: D(1:3,1:50), het_radi, het_radi1, length_scale, num_small, num_large
	REAL(8) :: total_surfaces



	total_surfaces = SUM(2*pi*D(3,:)*length_scale)


	surface_coverage = (num_small*het_radi + (num_large-1.0) * het_radi1)*SIZE(D(3,:))/total_surfaces*100.0

END



! *********************************************************************************************
SUBROUTINE het_domain(het_radi,num_small,het_radi1,num_large)

	USE SHARE, ONLY: case_num, RZOI
	IMPLICIT NONE
	REAL(8), INTENT(OUT) :: het_radi,num_small,het_radi1,num_large


	IF ( case_num == 1	) THEN
		het_radi = 240.0d-9
		num_small = 360.0
		het_radi1 = 0.0d-9
		num_large = 36.0
	END IF


	IF ( case_num == 2	) THEN
		het_radi = 120.0d-9
		num_small = 360.0
		het_radi1 = 0.0d-9
		num_large = 36.0
	END IF


	IF ( case_num == 3	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 240.0d-9
		num_large = 36.0
	END IF


	IF ( case_num == 4	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 120.0d-9
		num_large = 36.0
	END IF


	IF ( case_num == 5	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 240.0d-9
		num_large = 18.0
	END IF


	IF ( case_num == 6	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 120.0d-9
		num_large = 18.0
	END IF



	IF ( case_num == 7	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 240.0d-9
		num_large = 9.0
	END IF


	IF ( case_num == 8	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 120.0d-9
		num_large = 9.0
	END IF


	IF ( case_num == 9	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 240.0d-9
		num_large = 6.0
	END IF


	IF ( case_num == 10	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 120.0d-9
		num_large = 6.0
	END IF


	IF ( case_num == 11	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 240.0d-9
		num_large = 4.0
	END IF


	IF ( case_num == 12	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 120.0d-9
		num_large = 4.0
	END IF


	IF ( case_num == 13	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 240.0d-9
		num_large = 2.0
	END IF


	IF ( case_num == 14	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 120.0d-9
		num_large = 2.0
	END IF

	IF ( case_num == 15	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 60.0d-9
		num_large = 2.0
	END IF


	IF ( case_num == 16	) THEN
		het_radi = 0.0
		num_small = 360.0
		het_radi1 = 120.0d-9
		num_large = 2.0
	END IF





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



! *******************************************************************************************
SUBROUTINE linspace(from, to, Np_pp, array)
    REAL(8), INTENT(in) :: from, to
    INTEGER, INTENT(in) :: Np_pp
    REAL(8), INTENT(out) :: array(Np_pp)
    REAL(8) :: range
    INTEGER :: i, n
    n = Np_pp

    range = to - from


    IF (n == 0) RETURN

    IF (n == 1) THEN
        array(1) = from
        RETURN
    END IF


    DO i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    END DO


END


! ******************************************************************************************

SUBROUTINE Read_inlet_colloids_locations(yp_tot,xp_tot)
	
	USE SHARE, ONLY: Np, Nx, Ny
	IMPLICIT NONE
	REAL(8), INTENT(OUT) :: yp_tot(1:Np), xp_tot(1:Np)
	CHARACTER(LEN=20) :: FileName1, FileName2
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




	WRITE(FileName2,'(A)') 'xp.txt'

	OPEN(101, file = FileName2, status = 'old', IOSTAT = I)

	IF( I > 0 ) THEN
		PRINT*, 'BASHAR'
		PRINT*, 'Obstacles File "', FileName2, '" does not exist.'
		PRINT*, '... PROGRAM TERMINATED ...'
		STOP
	END IF




	READ (101,*) ( xp_tot(X), X = 1, Np)

        CLOSE(101)


END


! ********************************************************************************************
SUBROUTINE Read_Velocity(Uffx, Uffy)
	USE SHARE, ONLY: Nx, Ny, length_scale
	IMPLICIT NONE
	REAL(8), INTENT(OUT) :: Uffx(Nx,Ny), Uffy(Nx,Ny)
	REAL(8) :: Ngo, Dx(50),Dy(50), Dr(50), velocity_conversion
	INTEGER :: I, Imax, Jmax, X, Y, Ng
	CHARACTER(LEN=20) :: FileName1, FileName2, FileName3, FileName4, FileName5, FileName6, FileName7
 

	WRITE(FileName1,'(A)') 'Ufx.txt'
	WRITE(FileName2,'(A)') 'Ufy.txt'

!!!!!!!!!!!!!!!!!!
	OPEN(20, file = FileName1, status = 'old', IOSTAT = I)
	IF( I > 0 ) THEN
		PRINT*, 'Obstacles File "', FileName1, '" does not exist.'
		PRINT*, '... PROGRAM TERMINATED ...'
		STOP
	END IF


	READ (20,*) ( (Uffx(X,Y), X=1,Nx), Y=1,Ny)
        CLOSE(20)
!!!!!!!!!!!!!!!!!!

	OPEN(30, file = FileName2, status = 'old', IOSTAT = I)
	IF( I > 0 ) THEN
		PRINT*, 'Obstacles File "', FileName2, '" does not exist.'
		PRINT*, '... PROGRAM TERMINATED ...'
		STOP
	END IF

	READ (30,*) ( (Uffy(X,Y), X=1,Nx), Y=1,Ny)
        CLOSE(30)

!!!!!!!!!!!!!!!!!!!

END

! ******************************************************************************************


SUBROUTINE Read_Geometry(Ng)
	USE SHARE, ONLY: D
	IMPLICIT NONE
	INTEGER, INTENT(OUT) :: Ng
	INTEGER :: I, X, Y
	CHARACTER(LEN=20) :: FileName3

	WRITE(FileName3,'(A)') 'D.txt'

	OPEN(40, file = FileName3, status = 'old', IOSTAT = I)

	IF( I > 0 ) THEN
		PRINT*, 'Obstacles File "', FileName3, '" does not exist.'
		PRINT*, '... PROGRAM TERMINATED ...'
		STOP
	END IF
	READ(40,*) Ng
	ALLOCATE(	D(1:3,1:Ng)	)



	READ (40,*) ( (D(X,Y), Y=1,Ng), X=1,3)


        CLOSE(40)

END

! ******************************************************************************************

! THE END








