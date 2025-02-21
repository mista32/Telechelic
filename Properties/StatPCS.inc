!====Get_PCS: compute Distribution of Cluster Sizes
SUBROUTINE Get_PCS(Switch)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER Switch, NStat, I, J, K, IT, LK, LIT, NCluster
	INTEGER, ALLOCATABLE:: NIT(:), L(:), NumN(:)
	
	DOUBLE PRECISION SigmaA, SigmaC, Rcut, Rcut2, SumNC, NC, PCS, PCSNorm
	DOUBLE PRECISION dx, dy, dz, r2
	DOUBLE PRECISION, ALLOCATABLE:: PCluster(:), PNorm(:)
	DOUBLE PRECISION, ALLOCATABLE:: Sigma(:)
	
	SAVE NStat, Rcut2, SumNC, NIT, L, NumN, PCluster, PNorm, Sigma
!=========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
				
!		====find clusters
		DO I = 1, NIon
			L(I) = I
		END DO
				
		DO I = 1, NIon - 1
			IF (I==L(I)) THEN
				J = I
				DO K = I + 1, NIon
					LK = L(K)
					IF (LK==K) THEN
						dx = IX(J) - IX(K)
						dy = IY(J) - IY(K)
						dz = IZ(J) - IZ(K)
						dx = dx - Box*ANINT(dx/Box)
						dy = dy - Box*ANINT(dy/Box)
						dz = dz - Box*ANINT(dz/Box)
						r2 = dx**2 + dy**2 + dz**2
						IF (r2<Rcut2*(Sigma(J)+Sigma(K))**2/4) THEN
							L(K) = L(J)
							L(J) = LK
						END IF
					END IF
				END DO
				J = L(J)
10				IF (J .NE. I) THEN
					DO K = I + 1, NIon
						LK = L(K)
						IF (LK==K) THEN
							dx = IX(J) - IX(K)
							dy = IY(J) - IY(K)
							dz = IZ(J) - IZ(K)
							dx = dx - Box*ANINT(dx/Box)
							dy = dy - Box*ANINT(dy/Box)
							dz = dz - Box*ANINT(dz/Box)
							r2 = dx**2 + dy**2 + dz**2
							IF (r2<Rcut2*(Sigma(J)+Sigma(K))**2/4) THEN
								L(K) = L(J)
								L(J) = LK
							END IF
						END IF
					END DO
					J = L(J)
					GOTO 10
				END IF
			END IF
		END DO
				
!		====count number in a cluster containing particle IT
		DO IT = 1, NIon
			NIT(IT) = 1
			LIT = L(IT)
20			IF (LIT .NE. IT) THEN
				NIT(IT) = NIT(IT) + 1
				LIT = L(LIT)
				GOTO 20
			END IF
		END DO
				
!		====count number of total particles for a given cluster size  done
		DO I = 1, NIon
			NumN(I) = 0
		END DO
		DO IT = 1, NIon
			DO I = 1, NIon
				IF (NIT(IT)==I) THEN
					NumN(I) = NumN(I) + 1
				END IF
			END DO
		END DO

!		====distribution of cluster sizes
		DO I = 1, NIon
			PCluster(I) = PCluster(I) + NumN(I)/I
			!==========Sum number of clusters   12.14
		END DO
		
!		====average cluster size				
		NCluster = 0
		DO I = 1, NIon
			NCluster = NCluster + NumN(I)/I
		END DO
		SumNC = SumNC + DBLE(NIon)/DBLE(NCluster)
		
!		====normalized distribution of cluster sizes
		DO I = 1, NIon
			PNorm(I) = PNorm(I) + (NumN(I)/I)/DBLE(NCluster)
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		NStat = 0
		SumNC = 0.0
		
!		====read cutcoff
		OPEN (UNIT=21, FILE='Rcut.txt')
		READ (21, *) Rcut
		CLOSE (21)
		Rcut2 = Rcut**2
		
		OPEN (UNIT=23, FILE='../Conf/Part_Size')
		READ (23, *) SigmaA, SigmaC
		CLOSE (23)
		
!		===allocate arrays
		ALLOCATE(NIT(NIon), L(NIon), NumN(NIon), PCluster(NIon), PNorm(NIon), Sigma(NIon))
		
		DO I = 1, NIon
			IF (PCharge(I)==-1) THEN
				Sigma(I) = SigmaA
			ELSE
				Sigma(I) = SigmaC
			END IF
		END DO
		
!		====initilize distribution of cluster sizes
		DO I = 1, NIon
			PCluster(I) = 0.0
			PNorm(I) = 0.0 
		END DO

		WRITE (12, "('  To calculate PCS')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=35, FILE='PCluster')
		DO I = 1, NIon
			PCS = PCluster(I)/(DBLE(NStat))
			PCSNorm = PNorm(I)/(DBLE(NStat))
			WRITE (35, "(10E17.8)") DBLE(I), PCS, PCSNorm
		END DO
		CLOSE (35)
		
		NC = SumNC/(DBLE(NStat))
		OPEN (UNIT=36, FILE='Cluster')
		WRITE (36, "(10E17.8)") NC
		CLOSE (36)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NIT, L, NumN, PCluster, PNorm, Sigma)
	END IF

RETURN
END
