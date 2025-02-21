!====Get_PSP: compute Polymer-Specific Properties
SUBROUTINE Get_PSP(Switch)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Switch, NStat, I, J, P, B, B1, B2
	
	DOUBLE PRECISION dx, dy, dz, r2, L2, SumRe2, SumRg2
	DOUBLE PRECISION BCF, RE2, RG2
	DOUBLE PRECISION, ALLOCATABLE:: Acc(:), Sum(:), Fbb(:)
	DOUBLE PRECISION, ALLOCATABLE:: LX(:), LY(:), LZ(:)
	DOUBLE PRECISION, ALLOCATABLE:: UX(:), UY(:), UZ(:)
	
	SAVE NStat, Acc, Sum, Fbb, SumRe2, SumRg2
	SAVE LX, LY, LZ, UX, UY, UZ
!=========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
		DO P = 1, NPoly
!			====set bond vectors and unit bond vectors
			DO B = 1, NBead - 1
				I = (P - 1)*NBead + B
				LX(B) = BX(I+1) - BX(I)
				LY(B) = BY(I+1) - BY(I)
				LZ(B) = BZ(I+1) - BZ(I)
				L2 = LX(B)**2 + LY(B)**2 + LZ(B)**2
				UX(B) = LX(B)/L2**0.5
				UY(B) = LY(B)/L2**0.5
				UZ(B) = LZ(B)/L2**0.5
			END DO
			
!			====bond correlation function
			DO B = 1, NBead - 1
				Sum(B) = 0
				Acc(B) = 0
			END DO
			
			DO B1 = 1, NBead - 1
				DO B2 = B1, NBead - 1
					B = B2 - B1 + 1
					Sum(B) = Sum(B) + UX(B1)*UX(B2) + UY(B1)*UY(B2) + UZ(B1)*UZ(B2)
					Acc(B) = Acc(B) + 1
				END DO
			END DO
			
			DO B = 1, NBead - 1
				Fbb(B) = Fbb(B) + Sum(B)/Acc(B)
			END DO
			
!			====end-to-end distance
			dx = BX(P*NBead) - BX((P-1)*NBead + 1)
			dy = BY(P*NBead) - BY((P-1)*NBead + 1)
			dz = BZ(P*NBead) - BZ((P-1)*NBead + 1)
			r2 = dx**2 + dy**2 + dz**2
			SumRe2 = SumRe2 + r2
			
!			====radius of gyration
			DO B = 1, NBead
				I = (P - 1)*NBead + B
				dx = BX(I) - PX(P)
				dy = BY(I) - PY(P)
				dz = BZ(I) - PZ(P)
				r2 = dx**2 + dy**2 + dz**2
				SumRg2 = SumRg2 + r2
			END DO
		END DO
		
	ELSE IF (Switch==0) THEN
!		====initialization
		NStat = 0
		
		SumRe2 = 0.0
		SumRg2 = 0.0

!		====allocate arrays
		ALLOCATE(Acc(NBead), Sum(NBead), Fbb(NBead))
		
		DO I = 1, NBead
			Fbb(I) = 0.0
		END DO
		
!		====allocate arrays
		ALLOCATE(LX(NBead), LY(NBead), LZ(NBead))
		ALLOCATE(UX(NBead), UY(NBead), UZ(NBead))
		
		WRITE (12, "('  To calculate PSP')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=26, FILE='BCF')
		DO I = 1, NBead - 1
			BCF = Fbb(I)/(DBLE(NStat)*DBLE(NPoly))
			WRITE (26, "(10E17.8)") DBLE(I - 1), BCF
		END DO
		CLOSE (26)
		
		RE2 = SumRe2/(DBLE(NStat)*DBLE(NPoly))
		RG2 = SumRg2/(DBLE(NStat)*DBLE(NPoly)*DBLE(NBead))
		
		OPEN (UNIT=27, FILE='Size')
		WRITE (27, "(10E17.8)") RE2, RG2
		CLOSE (27)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(Acc, Sum, Fbb)
		DEALLOCATE(LX, LY, LZ)
		DEALLOCATE(UX, UY, UZ)
	END IF

RETURN
END
