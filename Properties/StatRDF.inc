!====Get_RDF: compute Radial Distribution Function
SUBROUTINE Get_RDF(Switch)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER Switch, NStat, I, J, NBin, Bin, Natom(0:3)
	
	DOUBLE PRECISION Rmax, DelR, dx, dy, dz, r2
	DOUBLE PRECISION RDF00, RDFii, RDFaa, RDFcc, RDFac
	DOUBLE PRECISION, ALLOCATABLE:: RBin(:), VBin(:), G00(:)
	DOUBLE PRECISION, ALLOCATABLE:: Gii(:), Gaa(:), Gcc(:), Gac(:)

	SAVE Natom, NStat, NBin, DelR, RBin, VBin, G00, Gii, Gaa, Gcc, Gac
!=========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
!		====whole system
		DO I = 1, NPart
			DO J = I + 1, NPart
				dx = RX(I) - RX(J)
				dy = RY(I) - RY(J)
				dz = RZ(I) - RZ(J)
				dx = dx - Box*ANINT(dx/Box)
				dy = dy - Box*ANINT(dy/Box)
				dz = dz - Box*ANINT(dz/Box)
				r2 = dx**2 + dy**2 + dz**2
				Bin = INT(r2**0.5/DelR)
				IF (Bin<=NBin) THEN
					G00(Bin) = G00(Bin) + 1
				END IF
			END DO
		END DO
		
!		====charged beads				
		DO I = 1, NIon
			DO J = I + 1, NIon
				dx = IX(I) - IX(J)
				dy = IY(I) - IY(J)
				dz = IZ(I) - IZ(J)
				dx = dx - Box*ANINT(dx/Box)
				dy = dy - Box*ANINT(dy/Box)
				dz = dz - Box*ANINT(dz/Box)
				r2 = dx**2 + dy**2 + dz**2
				Bin = INT(r2**0.5/DelR)
				IF (Bin<=NBin) THEN
!					====ion-ion
					Gii(Bin) = Gii(Bin) + 1
!					====anion-anion
					IF (PCharge(I)==-1 .AND. PCharge(J)==-1) THEN
						Gaa(Bin) = Gaa(Bin) + 1
					END IF
!					====cation-cation
					IF (PCharge(I)==1 .AND. PCharge(J)==1) THEN
						Gcc(Bin) = Gcc(Bin) + 1
					END IF
!					====anion-cation
					IF (PCharge(I) .NE. PCharge(J)) THEN
						Gac(Bin) = Gac(Bin) + 1
					END IF
				END IF
			END DO
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Natom(0) = NPart
		Natom(1) = NIon
		Natom(2) = NAnion
		Natom(3) = NCation
		
		NStat = 0
		
!		====set basic info
		Rmax = 0.50D0*Box
		DelR = 0.01D0
		
!		====number of bins
		NBin = INT(Rmax/DelR)
		
!		====allocate arrays
		ALLOCATE(RBin(0:NBin), VBin(0:NBin), G00(0:NBin))
		ALLOCATE(Gii(0:NBin), Gaa(0:NBin), Gcc(0:NBin), Gac(0:NBin))

!		====bin position and shell volume
		DO I = 0, NBin
			RBin(I) = (I + 0.5D0)*DelR
			VBin(I) = 4*PI*DelR**3*(DBLE(I+1)**3 - DBLE(I)**3)/3
			G00(I) = 0.0
			Gii(I) = 0.0
			Gaa(I) = 0.0
			Gcc(I) = 0.0
			Gac(I) = 0.0
		END DO
		
		WRITE (12, "('  To calculate RDF')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=25, FILE='RDF')
		DO I = 0, NBin - 1
			RDF00 = 2*Vol*G00(I)/(VBin(I)*DBLE(NStat)*Natom(0)*Natom(0))
			RDFii = 2*Vol*Gii(I)/(VBin(I)*DBLE(NStat)*Natom(1)*Natom(1))
			RDFaa = 2*Vol*Gaa(I)/(VBin(I)*DBLE(NStat)*Natom(2)*Natom(2))
			RDFcc = 2*Vol*Gcc(I)/(VBin(I)*DBLE(NStat)*Natom(3)*Natom(3))
			RDFac = 1*Vol*Gac(I)/(VBin(I)*DBLE(NStat)*Natom(2)*Natom(3))
			WRITE (25, "(10E17.8)") RBin(I), RDF00, RDFii, RDFaa, RDFcc, RDFac
		END DO
		CLOSE (25)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(RBin, VBin, G00)
		DEALLOCATE(Gii, Gaa, Gcc, Gac)
	END IF

RETURN
END
