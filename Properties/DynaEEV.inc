!====Get_EEV: compute End-to-End Vector relaxation
SUBROUTINE Get_EEV(Switch)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER Switch, Nt, Nt0, ttel, t, dt, I
	INTEGER, ALLOCATABLE:: NDyn(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, EProd, EEVb0, EEVb, EEVm0, EEVm
	DOUBLE PRECISION, ALLOCATABLE:: EXb(:), EYb(:), EZb(:), EXb0(:,:), EYb0(:,:), EZb0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: SumEEVb(:)
	
	SAVE Nt, Nt0, NDyn, ttv0
	SAVE EXb, EYb, EZb, EXb0, EYb0, EZb0, SumEEVb
!=========================================================
	IF (Switch==1) THEN
!		====sample
		Nt = Nt + 1
		
!		====end-to-end vectors
		DO I = 1, NPoly
			EXb(I) = BX(I*NBead) - BX((I-1)*NBead + 1)
			EYb(I) = BY(I*NBead) - BY((I-1)*NBead + 1)
			EZb(I) = BZ(I*NBead) - BZ((I-1)*NBead + 1)
		END DO

!		====set new t=0
		IF (MOD(Nt, IT0)==0) THEN
			Nt0 = Nt0 + 1
			ttel = MOD(Nt0 - 1, T0Max) + 1
			ttv0(ttel) = Nt
			DO I = 1, NPoly
				EXb0(I, ttel) = EXb(I)
				EYb0(I, ttel) = EYb(I)
				EZb0(I, ttel) = EZb(I)
			END DO
		END IF

		DO t = 1, MIN(Nt0, T0Max)
			dt = Nt - ttv0(t) + 1
			IF (dt>0 .AND. dt<TMax) THEN
				NDyn(dt) = NDyn(dt) + 1
				DO I = 1, NPoly
					EProd = EXb(I)*EXb0(I, t) + EYb(I)*EYb0(I, t) + EZb(I)*EZb0(I, t)
					SumEEVb(dt) = SumEEVb(dt) + EProd
				END DO
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		Nt0 = 0
		
!		===allocate arrays
		ALLOCATE(NDyn(TMax), ttv0(T0Max))
		ALLOCATE(EXb(NPoly), EYb(NPoly), EZb(NPoly))
		ALLOCATE(EXb0(NPoly, T0Max), EYb0(NPoly, T0Max), EZb0(NPoly, T0Max))
		ALLOCATE(SumEEVb(TMax))
		
		DO I = 1, TMax
			NDyn(I) = 0
			SumEEVb(I) = 0.0
		END DO
		
		WRITE (12, "('  To calculate EEV_', A1)") FType
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=81, FILE='EEV_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyn(I)>0) THEN
				EEVb = SumEEVb(I)/(DBLE(NDyn(I))*DBLE(NPoly))
				IF (I==1) EEVb0 = EEVb
				WRITE (81, "(10E17.8)") MDTime, EEVb/EEVb0
			END IF
		END DO
		CLOSE (81)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyn, ttv0)
		DEALLOCATE(EXb, EYb, EZb, EXb0, EYb0, EZb0, SumEEVb)
	END IF

RETURN
END
