      subroutine adapchem(igrid,igridTot,igm1,p,t,y,ickwrk,rckwrk,wdot,
     &ispecies,jsolv)
C *********************************************************************
C     Luwi Oluwole 6/9/2004 (completed 8/2005)
C     oluwoleo@mit.edu
C     Adaptive Chemistry package to be interfaced with any CFD solver
C     for efficient simulation of reacting flows using CHEMKIN-like
C     routines
C     To use, substitute call to
C     CHEMKIN routine CKWYP(P,T,Y,ICKWRK,RCKWRK,WDOT)
C     with a call to
C     adapchem2(igrid,igridTot,igm1,p,t,y,ickwrk,rckwrk,wdot,ispecies,jsolv)
C     where 'igrid' is the current grid point address/number
C           'igridTot' is the total number of grid points
C           'igm1' is set to a value different from igrid to distinguish differnt stages of the solver
C     The additional input "ispecies"
C     is an array of binary values 0/1 where ispecies(k) = 1 if species
C     k is chemically active or ispecies(k) = 0 if species k is chemically
C     inactive (i.e. no wdot(j) is changed with change in y(k))
C     This may be useful in Jacobian evaluation by finite differencing
C     where wdot is usually evaluated at Y and revaluated at Y+deltaY(k)
C     Don't need to re-evaluate wdot if ispecies(k) = 0; simply use
C     wdot@Y+deltaY(k) = wdot@Y
C     dimension ispecies = number of species
C     Also add to CFD solver, the common block
C     common/adpchm/resid_tol, icount
C     'icount' is the iteration count that lets adapchem know when CFD
C     solver has completed all operations for the previous step and is
C     starting a new iteration (or integration) step.
C     'resid_tol' is the steady-state convergence tolerance:
C     steady state when max(|dT/dt|,|dC/dt|)<resid_tol
C     This should give you an idea of where to insert the common block
C     in your CFD solver
C     resid_tol should ONLY be changed by adapchem routine
C *********************************************************************
C
C     Returns the molar production rates of the species given the
C     pressure, temperature and mass fractions;  see Eq. (49).
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     EPS0   - Error tolerance for species.
C
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
c       March 22/2010--Luwi Oluwole:
c                      jsolv= -1 ===> load reduced models only; then return
c                      jsolv = 0 ===> don't calculate reaction rates;
c                                     only return reduced model species (binary) array
c                      jsolv = 1 ===> calculate reaction rates
C
C  END PROLOGUE
C ************************

C*****precision > double
      implicit double precision (a-h, o-z), integer (i-n)
      dimension ickwrk(*), rckwrk(*), y(*), wdot(*)
      dimension ispecies(*)
C Including Chemkin common block for array address
      include 'ckstrt.i'
C Including array boundary size related parameters
      include 'ackparams.i'
      dimension cfull(kmax), yfull(kmax), yinert(kmax,nptmax)
      dimension nur(imax,nm+1), iModel(NPTMAX), iicount(nptmax)
      dimension iCall(1)
      dimension nus(kmax,nm+1)! species binary array
      dimension modelSize(nm+1,2)
      save igrTotm1, jfindmod
      common/adpchm1/ resid_tol
      common/adpchm2/ icount
      common/adpchm3/ ii, kk, kfuel, icountm1,
     & nsteps, jprint, nprint, icountTot,
     & iSizeTot
      common/adpchm4/ rmaxfull
      common/adpchm5/ nFull,nRed
      common/adpchm6/ iModel, iicount
      common/saveNUR/ nur, modelSize
      common/saveNUS/ nus
      common/inert/ yinert
        common/adpchm7/ inewlib, iCall
!       <Luwi===9/2/2010====>
c     common /mpi/ noprocs,nidproc
      common /adpmpi/ noprocs,nidproc
      common/adpchm8/ iadapnew
      if (jsolv .eq. -1) then
!     This call is supposed to load new model library
!       if no new library available, then exit without doing anything
         if( iadapnew .eq. 0 ) return
      endif
!       </Luwi===============>


C     default : do not reset count of nFull, nRed
      ireset = 0
C     Save inert species mass fractions
      if (jsolv .eq. 0) then
         do k=1,kk
            yinert(k,igrid) = y(k)
         enddo
C     convert mass fractions to concentrations
         call ackytcp(p,t,y,ickwrk,rckwrk,cfull)
           if( igm1 .ne. igridTot ) then
C     Also, CFD solver lets me know if grid cell indexes have been reshuffled since last call,Luwi--4/20/2010
              igrTotm1 = igm1
              igm1 = igridTot
C     reset count of nRed, nFull
              ireset = 1
           endif
      elseif(jsolv.eq.1) then !model already selected for these conditions
         imodNum = iModel(igrid)
         kr = 0
         do k=1,kk
            yfull(k) = 0.0d0
            if(nus(k,imodNum) .eq. 1) then
                kr = kr + 1
                yfull(k) = y(kr)
                yinert(k,igrid) = 0.0d0 !frees up some memory
            else
                yfull(k) = yinert(k,igrid)
            endif
         enddo
         call ackytcp(p,t,yfull,ickwrk,rckwrk,cfull)

      endif ! "if jsolv=0/1" conditional
C     <initialize parameters>
C     On first call, initialize, read inputs, and save data
      if(iCall(1).ne.2) call adapchinit
     & (igridTot,icount,igrTotm1,icountm1,ii,kk,kfuel,resid_tol,
     &  jprint,nprint,nsteps,iicount,iModel,nFull,nRed,rmaxfull)
C********NOTE: OLDER VERSIONS ASSUME ARRAYS (eg. ICAL()) ARE ALWAYS INITIALIZED
C              TO 0 BY YOUR COMPILER. I.E. ICAL(I) = 0 FORALL I
C     </initialize>
C
C     is a new model needed?
      if(jsolv .eq. -1) then
C     this call is to load reduced models only; then return
         jfindmod = 1
         iwritesize = 0
      else
         call needmodel(igrid,igridTot,igrTotm1,nsteps,icount,icountm1,
     & jfindmod,iwritesize,jprint,nprint,iicount)
      endif
C     record model usage statistics to date
      if(iwritesize.eq.1) then
!Luwi--10/15/10--don't write stats! call sizewrite(nFull,nRed,iSizeTot,icountTot,ii,kk)
C     <modelLibrary>
C     skip the following if a new model library is being constructed
         if(inewlib.gt.0) goto 007
  007    continue
C     </modelLibrary>
      endif
C     reset iSizeTot and icountTot
      if(ireset .eq. 1) then
        if(jfindmod.eq.1) then
            iSizeTot = 0
            icountTot = 0
            nFull = 0
            nRed = 0
         endif
         icountm1 = icount
      endif
C     <Update models: only if necessary>
C      Find smallest valid reduced model
      if(jfindmod.eq.1) then
      call findModel(jsolv,1,
     & iCall(1), kk,ii, kfuel, T, Cfull, imodNum,nFull,
     & nRed)
      if (jsolv .eq. -1) return !Initialization call only
         iModel(igrid) = imodNum
         iSizeTot = iSizeTot + modelSize(imodNum,2)!computes avg number of species
         icountTot = icountTot + 1
      endif
C     Once models obtained for all grid points, don't search for new models
C     again until retired models or new mesh.
C     </Update models>

C     <Calculate rates using selected model>
      if(jsolv .eq. 0) then
C     this call was for determining reduced model only; AND
C     it is possible reduced model already determined on prior call
C     e.g. ODE solver reducing stepsize and repeating integration
C     after prior failed attempt
         imodNum = iModel(igrid)
         do k=1,kk
            ispecies(k) = nus(k,imodNum)
         enddo
      elseif(jsolv.eq.1) then
         imodNum = iModel(igrid)
         imodSize = modelSize(imodNum,1)
         kmodSize = modelSize(imodNum,2)

C       chemkin call uses species from reduced-order ODE mapped into full concentration array
C
         call ackwc(imodSize,p,t,cfull,nur(1,imodNum),nus(1,imodNum),
     & ickwrk,rckwrk,ispecies,wdot)
C     return the reduced model species binary array to CFD code
         do k=1,kk
            ispecies(k) = nus(k,imodNum)
         enddo
      elseif(jsolv .gt. 1) then
         print*,'Error: Invalid jsolv value:',jsolv
         print*,'jsolv must be -1, 0 or 1'
         stop
      endif
C     </Calc. rates>

      return
      end

C*************************************************************************
      subroutine needmodel(igrid,igridTot,igrTotm1,nsteps,icount,
     & icountm1,jfindmod,iwritesize,jprint,nprint,iicount)
      implicit double precision (a-h, o-z), integer (i-n)
      dimension iicount(*)
C     Keep track of age of current model distribution
C     icount changes on a new iteration, remains constant for all calls
C     on the same iteration (e.g calls for jacobian calculations)
C     Force update if number of gridpoints has changed
      if(igridTot.ne.igrTotm1) then
         do j = 1, igridTot
            iicount(j) = -nsteps
         enddo
C     <Update>
         igrTotm1 = igridTot
C     </Update>
      endif
C     default: don't need to update models
      jfindmod = 0
C     default: don't write model use stats to file
      iwritesize = 0
C     iicount(igrid) = iteration number of most recent model update for igrid
      iicountI = iicount(igrid)
      if(iicountI.ne.icount) then
C     new iteration/integration step at this grid point;
         if(icount.eq.0) then
C     user wishes to force model update
            istepI = 0
         else
C     determine whether update is needed
            istepI = mod(iicountI,nsteps)
         endif
         if(istepI.eq.0) then
C     Models at retirement age, or first call, or new mesh
            jfindmod = 1
         endif
C     update iicount(igrid)
         iicount(igrid) = icount
      endif
C     <modelUsage statistics>
      if(icount.ne.icountm1) then
C     new step
C     increment count of number of steps between 'modelUsage' writes
         jprint = jprint + 1
         if(jprint.eq.nprint) then
            jprint = 0
            iwritesize = 1
         endif
      endif
C     </modelUsage>

      return
      end

C************************************************************************
      subroutine adapchinit(igridTot,icount,igrTotm1,icountm1,ii,kk,
     & kfuel,resid_tol,jprint,nprint,nsteps,iicount,iModel,nFull,nRed,
     & rmaxfull)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension iModel(*), iicount(*), iCall(1)
      character line*5

C       initialize and set defaults on number of steps between modelUsage write
      jprint = 0
      nprint = 1
C     for calculating average model size; printed to avgModelSize
      icountTot = 0
      iSizeTot = 0
C     modelUsage:
      nFull = 0
      nRed = 0
      rmaxfull = 1.
      igrTotm1 = igridTot
      icountm1 = icount
      do j=1, igridTot
         iModel(j) = 0
         iicount(j) = 0
      end do
      print*,'reading adapchem.in...'
      open(7, file='adapchem.in', status='old')
      do while(7.gt.0)
         read(7,*) line
         if(line(:3).eq.'END') goto 005
         if(line.eq.'TOLSS') read(7,*) tolss
         if(line.eq.'TOLMR') read(7,*) tolmr
         if(line.eq.'NSTEP') read(7,*) nsteps
         if(line.eq.'NPRIN') read(7,*) nprint
C     <modelLib>
         if(line.eq.'MAXFU') read(7,*) rmaxfull
C     </modelLib>
      enddo
 005  continue
      close(7)
      resid_tol = tolss - tolmr
      print*,'resid_tol= ',resid_tol
      print*,'nsteps= ',nsteps
      print*,'rmaxfull= ',rmaxfull
      print*,'reading models.chem...'

      open(1, file='models.chem', status='old')
      read (1,*) ii, kk!Luwi-9/2/10-kfuel not needed/used!, kfuel
C     closed in subroutine findModel
!       <Luwi--9/3/2010-- avoiding multiple file copies in parallel code>
CC     Write headers for model usage documentation
C     open(unit=2,file='avgModelSize',status='unknown')
C     write(2,'(A14)') '======================'
C     write(2,'(A4,X,A9)') 'nRxns','% Of Full'
C     write(2,'(A14)') '======================'
CC
C     open(unit=3,file='modelUsage',status='unknown')
C     write(3,'(A22)')        '======================'
C     write(3,'(A10,2X,A10)') '     nFull','  nReduced'
C     write(3,'(A22)')        '======================'
!       <Luwi>
      return
      end

C**********************************************************************
      subroutine sizewrite(nFull,nRed,iSizeTot,icountTot,ii,kk)
      implicit none
      integer iSizeTot, icountTot, ii, nFull, nRed,kk
      write(2,1) nint(dble(iSizeTot)/dble(icountTot)),
     &dble(iSizeTot)/dble(icountTot*kk)*100.
      write(3,2) nFull, nRed

    1 format(I4,2X,F5.1)
    2 format(I10,2X,I10)

      return
      end

C**********************************************************************
      subroutine readModels(lin,kk,ii,rangeLo,rangeHi,nModels)
C     O.Oluwole 6/9/2004
C     This subroutine reads in the reduced models in the model library
C     for use during the adaptive chemistry calculation
      implicit double precision (a-h,o-z), integer (i-n)
      include 'ackparams.i'
      dimension rangeLo(kmax+1,nm+1), rangeHi(kmax+1,nm+1)
      dimension modelSize(nm+1,2)
      dimension NUR(imax,nm+1), NUS(kmax,nm+1)
      common /saveNUS/ NUS
      common /saveNUR/ NUR, modelSize
      character line*12
C     Columns of NUR(,) contain all models including full model
C     Columns of rangeLo and rangeHi contain valid ranges of reduced
C     models beginning in column 2
C     all reduced models are stored in order of increasing size
C     <O.Oluwole--6/26/2008>
C     Aerodyne Research, Inc.
C     For NASA Phase I SBIR with MIT
C
C     First store default model to be used whenever no valid reduced model
C     in library--usually full model, but another (smaller) model may be
C     specified (a subset of the full model). This is to be used for cases
C     where user wishes to specify a maximum model size to be used over the
C     course of the Adaptive Chemistry simulation. Adaptive Chemistry should
C     evaluate and report the error of the skeletal model at the end of the
C     steady-state simulation
      open(4,file='default.model',status='old')
      read(4,*) ivalr, ivals
      if(ivalr.eq.ii) then
C     default is full model
         do i=1,ii
            NUR(i,1) = 1
         enddo
         do k=1,kk
            NUS(k,1) = 1
         enddo
      else
C     default is subset of full model
         read(4,*) ( NUR(i,1), i=1, ii )
         read(4,*) ( NUS(k,1), k=1, kk )
      endif
      close(4)
      modelSize(1,1) = ivalr
      modelSize(1,2) = ivals ! number of species in default model
      nModels = 1
      ireadSize = 1
      ifoundSize = 0
      imodelSize = ii*10
      kmodelSize = kk*10
      do 100 while(7.gt.0)
         read(lin,*) line
         if(line(:3).eq.'END') return
         if(ireadSize.eq.1) then
            if (ifoundSize.eq.0) then
C     look for model size
              if(line(2:9).eq.'Contains'.or.line(2:9).eq.'Includes')then
                backspace(lin)
                read(lin,*) line, kmodelSize
                ifoundSize = 1
              endif
            elseif(ifoundSize.eq.1) then !found number of reactions; need number of species
C     look for model species size
          if(line(2:9).eq.'Contains'.or.line(2:9).eq.'Includes') then
              backspace(lin)
              read(lin,*) line, imodelSize
              ifoundSize = 2
          elseif(line(:6).eq.'REGION') then
C     model size not specified; never look for model size
              ireadSize = 0
          endif
            endif
         endif
         if(line(:6).eq.'REGION') then
            nModels = nModels + 1
C     Store the model
            read(lin,*) (NUS(k,nModels), k=1,kk)!species
            read(lin,*) (NUR(i,nModels), i=1,ii)!reactions
            modelSize(nModels,1) = imodelSize
            modelSize(nModels,2) = kmodelSize
C     Valid range
            read(lin,*) line
            do k=1, kk+1
         read(lin,*) rangeLo(k,nModels), rangeHi(k,nModels)
C     <approx. for zero>
         if(rangeLo(k,nModels).lt.1.d-15) rangeLo(k,nModels) = 0.0d0
         rangeHi(k,nModels) = max(rangeHi(k,nModels),1.0d-15)
C     </approx. for zero>
            enddo
            ifoundSize = 0
         endif
  100 continue

      return
      end

C**********************************************************************
      subroutine findModel(jsolv,lin,iCall,kk,ii,kfuel,T,C,iModNum,
     & nFull,nRed)
C     This subroutine finds the smallest valid model within the model library
C     at a given set of species concentrations and temperature
C     Luwi Oluwole 6/10/2004
      implicit double precision (a-h,o-z), integer (i-n)
      include 'ackparams.i'
      dimension C(kk)
      dimension rangeLo(kmax+1,nm+1), rangeHi(kmax+1,nm+1)
C     save data for use on subsequent calls
      common /saveV/ rangeLo,rangeHi
      common /saveV1/ nModels
      if(iCall.eq.2) goto 007
C     On first call, read and save model library
      call readModels(lin,kk,ii,rangeLo,rangeHi,nModels)
      close(lin)
      iCall = 2
C     <modelLib: signal completion of model library upload>
      open(78,file='donemodels',status='new')
      close(78)
      if(jsolv .eq. -1) then
        return ! Luwi--4/23/2010--this call was for loading models only
        endif
C     </modelLib>
  007 continue
C
C     If no valid reduced model is found use full model
      iModNum = 1
C     Search model library for smallest valid model
      do 001 j = 2, nModels
         if(T.lt.rangeLo(1,j).or.T.gt.rangeHi(1,j)) goto 001
         do k = 1, kk
            if(abs(C(k)).lt.rangeLo(k+1,j).or.C(k).gt.rangeHi(k+1,j))
     &      goto 001
         enddo
         iModNum = j
         nRed = nRed + 1
         return
  001 continue
C     Increment count of number of times full model is used
      nFull = nFull + 1

      return
      end

C**********************************************************************
      subroutine closeadcfiles
      implicit double precision (a-h, o-z)
C     close model Usage statistics files
      close (2) !filename: avgModelSize
      close (3) !filename: modelUsage

      return
      end

C**********************************************************************
      SUBROUTINE ACKWC(nrxn,p,t,c,nur,nus,ickwrk,rckwrk,ispecies,wdotr)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C     Returns the molar production rates of the species given the
C     pressure, temperature and mass fractions;  see Eq. (49).
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), wdotr(*), c(*), nur(*)
      dimension ispecies(*)
      dimension nus(*), wdot(nkk)!reduced species
C
      INCLUDE 'ckstrt.i'
C
      DO 50 K = 1, NKK
       WDOT(K) = 0.0
 50   CONTINUE
C     no further operations necessary for 0-reaction model
      if(nrxn.eq.0) return

      CALL ACKRATT (nur,RCKWRK, ICKWRK, NII, MXSP, RCKWRK(NcRU),
     1             RCKWRK(NcPA), T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO), NREV,
     3             ICKWRK(IcRV), RCKWRK(NcRV), NLAN, NLAR, ICKWRK(IcLT),
     4             RCKWRK(NcLT), NRLT, ICKWRK(IcRL), RCKWRK(NcRL),
     5             RCKWRK(NcK1), RCKWRK(NcKF), RCKWRK(NcKR),
     6             RCKWRK(NcI1), NRNU, ICKWRK(IcRNU), RCKWRK(NcRNU))
C
C
      CALL ACKRATX (p,rckwrk(NcRU),nur, NII, NKK, MXSP, MXTB, T, C,
     $             ICKWRK(IcNS),
     1             ICKWRK(IcNU), ICKWRK(IcNK), NPAR+1, RCKWRK(NcCO),
     2             NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF), NFAR,
     3             RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),
     4             RCKWRK(NcKT), ICKWRK(IcKT), RCKWRK(NcKF),
     5             RCKWRK(NcKR), RCKWRK(NcI1), RCKWRK(NcI2),
     6             RCKWRK(NcI3), NRNU, ICKWRK(IcRNU), RCKWRK(NcRNU),
     7             NORD, ICKWRK(IcORD), MXORD, ICKWRK(IcKOR),
     8             RCKWRK(NcKOR),
     $             ispecies)
C
      DO 101 N = 1, MXSP
         DO 100 I = 1, NII
      if(nur(i).eq.0) goto 100
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            ROP = RCKWRK(NcI1+I-1) - RCKWRK(NcI2+I-1)
            IF (K .NE. 0) THEN
               RNU =
C*****precision > double
     1         DBLE(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
C*****END precision > double
C*****precision > single
C     1         REAL (ICKWRK(IcNU + (I-1)*MXSP + N - 1))
C*****END precision > single
               WDOT(K) = WDOT(K) + RNU * ROP
            ENDIF
  100 CONTINUE
  101 continue
C
      IF (NRNU .LE. 0) goto 007
      DO 201 L = 1, NRNU
         I = ICKWRK(IcRNU + L - 1)
      if(nur(i).eq.0) goto 201
         DO 200 N = 1, MXSP
            K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
            ROP = RCKWRK(NcI1+I-1) - RCKWRK(NcI2+I-1)
            IF (K .NE. 0) THEN
               RNU = RCKWRK(NcRNU + (L-1)*MXSP + N - 1)
               WDOT(K) = WDOT(K) + RNU * ROP
            ENDIF
  200 CONTINUE
  201 continue

C     <Luwi Oluwole: 3/22/2010--reduced species wdot>
  007 kr = 0
      do k=1,nkk
         if(nus(k).eq.1) then
            kr = kr+1
            wdotr(kr) = wdot(k)
         endif
      enddo
C     </Luwi>

      RETURN
      END

C**********************************************************************
      SUBROUTINE ACKRATT(nur,RCKWRK, ICKWRK, II, MAXSP, RU, PATM, T,
     1              NSPEC, NU, NUNK, NPAR, PAR, NREV, IREV, RPAR, NLAN,
     2              NLAR, ILAN, PLT, NRLT, IRLT, RPLT, SMH, RKFT,
     3              RKRT, EQK, NRNU, IRNU, RNU)
C
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RCKWRK(*), ICKWRK(*), NSPEC(*), NU(MAXSP,*),
     1          NUNK(MAXSP,*), PAR(NPAR,*), IREV(*), RPAR(NPAR,*),
     2          ILAN(*), IRLT(*), PLT(NLAR,*), RPLT(NLAR,*), SMH(*),
     3          RKFT(*), RKRT(*), EQK(*), IRNU(*), RNU(MAXSP,*)
      dimension nur(*)
      COMMON /MACH/ SMALL,BIG,EXPARG
      ALOGT = LOG(T)
      DO 20 I = 1, II
       if(nur(i).eq.0) goto 20
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)/T)
   20 CONTINUE
C
C        Landau-Teller reactions
C
      DO 25 N = 1, NLAN
         I = ILAN(N)
       if(nur(i).eq.0) goto 25
         TFAC = PLT(1,N)/T**(1.0/3.0) + PLT(2,N)/T**(2.0/3.0)
         RKFT(I) = RKFT(I) * EXP(TFAC)
   25 CONTINUE
C
      CALL ACKSMH (T, ICKWRK, RCKWRK, SMH)
      DO 50 I = 1, II
       if(nur(i).eq.0) goto 50
          SUMSMH = 0.0
          DO 40 N = 1, MAXSP
             IF (NUNK(N,I).NE.0)
C*****precision > double
     1         SUMSMH = SUMSMH + DBLE(NU(N,I))*SMH(NUNK(N,I))
C*****END precision > double
C*****precision > single
C     1           SUMSMH = SUMSMH + REAL(NU(N,I))*SMH(NUNK(N,I))
C*****END precision > single
   40     CONTINUE
          IF (SUMSMH .NE. 0.0) EQK(I) = EXP(MIN(SUMSMH,EXPARG))
   50 CONTINUE
C
      DO 55 N = 1, NRNU
         SUMSMH = 0.0
         I = IRNU(N)
       if(nur(i).eq.0) goto 55
         DO 45 L = 1, MAXSP
            IF (NUNK(L,I).NE.0) SUMSMH=SUMSMH+RNU(L,N)*SMH(NUNK(L,I))
   45    CONTINUE
         IF (SUMSMH .NE. 0.0) EQK(I) = EXP(MIN(SUMSMH,EXPARG))
   55 CONTINUE
C
      PFAC = PATM / (RU*T)
      DO 60 I = 1, II
       if(nur(i).eq.0) goto 60
         NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
         EQK(I) = EQK(I) * PFAC**NUSUMK
   60 CONTINUE
      DO 65 N = 1, NRNU
         I = IRNU(N)
       if(nur(i).eq.0) goto 65
         RNUSUM = RNU(1,N)+RNU(2,N)+RNU(3,N)+RNU(4,N)+RNU(5,N)+RNU(6,N)
         PFR = PFAC ** RNUSUM
         EQK(I) = EQK(I) * PFR
   65 CONTINUE
C
      DO 68 I = 1, II
C
C     RKR=0.0 for irreversible reactions, else RKR=RKF/MAX(EQK,SMALL)
C
         RKRT(I) = 0.0
       if(nur(i).eq.0) goto 68
         IF (NSPEC(I).GT.0) RKRT(I) = RKFT(I) / MAX(EQK(I),SMALL)
   68 CONTINUE
C
C     if reverse parameters have been given:
C
      DO 70 N = 1, NREV
         I = IREV(N)
       if(nur(i).eq.0) goto 70
         RKRT(I) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT - RPAR(3,N)/T)
         EQK(I)  = RKFT(I)/RKRT(I)
   70 CONTINUE
C
C     if reverse Landau-Teller parameters have been given:
C
      DO 75 N = 1, NRLT
         I = IRLT(N)
       if(nur(i).eq.0) goto 75
         TFAC = RPLT(1,N)/T**(1.0/3.0) + RPLT(2,N)/T**(2.0/3.0)
         RKRT(I) = RKRT(I) * EXP(TFAC)
         EQK(I) = RKFT(I)/RKRT(I)
   75 CONTINUE
C
      RETURN
      END

C**********************************************************************
      SUBROUTINE ACKRATX
     &  (p,ru,nur,II, KK, MAXSP, MAXTB, T, C, NSPEC, NU,
     1     NUNK, NPAR, PAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR,
     2     NTHB, ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF,
     3     RKR, CTB, NRNU, IRNU, RNU, NORD, IORD, MXORD,
     4     KORD, RORD, ispecies)
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION C(*), NSPEC(*), NU(MAXSP,*), NUNK(MAXSP,*), PAR(NPAR,*),
     1          IFAL(*), IFOP(*), KFAL(*), FPAR(NFAR,*), ITHB(*),
     2          NTBS(*), AIK(MAXTB,*), NKTB(MAXTB,*), RKFT(*),
     3          RKRT(*), RKF(*), RKR(*), CTB(*), IRNU(*), RNU(MAXSP,*),
     4          IORD(*), KORD(MXORD,*), RORD(MXORD,*)
C
      dimension nur(*)
      dimension ispecies(*)
      COMMON /MACH/ SMALL,BIG,EXPARG
C
      DO 20 I = 1, II
       if(nur(i).eq.0) goto 20
         CTB(I) = 1.0
         RKF(I) = 0.0
         RKR(I) = 0.0
   20 CONTINUE
C
C     third-body reactions
C
      IF (NTHB .GT. 0) THEN
         CTOT = 0.0
         DO 10 K = 1, KK
            CTOT = CTOT + C(K)
   10    CONTINUE
c         ctot = p/(ru*t)
         DO 81 N = 1, NTHB
          i = ithb(n)
          if(nur(i).eq.0) goto 81
            CTB(ITHB(N)) = CTOT
            DO 80 L = 1, NTBS(N)
               CTB(ITHB(N)) = CTB(ITHB(N)) + (AIK(L,N)-1.0)*C(NKTB(L,N))
             ispecies(nktb(l,n)) = 1
   80       CONTINUE
   81    continue
      ENDIF
C
C     If fall-off (pressure correction):
C
      IF (NFAL .GT. 0) THEN
         ALOGT = LOG(T)
C
         DO 90 N = 1, NFAL
          i = ifal(n)
          if(nur(i).eq.0) goto 90
C
            RKLOW = FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)/T)
C
C        CONCENTRATION OF THIRD BODY
C
            IF (KFAL(N) .EQ. 0) THEN
               PR = RKLOW * CTB(IFAL(N)) / RKFT(IFAL(N))
               CTB(IFAL(N)) = 1.0
            ELSE
               PR = RKLOW * C(KFAL(N)) / RKFT(IFAL(N))
             ispecies(kfal(n)) = 1
            ENDIF
C
            PCOR = PR / (1.0 + PR)
C
            IF (IFOP(N) .GT. 1) THEN
               PRLOG = LOG10(MAX(PR,SMALL))
C
               IF (IFOP(N) .EQ. 2) THEN
C
C              8-PARAMETER SRI FORM
C
                  XP = 1.0/(1.0 + PRLOG**2)
                  FC = ((FPAR(4,N)*EXP(-FPAR(5,N)/T)
     1                   + EXP(-T/FPAR(6,N))) **XP)
     2                  * FPAR(7,N) * T**FPAR(8,N)
C
               ELSE
C
C              6-PARAMETER TROE FORM
C
                  FCENT = (1.0-FPAR(4,N)) * EXP(-T/FPAR(5,N))
     1                  +       FPAR(4,N) * EXP(-T/FPAR(6,N))
C
C              7-PARAMETER TROE FORM
C
                  IF (IFOP(N) .EQ. 4) FCENT = FCENT + EXP(-FPAR(7,N)/T)
C
                  FCLOG = LOG10(MAX(FCENT,SMALL))
                  XN    = 0.75 - 1.27*FCLOG
                  CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
                  FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
                  FC = 10.0**FLOG
               ENDIF
               PCOR = FC * PCOR
            ENDIF
C
            RKFT(IFAL(N)) = RKFT(IFAL(N)) * PCOR
            RKRT(IFAL(N)) = RKRT(IFAL(N)) * PCOR
   90    CONTINUE
      ENDIF
C
C     Multiply by the product of reactants and product of products
C     PAR(4,I) is a perturbation factor
C
      DO 150 I = 1, II
       if(nur(i).eq.0) goto 150
         RKFT(I) = RKFT(I) * CTB(I) * PAR(4,I)
         RKRT(I) = RKRT(I) * CTB(I) * PAR(4,I)
C
         IF (NU(1,I) .NE. 0) THEN
            RKF(I) = RKFT(I)*C(NUNK(1,I))**IABS(NU(1,I))
          ispecies(nunk(1,i)) = 1

            RKR(I) = RKRT(I)*C(NUNK(4,I))**NU(4,I)
          ispecies(nunk(4,i)) = 1

            IF (NUNK(2,I) .NE. 0) THEN
               RKF(I)= RKF(I) * C(NUNK(2,I))**IABS(NU(2,I))
             ispecies(nunk(2,i)) = 1

               IF (NUNK(3,I) .NE. 0) then
          RKF(I) = RKF(I) * C(NUNK(3,I))**IABS(NU(3,I))
          ispecies(nunk(3,i)) = 1
             endif

            ENDIF
            IF (NUNK(5,I) .NE. 0) THEN
               RKR(I) = RKR(I) * C(NUNK(5,I))**NU(5,I)
             ispecies(nunk(5,i)) = 1

               IF (NUNK(6,I) .NE. 0) then
          RKR(I)=RKR(I)*C(NUNK(6,I))**NU(6,I)
          ispecies(nunk(6,i)) = 1
             endif
            ENDIF
         ENDIF
  150 CONTINUE
C
      DO 160 N = 1, NRNU
         I = IRNU(N)

       if(nur(i).eq.0) goto 160

         C1 = C(NUNK(1,I)) ** ABS(RNU(1,N))
       ispecies(nunk(1,i)) = 1

         C4 = C(NUNK(4,I)) ** RNU(4,N)
       ispecies(nunk(4,i)) = 1

         RKF(I) = RKFT(I) * C1
         RKR(I) = RKRT(I) * C4
         IF (NUNK(2,I) .NE. 0) THEN
            C2 = C(NUNK(2,I)) ** ABS(RNU(2,N))
          ispecies(nunk(2,i)) = 1

            RKF(I) = RKF(I) * C2
            IF (NUNK(3,I) .NE. 0) THEN
               C3 = C(NUNK(3,I)) ** ABS(RNU(3,N))
             ispecies(nunk(3,i)) = 1

               RKF(I) = RKF(I) * C3
            ENDIF
         ENDIF
         IF (NUNK(5,I) .NE. 0) THEN
            C5 = C(NUNK(5,I)) ** RNU(5,N)
          ispecies(nunk(5,i)) = 1

            RKR(I) = RKR(I) * C5
            IF (NUNK(6,I) .NE. 0) THEN
               C6 = C(NUNK(6,I)) ** RNU(6,N)
             ispecies(nunk(6,i)) = 1

               RKR(I) = RKR(I) * C6
            ENDIF
         ENDIF
  160 CONTINUE
C
      DO 200 N = 1, NORD
         I = IORD(N)

       if(nur(i).eq.0) goto 200

         RKF(I) = RKFT(I)
         RKR(I) = RKRT(I)
C
         DO 190 L = 1, MXORD
            NK = KORD(L,N)
            IF (NK .LT. 0) THEN
               NK = IABS(NK)
               CNK = C(NK) ** RORD(L,N)
             ispecies(nk) = 1

               RKF(I) = RKF(I) * CNK
            ELSEIF (NK .GT. 0) THEN
               CNK = C(NK) ** RORD(L,N)
             ispecies(nk) = 1

               RKR(I) = RKR(I) * CNK
            ENDIF
  190    CONTINUE
  200 CONTINUE
C
      RETURN
      END

C**********************************************************************
      SUBROUTINE ACKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
C     Returns the molar concentrations given the pressure,
C     temperature and mass fractions;  see Eq. (7).
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     C      - Molar concentrations of the species.
C                   cgs units - mole/cm**3
C                   Data type - real array
C                   Dimension C(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), C(*)
C
      INCLUDE 'ckstrt.i'
C
      SUMYOW = 0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K) / RCKWRK(NcWT + K - 1)
150   CONTINUE
      SUMYOW = SUMYOW * T * RCKWRK(NcRU)
      DO 200 K = 1, NKK
         C(K) = P * Y(K) / (SUMYOW * RCKWRK(NcWT + K - 1))
200   CONTINUE
      RETURN
      END

C**********************************************************************
      SUBROUTINE ACKSMH  (T, ICKWRK, RCKWRK, SMH)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)*
C     Returns the array of entropies minus enthalpies for the species.
C     It is normally not called directly by the user.
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     SMH    - Entropy minus enthalpy for the species,
C              SMH(K) = S(K)/R - H(K)/RT.
C                   cgs units - none
C                   Data type - real array
C                   Dimension SMH(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(*), RCKWRK(*), SMH(*), TN(10)
C
      INCLUDE 'ckstrt.i'
C
      TN(1) = LOG(T) - 1.0
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)/((N-1)*N)
 150  CONTINUE
C
      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         SUM = 0.0
         DO 225 N = 1, NCP
            SUM = SUM + TN(N)*RCKWRK(NA1 + N - 1)
  225    CONTINUE
         SMH(K) = SUM + RCKWRK(NA1 + NCP2 - 1)
     1                - RCKWRK(NA1 + NCP1 - 1)/T
C
  250 CONTINUE

      RETURN
      END

C**********************************************************************
C Initialize variables
      subroutine adapchemINI()
      dimension icall(1)
      common/adpchm7/ inewlib, iCall
        common /adpchm8/iadapnew
        iadapnew = 1
!       <Luwi--9/2/2010-->
      icall(1) = 0
!       </Luwi>
      end


C**********************************************************************
C create an empty file to signal terminating model library generation script
      subroutine adapchemEND()
         open(unit=78, file='full.rstrt',status='new')
           close(78)
      end
