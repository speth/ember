        subroutine adapbox(igridtot,nmaxbox,y,p,temp,rwork,iwork)
C
C	Divide the computational domain into sub-domains based on the similarity
C       of state variables
C       To generate input files for RIOT
C       Author: Yu Shi
C       Date: 06/24/2010
C       Input file: boxDomain.in
C       Inputs: igridtot,nmaxbox,y(), p(), temp(),rwork,iwork
C       	igridtot: The total number of CFD grids
C       	nmaxbox: The maximum number of sub-domains to be divided
C       	y(igridtot,NKK): mass fraction of species for each CFD grid
C       	p(igridtot): pressure of each CFD grid, unit in dyne
C       	temp(igridtot): temperature of each CFD grid, unit in K
C       	rwork(): Chemkin parameter real array
C       	iwrok(): Chemkin parameter integer array
C               common block /box1/
C		commom block /adpchm4/
C		common block /adpchm7/
C       Outputs: Files: save#.dat, numpts#, statmodLib, newmodLib
C               commom block /adpchm7/
C
	implicit none
        include 'ckstrt.i'

	common /adpchm4/rmaxfull
	common /adpchm5/nFull,nRed
        common/adpchm7/ inewlib, iCall
        common/adpchm8/ iadapnew
	common /counter/icount
        integer icount(1)
        integer inewlib,iCall(1)
        integer iadapnew
	integer nFull, nRed
	double precision rmaxfull
	logical flagf1, flagf2

        integer igridtot,nmaxbox, i, j, k, jj,nsub,ntot
        integer iwork(*)

        double precision y(igridtot,NKK),p(igridtot),temp(igridtot)
        double precision rwork(*)


        integer isort(igridtot)
        integer isave, igridsub, igridsub1, igridsub2
        integer idxsavedom(nmaxbox,igridtot+1),idxtempdom(igridtot)
        integer numpts, kerr,istrlen

        character*16 string, filename1, filename2

        double precision tmax,tmin,deltaT(nmaxbox)
        double precision tratemax, tratemin, tratediff, tup,tlo
        double precision ytmp(NKK),tmpwdot(NKK),wdot(igridtot,NKK),ymin
        double precision temp_rate(igridtot), sum, cpb, rho
        double precision SYS_MIN,SYS_MAX
        SYS_MIN = 1.d-100
        SYS_MAX = 1.d100
        ymin = 1.e-15

C To detect if new model library is requested
        icount(1) = icount(1) + 1
	if ((nRed+nFull).ne.0 .and. mod(icount(1),100).eq.0) then
	   write(*,*)'Adap: Red ', nRed, 'Full ', nFull
     &, 'percent:', nRed*1.0/(nFull+nRed)*100.0
	endif
	if(inewlib.eq.0 .and. dble(nFull)/dble(nRed+nFull)
     &.gt. rmaxfull) inewlib=1

        if(inewlib.eq.2) then
c     if "in process", check status of new model library; based on feedback
c     from model library tool
        	open(76,file='statmodLib',status='old')
             	read(76,*) inewlib
             	close(76)
c     if "done", signal adapchem of completion of new library
             	if(inewlib.eq.3) then
			iCall(1) = 0
			iadapnew = 1
                else
                        iadapnew = 0
		endif
		return
        endif

        if(inewlib.eq.3) then
c     after adapchem notified and one macro step taken by solver
c     can assume that new library has been loaded,
c     so reset flag
                inewlib = 0
                iadapnew = 0
                close(76)
		return
	endif

C If these two files exist, no new boxes will be generated as
C they indicate the reduced models are still being created

	INQUIRE(file='save1.dat',EXIST=flagf1)
      	INQUIRE(file='makingLib',EXIST=flagf2)
      	if(inewlib.ne.1 .or. flagf1 .or. flagf2) return

C If all above conditions are passed, this indicates new boxes need to be created
C Calculate species rate and temperature rate
	do i = 1, igridtot
            ytmp(:)=y(i,:)
	    call CKWYP (p(i), temp(i),ytmp, iwork, rwork,tmpwdot)
            call CKRHOY (p(i), temp(i), ytmp, iwork, rwork, rho)
            wdot(i,:) = tmpwdot(:)
            call ckcvbs (temp(i),ytmp, iwork, rwork, cpb)
            call ckums (temp(i), iwork, rwork, rwork(NcK1))
            sum = 0.0
            do j = 1, NKK
                sum = sum + rwork(NcK1+j-1)*tmpwdot(j)*rwork(NcWT+j-1)
            enddo
            temp_rate(i)=-sum/(rho*cpb)
	enddo


C Sort grids in tempeature ascending order, index are stored in isort

        Call SORTRX(igridtot,temp,isort)
        tmax = temp(isort(igridtot))
        tmin = temp(isort(1))

C Determin deltaT in order to divide grids into sub-domains based on temperature only
C Each temperature sub-domain allows for 2 sub-domains based on temperature rate
        nsub = int(nmaxbox/2.0)
        if (nsub .gt. 1) then
          if (tmin .lt. 400.0) then
                deltaT(1) = tmin            !The first zone from tmin to 400 K
                do i = 2, nsub
                  deltaT(i) = 400.0+(i-2)*(tmax-400.0)/(nsub-1)
                enddo
          else
                do i = 1, nsub
                  deltaT(i) = tmin+(i-1)*(tmax-tmin)/nsub
                enddo
          endif
        else
                do i = 1, nsub
                  deltaT(i) = tmin+(i-1)*(tmax-tmin)/nsub
                enddo
        endif
        deltaT(nsub+1)=tmax !The last upper bound in temperature
C Create groups
        isave = 0
	do i =1, nsub
	    igridsub=0
            tratemax = SYS_MIN
            tratemin = SYS_MAX
            do j = 1, igridtot
                k = isort(j)
                tup = deltaT(i+1)
                tlo = deltaT(i)
                if ((temp(k).ge.tlo .and. temp(k).lt.tup)
     & .or. (i.eq.nsub.and.j.eq.igridtot)) then  !special treatment for the last cell
                    igridsub = igridsub+1
                    idxtempdom(igridsub) = k
                    tratemax = max(tratemax,temp_rate(k))
                    tratemin = min(tratemin,temp_rate(k))
                endif
            enddo
            if (igridsub .le. 2) then    !if only find 2 cells or less, they are grouped without considering their temperature rates
                isave = isave + 1
                idxsavedom(isave,1)=igridsub
                do j = 1, igridsub
                    idxsavedom(isave,j+1)=idxtempdom(j)
                enddo
            else                        !subdivide the temperature zone into 2 zones based on temperature rate
                igridsub1 = 0
                igridsub2 = 0
                do j = 1, igridsub
                    k = idxtempdom(j)
                    tratediff = tratemax-tratemin
                    if (temp_rate(k).lt.(tratemin+tratediff/2.0)) then !first sub-zone
                        igridsub1 = igridsub1+1
                    endif

                    if (temp_rate(k).ge.(tratemin+tratediff/2.0)) then
                        ! secondsub-zone
                        igridsub2 = igridsub2+1
                    endif
                enddo

                if (igridsub1.ne.0) then
                    isave = isave + 1
                    idxsavedom(isave,1)=igridsub1
                    jj = 1
                    do j = 1, igridsub
                        k = idxtempdom(j)
                    if (temp_rate(k).lt.(tratemin+tratediff/2.0)) then
                        jj = jj + 1
                        idxsavedom(isave,jj)=k
                    endif
                    enddo
                endif
                if (igridsub2.ne.0 ) then
                    ! if both subdomains contain more than 10 cells, do not merge them
                    if (igridsub1 .gt. 10 .and. igridsub2 .gt. 10)then
                        isave = isave + 1
                        idxsavedom(isave,1)=igridsub2
                        jj = 1
                        do j = 1, igridsub
                            k = idxtempdom(j)
                            if (temp_rate(k).ge.
     &(tratemin+tratediff/2.0)) then
                                jj = jj + 1
                                idxsavedom(isave,jj)=k
                             endif
                        enddo
                    else  !merge the subdomains
                        if(igridsub1 .eq. 0) then   !if the first subdomain is not created
                                isave = isave + 1
                                idxsavedom(isave,1)= igridsub2
                        else
                                idxsavedom(isave,1)= igridsub !merge to the first subdomain
                        endif
                        jj = igridsub1+1
                        do j = 1, igridsub
                            k = idxtempdom(j)
                            if (temp_rate(k).ge.
     &(tratemin+tratediff/2.0)) then
                                jj = jj + 1
                                idxsavedom(isave,jj)=k
                            endif
                        enddo
                    endif
                endif
           endif
       enddo
C scrutiny check
      ntot = 0
      do i = 1, isave
         ntot = ntot + idxsavedom(i,1)
      enddo
      if (ntot.ne.igridtot)then
          write(*,*)'Total number of grids in all boxes is not equal
     &to igridtot, program exits in box.f!'
          write(*,*)'ntot:',ntot,'igridtot:',igridtot
          call exit(0)
      endif

      open(1,file='nbox',status='new')
      write(1,'(I16)')isave
      close(1)

      do i = 1, isave
          call cki2ch(i,string,istrlen,kerr)
          filename1 = 'numpts'//string(1:istrlen)
          filename2 = 'save'//string(1:istrlen)//'.dat'
          open(2,file=filename1,status='new')
          numpts = idxsavedom(i,1)
          write(2,'(A6,I16)') 'numpts',numpts
          close(2)
          open(3,file=filename2,status='new')
          do j = 1, numpts
             k=idxsavedom(i,j+1)
             write(3,*)p(k)
             ytmp(:)=y(k,:)
             write(3,*)temp(k),(max(ytmp(jj),ymin),jj=1,NKK)
          enddo
          close(3)
      enddo

c     and signal request to model library tool
      open(77,file='newmodLib',status='new')
      close(77)
c     update status of model library to "in process"
      inewlib = 2
c     open status file
      open(76,file='statmodLib',status='unknown')
      write(76,*) inewlib
      close(76)
      return
      end subroutine
