c 
c.......................................................................
c
      subroutine iqaom (aom,coef,ncent,nbasins,nmo,nprims,iover,lw,lr,
     &  lu18,wfnfile,qtaim,mulliken,lowdin,mindef,becke,mindefrho,
     &  netrho,promrho,hesel,largwr)
      implicit none
      real    (kind=8) coef(nmo+nmo,nprims)
      real    (kind=8) aom(ncent,nmo,nmo)
      integer (kind=4) iaom(ncent)
      real    (kind=8) rhopow,alphahesel,betahesel
      integer (kind=4) i,j,k,l,m,n,lr,lw,irho
      integer(kind=4)  iover,lu18,ii,nbasins,natoms,ncent,nmo,nprims
      real(kind=8),    allocatable,dimension (:,:,:) :: aomtopmod
      real(kind=8),    allocatable,dimension (:,:)   :: cc
      real(kind=8),    allocatable,dimension (:,:)   :: eta
      logical          mulliken,lowdin,qtaim,mindef,becke
      logical          mindefrho,netrho,promrho,hesel,notaom
      logical          largwr,okaom,warn,warno
      character*16     namebas(100)
      real(kind=8),    parameter ::   deftolaom     =  0.05D0
      real(kind=8),    parameter ::   defcritics    =  0.02D0
      real(kind=8),    parameter ::   defahesel     =  2.00D0
      real(kind=8),    parameter ::   defbhesel     =  2.00d0
      real(kind=8),    parameter ::   defdamps      =  0.95d0
      real(kind=8),    parameter ::   defcovx       =  1.20d0
      real(kind=8),    parameter ::   ranmindef     = -5.00d0
      real(kind=8),    parameter ::   ranmaxdef     = +5.00d0
      real(kind=8),    parameter ::   defepsproba   =  1.00D-14
      real(kind=8),    parameter ::   defepsneg     = -1.00D-6
      real(kind=8),    parameter ::   defepseigen   = 1.00D-4
      integer(kind=4)  pow,nrad,nang,irmesh
      character*(*)    wfnfile
c
      allocate (aomtopmod(ncent,nmo,nmo))
      allocate (eta(nmo,nmo))
      if (qtaim) then
        if (iover.lt.2.or.iover.eq.4) then
          iaom(1:ncent)=0
          do i=1,ncent
            read (lu18,*) iaom(i)
            if (iaom(i).gt.ncent) write (lr,17) iaom(i),ncent
            if (iaom(i).gt.ncent) stop
            if (largwr) write (lw,16) iaom(i)
            ii=iaom(i)
            if (iover.eq.0) then
              read (lu18,80) ((aom(ii,m,j),m=1,j),j=1,nmo) ! PROMOLDEN
            elseif (iover.eq.4) then
              do m=1,nmo
                read (lu18,81) (aom(iaom(i),j,m),j=1,m)     ! AIMALL
              enddo
            else
              do m=1,nmo
                read (lu18,82) (aom(iaom(i),j,m),j=1,m)     ! AIMPAC
              enddo
            endif
          enddo
c
c.........topmod
c
        elseif (iover.eq.2.or.iover.eq.3) then    !topmod
          ncent=nbasins
          read (lu18) (namebas(i),i=1,ncent)
          if (iover.eq.2) then
            write(lw,*) '# USING TOPMOD ELF BASINS'
            write (lw, *) '# Number of elf basins, names', ncent
            write (lw,*) (namebas(i),i=1,ncent)
          endif
          do i=1,ncent
            iaom(i)=i
            read (lu18) ((aomtopmod(iaom(i),j,m),m=j,nmo),j=1,nmo)
          enddo
          if (iover.eq.3) then    !aim
            ncent=natoms
            read (lu18) (namebas(i),i=1,ncent)
            write(lw,*) '# USING TOPMOD AIM BASINS'
            write (lw, *) '# Number of qtaim basins, names', ncent
            write (lw,*) (namebas(i),i=1,ncent)
            do i=1,ncent
              iaom(i)=i
              read (lu18) ((aomtopmod(iaom(i),j,m),m=j,nmo),j=1,nmo)
            enddo
          endif
          aom=aomtopmod
        endif
        close (lu18)
      else
c
c-------Orthogonalize exactly canonical and natural MOs using OFMO
c
        write (lw,10)
        call ofmortho (coef,eta,nmo,nprims,warno,lw)
        if (warno) stop ' # edf.f: Singular matrix in ofmortho.f'
        forall (i=1:ncent) iaom(i)=i
        if (mulliken) then
          call aomulliken (aom,0,.false.,lw,wfnfile)
        elseif (mindef) then
          call aomindef (aom,0,.false.,lw,wfnfile)
        elseif (mindefrho.or.netrho.or.promrho.or.hesel) then
          if (mindefrho) irho=1
          if (netrho)    irho=2
          if (promrho)   irho=3
          if (hesel)     irho=4
          allocate (cc(nmo+nmo,nprims))
          cc=coef
          call aomindefrho (aom,alphahesel,betahesel,
     &         irho,nrad,nang,irmesh,rhopow,.false.,cc,lw,wfnfile)
          deallocate (cc)
        elseif (lowdin) then
          call aomlowdin (aom,0,.false.,lw,wfnfile)
        elseif (becke) then
          allocate (cc(nmo,nprims))
          cc(1:nmo,1:nprims)=coef(nmo+1:nmo+nmo,1:nprims)
          call aombecke 
     &      (aom,nrad,nang,pow,irmesh,.false.,cc,lw,wfnfile)
          deallocate (cc)
        endif
      endif
c
c.....Test if the AOM has been read in or computed for all the basins.
c
      notaom=.false.
      do i=1,ncent
        j=1
        okaom=.false.
        do while (j.le.ncent.and.(.not.okaom))
          if (iaom(j).eq.i) okaom=.true.
          j=j+1
        enddo
        if (.not.okaom) then
          notaom=.true.
          write (lr,15) i
          stop
        endif
 14   enddo
      if (notaom) stop ' # edf.f: !! Bad definition of the AOM'
c
c.....Fill in the lower part of the aom() array.
c
      if (qtaim) then
        do m=1,nmo
          do j=1,m
            aom(1:ncent,m,j)=aom(1:ncent,j,m)
          enddo
        enddo
      endif

      deallocate (aomtopmod)
      deallocate (eta)
c
c.....Formats
c
 33   format (1x,'# AOM file                          = ',a)
 34   format (1x,'# Wave Function file                = ',a)
 55   format (1x,'# AOM read in with ',a,' format')
 103  format (1x,A4,4x,I3,3(2x,F15.8),4X,F5.1)
 80   format (6(1x,e16.10))
 81   format (6(1x,f13.10))
 82   format (8f10.6)
 501  format (1x,69('-'))
 432  format (1x,'# Atom number ',I3,' of group ',I2,
     &        1x,'is equal to atom number ',I3,' of group ',I2)
 15   format (/1x,'# Atomic overlap matrix not read for atom ',I3,/)
 16   format (1x,'# Reading AOM for atom ',I3)
 17   format (1x,'# Fatal error reading aom() array:',/,
     &        1x,'# Trying to read aom() for atom ',I3,
     &        ', and there are only ',I3,' atoms')
 21   format (1x,'# Number of electrons               = ',I10,/,
     &        1x,'# Number of molecular orbitals      = ',I10)
 444  format (6(2x,F12.6))
 30   format (10(1x,F12.6))
 335  format (1x,'# AOM is renormalized, the reference atom is ',I5)
 222    format (1x,'# FRAGMENT ',I3,' has ',I3,
     &             ' Localized CORE MOs: ',40I3)
 225    format (1x,'# <LOC_',I3,'|LOC_',I3,'>_',I2,' = ',F16.10)
 22   format (' # Probability = ',F12.6,4x,'for RSRS ',20I4)
 23   format (' # Approximate EDF assuming that electron populations')
 544  format (' # !! AOM is not Ok !! : SUM_n AOM(n,',2I3,') = ',E15.8)
 334  format (/,' # Second & fourth moment orbital spread of ',a)
 10   format (' # MOs are orthogonalized using the OFMO method')
 100  format (' # edf.f: File ',a,' NOT FOUND')
 101  format (' # edf.f: AOMNORM can not be used in UHF calculations',
     & /,' # edf.f: aomnorm ---> false')
 223  format (//1x,'# ',/' # ',I3,' Valence MOs',/,100(1x,'#',20I3,/))

      return
      end
