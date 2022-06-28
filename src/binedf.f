      SUbroutine binedf
     &    (ranmin,ranmax,
     &     epsdet,probcut,epsproba,npra,nprb,nprob,ngroup,
     &     ndets,nmo,ncore,nelact,nel,moval,ival,mocogrp,nelv,malv,
     &     mbev,ifilc,stdout,stderr,wfnfile,sg,resnca,resncb,largwr,
     &     doentropy,orderp,mem)
c
c.......................................................................
c
      USE        space_for_cidet
      include   'implicit.inc'
      include   'param.inc'
      include   'constants.inc'
      include   'lengrec.inc'
      include   'mline.inc'
c
      real   (kind=8), allocatable,dimension (:,:,:) :: ovea,oveb
      real   (kind=8), allocatable,dimension (:,:)   :: am,bm,oveaa
      real   (kind=8), allocatable,dimension (:,:)   :: ovebb
      real   (kind=8), allocatable,dimension (:,:)   :: tpowa,tpowb
      real   (kind=8), allocatable,dimension (:,:)   :: xlis
      real   (kind=8), allocatable,dimension (:)     :: proba,probb
      real   (kind=8), allocatable,dimension (:)     :: probt,probx
      real   (kind=8), allocatable,dimension (:)     :: pnew
      real   (kind=8), allocatable,dimension (:)     :: w,ww,xli
      real   (kind=8), allocatable,dimension (:)     :: p1,p1a,p1b,p2
      real   (kind=8), allocatable,dimension (:)     :: d1,p3,d2
      real   (kind=8), allocatable,dimension (:)     :: p2aa,p2bb
      real   (kind=8), allocatable,dimension (:)     :: p2ab,p2ba
      real   (kind=8), allocatable,dimension (:)     :: d1aa,d1bb
      real   (kind=8), allocatable,dimension (:)     :: d1ab,d1ba
      real   (kind=8), allocatable,dimension (:)     :: probord
      real   (kind=8), allocatable,dimension (:,:)   :: covmat
      real   (kind=8), allocatable,dimension (:,:)   :: covvec
      real   (kind=8), allocatable,dimension (:)     :: coveig
      real   (kind=8), allocatable,dimension (:)     :: wcov
      integer(kind=4), allocatable,dimension (:,:)   :: resnc,resncord
      integer(kind=4), allocatable,dimension (:)     :: ipvta,ipvtb,ind
      integer(kind=4), allocatable,dimension (:)     :: indxa,indxb
      integer(kind=4), allocatable,dimension (:)     :: iord,sob,sod
      integer(kind=4), allocatable,dimension (:)     :: idif,iorda
      integer(kind=4), allocatable,dimension (:)     :: npop,sodaux
      integer(kind=4), allocatable,dimension (:)     :: ordia,ordib
      integer(kind=4), allocatable,dimension (:)     :: ordja,ordjb
      integer(kind=4), allocatable,dimension (:)     :: kwa,kwb
      integer(kind=4), allocatable,dimension (:)     :: ioprob
      integer(kind=4), allocatable,dimension (:)     :: ndiffe
      integer(kind=4), allocatable,dimension (:)     :: edfa,edfb
      integer(kind=4), allocatable,dimension (:)     :: maxelec
      integer(kind=4), allocatable,dimension (:)     :: minelec
      real   (kind=8), allocatable,dimension (:,:)   :: pgroup
      real   (kind=8), allocatable,dimension (:,:,:) :: pg2
      integer(kind=8), allocatable,dimension (:)     :: deta,detb
      integer(kind=8), allocatable,dimension (:)     :: det1,det2
      integer(kind=8), allocatable,dimension (:,:,:) :: det
      integer(kind=4), allocatable,dimension(:,:)    :: nconfa,nconfb
      integer(kind=4), allocatable,dimension(:)      :: configa
      real(kind=8),    allocatable,dimension(:,:)    :: pacum,pbcum
      real(kind=8),    allocatable,dimension(:)      :: codet,codaux
      character(len=nmo),     allocatable, dimension(:) :: chdeta
      character(len=nmo),     allocatable, dimension(:) :: chdetb
      character(len=nmo+nmo), allocatable, dimension(:) :: chdet
c
      real   (kind=8)  sg(ngroup,nmo,nmo)
      real   (kind=8)  dumi,random,deter(2)
      integer(kind=4)  resnca(npra,ngroup),resncb(nprb,ngroup)
      integer(kind=4)  stdout,stderr
      integer(kind=4)  mocogrp(ngroup)
      integer(kind=4)  ival(moval)
      real   (kind=8)  indet
      integer(kind=4)  idum
      logical          inlist,orderp,largwr,doentropy,mem,found
      character(len=mline) wfnfile

      common /maxmin/ minimo,maximo

      integer(kind=8)  :: label
      character(len=nmo) :: chaux1,chaux2
      integer(kind=4)  nmobin,ishift
      character(len=1) dig(0:9)
      character(len=64) nchar
      integer(kind=4) :: mpab = 1000
      data dig/'0','1','2','3','4','5','6','7','8','9'/
!
!-----------------------------------------------------------------------
!
      call timer (2,ibinedf,'_binedf   ',-1)
!
!     ordia() and ordib() are the order numbers of the ALPHA and BETA
!     spin-orbitals (SO) occupied in each Slater determinant (detS).
!
!     kwa() and kwa() signal the ALPHA or BETA configuration to which
!     each detS belongs.
!
!
      if (.not.allocated(ordia)) then
        allocate (ordia(nmo),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ordia()'
      endif
      if (.not.allocated(ordib)) then
        allocate (ordib(nmo),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ordib()'
      endif
      if (.not.allocated(kwa)) then
        allocate (kwa(ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate kwa()'
      endif
      if (.not.allocated(kwb)) then
        allocate (kwb(ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate kwb()'
      endif
!
!     The core SOs are assumed to be occupied in all the detS.
!
      do mm=1,ncore
        ordia(mm)=mm
        ordib(mm)=mm
      enddo
!
!     INTEGER*8 numbers necessary to codify an ALPHA or BETA config
!
      nmobin = (nmo/64 + 1)
!
!     det(1:nmobin,i,j) are the integer*8 numbers that codify the ALPHA 
!     (i=1) or BETA (i=2) configuration of the jth detS.
!
      if (.not.allocated(det)) then
        allocate (det(nmobin,2,ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate det()'
      endif
!
!     chdeta(j) is a character variable used to codify the ALPHA
!     configuration of jth detS. The same meaning for chdetb(j) 
!     in the case of the BETA configuration. chdet(j) is simply
!     chdet(j)=chdeta(j)//chdetb(j)
!
      if (.not.allocated(chdeta)) then
        allocate (chdeta(ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate chdeta()'
      endif
      if (.not.allocated(chdetb)) then
        allocate (chdetb(ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate chdetb()'
      endif
      if (.not.allocated(chdet)) then
        allocate (chdet(ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate chdet()'
      endif
!
!     Mixing coefficients of all the detS
!
      if (.not.allocated(codet)) then
        allocate (codet(ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate codet()'
      endif
!
!     Some auxiliary arrays
!
      if (.not.allocated(deta)) then
        allocate (deta(nmobin),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate deta()'
      endif
      if (.not.allocated(detb)) then
        allocate (detb(nmobin),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate detb()'
      endif
      if (.not.allocated(det1)) then
        allocate (det1(nmobin),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate det1()'
      endif
      if (.not.allocated(det2)) then
        allocate (det2(nmobin),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate det2()'
      endif
      if (.not.allocated(codaux)) then
        allocate (codaux(ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate codaux()'
      endif
      if (.not.allocated(sob)) then
        allocate (sob(nelact),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate sob()'
      endif
      if (.not.allocated(sod)) then
        allocate (sod(nelact),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate sod()'
      endif
      if (.not.allocated(iord)) then
        allocate (iord(nelact),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate iord()'
      endif
      if (.not.allocated(idif)) then
        allocate (idif(nelact),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate idif()'
      endif
      if (.not.allocated(sodaux)) then
        allocate (sodaux(nelact),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate sodaux()'
      endif
!
!     Determine the number of ALPHA and BETA valence electrons (MAL,MBE)
!     (these numbers must be the same in all the detS)
!
      read (ifilc,rec=1) (cidet(m),m=0,nelact)
      mal=count(cidet(1:nelact)>0)
      mbe=count(cidet(1:nelact)<0)
!
!     Slater determinants are ordered
!
      if (ndets.gt.10) then
        iper=0
        shundred=dble(ndets)/10D0
        xhundred=shundred
        indet=0d0
      endif
      do i=1,ndets
        indet=indet+1d0
        if (ndets.gt.10.and.indet.ge.xhundred) then
          iper=iper+10
          write (stderr,475) iper
          call flush (stderr)
          xhundred=xhundred+shundred
        endif
        read (ifilc,rec=i) codet(i),(cidet(m),m=1,nelact)
        forall (m=1:nelact) iord(m)=m
        do m=1,nelact
          mm=int(cidet(m))
          if (mm.gt.0) then
            sob(m)=mm
          else
            sob(m)=mm+9999
          endif
        enddo
        call iqcksort (sob,iord,nelact,1,nelact)
        forall (m=1:mal) sod(m)=sob(iord(m))
        nn=0
        do m=mal+1,nelact
          iom=iord(m)
          sob(iom)=sob(iom)-9999
          sod(nelact-nn)=sob(iom)
          nn=nn+1
        enddo
        sodaux(1:nelact)=sod(1:nelact)
        call colord (sob,sodaux,nelact,nelact,nper,nsig,ndif,idif)
        if (ndif.gt.0) then
          stop ' # binedf.f: NDIF should be 0 at this point'
        endif
!
!       The phase of this detS has changed if nsig=-1
!
        codet(i) = codet(i) * dble(nsig)
        ordia(ncore+1:ncore+mal)=ncore+sod(1:mal)
        ordib(ncore+1:ncore+mbe)=ncore-sod(mal+1:nelact)
!
!       Codifying the detS using character variables
!
        forall (ic=1:nmo) chdeta(i)(ic:ic)='0'
        forall (ic=1:nmo) chdetb(i)(ic:ic)='0'
        do ic=1,ncore+mal
          chdeta(i)(nmo-ordia(ic)+1:nmo-ordia(ic)+1)='1'
        enddo
        do ic=1,ncore+mbe
          chdetb(i)(nmo-ordib(ic)+1:nmo-ordib(ic)+1)='1'
        enddo
        chdet(i)=chdeta(i)//chdetb(i)
      enddo
      xnorm = 1d0/sqrt(dot_product(codet,codet))

      if (allocated(chdeta)) then
        deallocate (chdeta,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate chdeta()'
      endif
!
!     Slater determinants are ordered
!
      if (.not.allocated(iorda)) then
        allocate (iorda(ndets),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate iorda()'
      endif
      forall (i=1:ndets) iorda(i)=i
      call chsort (chdet, iorda, nmo+nmo, 1, ndets, nstack)
!
!     Determining the number of different ALPHA configurations (NPA)
!
      npa=1
      do io=2,ndets
        icur=iorda(io)
        iant=iorda(io-1)
        chaux1(1:nmo)=chdet(icur)(1:nmo)
        chaux2(1:nmo)=chdet(iant)(1:nmo)
        if (chaux1.ne.chaux2) npa=npa+1
      enddo
!
!     Determining the array configa(1:NPA), that contains the number of
!     times that a given ALPHA configuration is repeated in the full
!     list of detS
!
      if (.not.allocated(configa)) then
        allocate (configa(npa),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate configa()'
      endif
      k=1
      configa(k)=1
      forall (io=1:ndets) codaux(io)=codet(iorda(io))
      forall (io=1:ndets) codet(io)=codaux(io)
      if (allocated(codaux)) then
        deallocate (codaux,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate codaux()'
      endif
      do io=1,ndets
        icur=iorda(io)
        if (io.gt.1) then
          iant=iorda(io-1)
          chaux1(1:nmo)=chdet(icur)(1:nmo)
          chaux2(1:nmo)=chdet(iant)(1:nmo)
          if (chaux1.ne.chaux2) then
            k=k+1
            configa(k)=1
          else
            configa(k)=configa(k)+1
          endif
        endif
        ka=0
        kb=0
        do ic=1,nmo
          nic=nmo-ic+1
          if (chdet(icur)(nic:nic).eq.'1') then
            ka=ka+1
            ordia(ka)=ic
          endif
          if (chdet(icur)(nic+nmo:nic+nmo).eq.'1') then
            kb=kb+1
            ordib(kb)=ic
          endif
        enddo
!
!       Re-Writing the detS in the original file, now ordered as follows: 
!
!       First the configa(1) detS with the ALPHA configuration 1, 
!       then  the configa(2) detS with the ALPHA configuration 2, etc
!       Within each subset of detS with the same ALPHA configuration,
!       the detS are ordered as well according to the BETA config.
!
        write (ifilc,rec=io) codet(io),
     .    (+dble(ordia(ncore+ka)-ncore),ka=1,mal),
     .    (-dble(ordib(ncore+kb)-ncore),kb=1,mbe)
      enddo

      if (allocated(iorda)) then
        deallocate (iorda,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate iorda()'
      endif
      if (allocated(chdet)) then
        deallocate (chdet,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate chdet()'
      endif
!
!.....Renormalize the WFN
!
      codet(1:ndets) = codet(1:ndets) * xnorm
!
!     Codifying strictly different ALPHA configurations
!
      det = 0_8
      ncumconfa = 0
      do k=1,npa
        i=ncumconfa+1
        read (ifilc,rec=i) dummy,(cidet(m),m=1,mal)
        ordia(ncore+1:ncore+mal)=+int(+cidet(1:mal))+ncore
        deta=0_8
        do ia=1,ncore+mal
          ibin=(ordia(ia)/64 + 1)
          deta(ibin)=ibset(deta(ibin),mod(ordia(ia)-1,63))
        enddo
        det(1:nmobin,1,k)=deta(1:nmobin)
        nkwa1=ncumconfa+1
        nkwa2=ncumconfa+configa(k)
        kwa(nkwa1:nkwa2) = k
        ncumconfa=ncumconfa+configa(k)
      enddo

      if (largwr) then 
        write (stdout,1002)
        if (.not.allocated(ndiffe)) then
          allocate (ndiffe(1:npa),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate ndiffe()'
        endif
        do k1=1,npa
          ndiffe=0
          do k2=1,k1
            nexci=popcnt(xor(det(1,1,k1),det(1,1,k2)))
            do l=2,nmobin
              nexci=nexci+popcnt(xor(det(l,1,k1),det(l,1,k2)))
            enddo
            ndiffe(k2)=ishft(nexci,-1)
          enddo
          write (stdout,1001) (ndiffe(k2),k2=1,k1)
        enddo
        if (allocated(ndiffe)) then
          deallocate (ndiffe,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ndiffe()'
        endif
      endif

      if (largwr) then
        write (stdout,*) '#'
        write (stdout,*) '# ALPHA AND BETA ACTIVE CONFIGURATIONS'
      endif
!
!     Codifying strictly different BETA configurations
!
      if (ndets.gt.10) then
        iper=0
        shundred=dble(ndets)/10D0
        xhundred=shundred
        indet=0d0
      endif

      call timer (2,icody,'_codifydet',-1)
      npb = 0
      do i=1,ndets
        indet=indet+1d0
        if (ndets.gt.10.and.indet.ge.xhundred) then
          iper=iper+10
          write (stderr,473) iper
          call flush (stderr)
          xhundred=xhundred+shundred
        endif
        read (ifilc,rec=i) dummy,(cidet(m),m=1,nelact)
        ordib(ncore+1:ncore+mbe)=+int(-cidet(mal+1:nelact))+ncore
        detb=0_8
        do ib=1,ncore+mbe
          ibin=(ordib(ib)/64 + 1)
          detb(ibin)=ibset(detb(ibin),mod(ordib(ib)-1,63))
        enddo
!
!       The first configa(1) detS have a different BETA configuration
!
        if (i.le.configa(1)) then
          npb=npb+1
          kwb(i)=npb
          det(1:nmobin,2,npb)=detb(1:nmobin)
          forall (ic=1:nmo) chdetb(npb)(ic:ic)='0'
          ishift = 0
          do ibin = 1,nmobin
            label = detb(ibin)
            do while (label /= 0_8)
              ic=trailz(label) + ishift + 1
              chdetb(npb)(nmo-ic+1:nmo-ic+1)='1'
              label = iand(label,label-1_8)
            enddo
            ishift = ishift + 63
          enddo
        else
          forall (ic=1:nmo) chaux1(ic:ic)='0'
          ishift = 0
          do ibin = 1,nmobin
            label = detb(ibin)
            do while (label /= 0_8)
              ic = trailz(label) + ishift + 1
              chaux1(nmo-ic+1:nmo-ic+1)='1'
              label = iand(label,label-1_8)
            enddo
            ishift = ishift + 63
          enddo
!
!         Determine if this BETA configuration is already in  the list
!
          call searchstring (chdetb,chaux1,1,npb,ifound,found)
          if (found) then
            kwb(i)=ifound
            cycle
          endif
          npb=npb+1
          kwb(i)=maximo
          do j=1,i-1
            if (kwb(j).ge.maximo) kwb(j)=kwb(j)+1
          enddo
          do ic=npb,maximo+1,-1
            det(1:nmobin,2,ic)=det(1:nmobin,2,ic-1)
            chdetb(ic)=chdetb(ic-1)
          enddo
          det(1:nmobin,2,maximo)=detb(1:nmobin)
          chdetb(maximo)=chaux1
        endif
      enddo
      call timer (4,icody,'_codifydet',-1)
!
!     Deallocate some arrays
!
      if (allocated(chdet)) then
        deallocate (chdet,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate chdet()'
      endif
      if (allocated(chdeta)) then
        deallocate (chdeta,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate chdeta()'
      endif
      if (allocated(chdetb)) then
        deallocate (chdetb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate chdetb()'
      endif
      if (allocated(iord)) then
        deallocate (iord,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate iord()'
      endif
      if (allocated(det1)) then
        deallocate (det1,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate det1()'
      endif
      if (allocated(det2)) then
        deallocate (det2,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate det2()'
      endif
      if (allocated(deta)) then
        deallocate (deta,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate deta()'
      endif
      if (allocated(detb)) then
        deallocate (detb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate detb()'
      endif

      do ispin=1,2
        if (ispin.eq.1) then
          write (stdout,*) '# NUMBER OF alpha CONFIGS = ',npa
          npis=npa
        else
          write (stdout,*) '# NUMBER OF beta  CONFIGS = ',npb
          npis=npb
        endif
        if (nmo.lt.64) then
          do k=1,npis
            do ibin = 1,nmobin
              forall (l=1:nmo) nchar(l:l)='0'
              label = det(ibin,ispin,k)
              i=1
              do while (label.ge.2)
                nchar(i:i)=dig(mod(label,2))
                i=i+1
                label = ishft(label,-1)
              enddo
              nchar(i:i)=dig(label)
              write (stdout,2323) k,(nchar(j:j),j=nmo,1,-1)
            enddo
          enddo
        endif
      enddo
      if (allocated(ordia)) then
        deallocate (ordia,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ordia()'
      endif
      if (allocated(ordib)) then
        deallocate (ordib,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ordib()'
      endif
!
!-----------------------------------------------------------------------
!
!-----Memory cannot be used in case that npa*(npa+1)/2 or npb*(npb+1)/2
!     is greater than mpab. If that happens, the mpab value increased
!     by 1000 until both npa*(npa+1)/2 and npb*(npb+1)/2 are smaller
!     than MPAB. If after this, the pacum() or pbcum() arrays can not
!     be allocated, the mem=.true. option can not be used and men is set
!     to .false.
!
!-----Redardless the input value, the mem variable is always set to 
!     .true. at the beginning
!-----------------------------------------------------------------------
!
      mem = .true.
 333  if (npa*(npa+1)/2.gt.mpab .or. npb*(npb+1)/2.gt.mpab) then
        mpab = mpab + 1000
        goto 333
      else
        allocate (pacum(mpab,npra),stat=ier)
        if (ier.ne.0) mem = .false.
        allocate (pbcum(mpab,nprb),stat=ier)
        if (ier.ne.0) mem = .false.
        if (mem) write (stdout,334) mpab
      endif      
 334  format (' # Actual value use for MPAB = ',I8)
!
!.....Random numbers generation for alpha block.
!
      call semilla (idum)
      idum=-idum
      dumi=random(idum)
      if (.not.allocated(am)) then
        allocate (am(npra,npra),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate am()'
      endif
      if (.not.allocated(tpowa)) then
        allocate (tpowa(npra,ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate tpowa()'
      endif
      if (.not.allocated(ipvta)) then
        allocate (ipvta(npra),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ipvta()'
      endif
      if (.not.allocated(w)) then
        allocate (w(npra),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate w()'
      endif

      rmami=ranmax-ranmin
      do
        am=zero
        tpowa=zero
        do i=1,npra
          do k=1,ngroup
            tpowa(i,k)=rmami*random(idum)+ranmin
          enddo
        enddo
        do i=1,npra
          do j=1,npra
            aco=one
            do k=1,ngroup-1
              aco=aco*tpowa(i,k)**resnca(j,k)
            enddo
            am(i,j)=aco
          enddo
        enddo
        call lulapack (am,ipvta,npra,info)
        if (info.eq.0) exit
      enddo
      if (allocated(w)) then
        deallocate (w,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate w()'
      endif
!
!.....run over all pairs of different alpha configs
!
      if (.not.allocated(ordia)) then
        allocate (ordia(nmo),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ordia()'
      endif
      if (.not.allocated(ordja)) then
        allocate (ordja(nmo),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ordja()'
      endif
      iper=0
      shundred=dble(npa)*(dble(npa)+1)/20D0
      xhundred=shundred
      indet=0d0

      if (.not.mem) then    
        indproba=59
        open (unit   = indproba,
     .        file   = wfnfile(1:leng(wfnfile))//'.probala.dat',
     .        access ='direct',
     .        recl   = RecLength*2*npra,
     .        form   = 'unformatted')
      endif

      if (.not.allocated(proba)) then
        allocate (proba(npra),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate proba()'
      endif
      if (.not.allocated(ovea)) then
        allocate (ovea(ngroup,malv,malv),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ovea()'
      endif
      if (.not.allocated(oveb)) then
        allocate (oveb(ngroup,mbev,mbev),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate oveb()'
      endif
      if (.not.allocated(oveaa)) then
        allocate (oveaa(malv,malv),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate oveaa()'
      endif
      if (.not.allocated(ovebb)) then
        allocate (ovebb(mbev,mbev),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ovebb()'
      endif
      if (.not.allocated(indxa)) then
        allocate (indxa(malv),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate indxa()'
      endif
      if (.not.allocated(ww)) then
        allocate (ww(moval),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ww()'
      endif

      call timer (2,ilinears,'_linearsys',-1)
      do i=1,npa
        ishift = 0
        k=0
        do ibin = 1,nmobin
          label = det(ibin,1,i)
          do while (label /= 0_8)
            k=k+1
            ordia(k)=trailz(label) + ishift + 1
            label = iand(label,label-1_8)
          enddo
          ishift = ishift + 63
        enddo
        do j=1,i
          ishift = 0
          k=0
          do ibin = 1,nmobin
            label = det(ibin,1,j)
            do while (label /= 0_8)
              k=k+1
              ordja(k)=trailz(label) + ishift + 1
              label = iand(label,label-1_8)
            enddo
            ishift = ishift + 63
          enddo
          indet=indet+1d0
          if (indet.ge.xhundred) then
            iper=iper+10
            if (npa.gt.10) write (stderr,470) '(alpha)',iper
            call flush (stderr)
            xhundred=xhundred+shundred
          endif
!
!.........Overlap matrices in ALPHA block.
!
          do m=1,malv
            ivm=ordia(ival(m))
            do k=1,malv
              ivk=ordja(ival(k))
              ovea(1:ngroup,m,k)=sg(1:ngroup,ivm,ivk)
            enddo
          enddo
 
          if (malv.gt.izero) then
            do n=1,npra
              oveaa(1:malv,1:malv)=ovea(ngroup,1:malv,1:malv)
              do igr=1,ngroup-1
                oveaa(1:malv,1:malv)=oveaa(1:malv,1:malv) + 
     .             tpowa(n,igr) * ovea(igr,1:malv,1:malv)
              enddo
              call detlapack (oveaa,indxa,malv,info,proba(n))
            enddo
            call dgetrs ('N',npra,1,am,npra,ipvta,proba,npra,info)
          else
            npra=1
            proba(1)=one
          endif
          irec=i*(i-1)/2+j

          if (.not.mem) then
            write (indproba,rec=irec)(proba(n),n=1,npra)
          else
            if (irec.gt.mpab) then
              write (stderr,468)
              stop 
            endif
            pacum(irec,1:npra)=proba(1:npra)
          endif
        enddo
      enddo

      if (allocated(tpowa)) then
        deallocate (tpowa,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate tpowa()'
      endif
      if (allocated(am)) then
        deallocate (am,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate am()'
      endif
      if (allocated(ipvta)) then
        deallocate (ipvta,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ipvta()'
      endif
      if (allocated(ordia)) then
        deallocate (ordia,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ordia()'
      endif
      if (allocated(ordja)) then
        deallocate (ordja,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ordja()'
      endif
      if (allocated(indxa)) then
        deallocate (indxa,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate indxa()'
      endif
      if (allocated(ovea)) then
        deallocate (ovea,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ovea()'
      endif
!
!.....Random numbers generation for beta block.
!
*     write (stdout,44) '# Linear System Solver'
      if (.not.allocated(bm)) then
        allocate (bm(nprb,nprb),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate bm()'
      endif
      if (.not.allocated(tpowb)) then
        allocate (tpowb(nprb,ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate tpowb()'
      endif
      if (.not.allocated(ipvtb)) then
        allocate (ipvtb(nprb),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ipvtb()'
      endif
      if (.not.allocated(w)) then
        allocate (w(nprb),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate w()'
      endif
      do
        bm=zero
        tpowb=zero
        do i=1,nprb
          do k=1,ngroup
            tpowb(i,k)=rmami*random(idum)+ranmin
          enddo
        enddo
        do i=1,nprb
          do j=1,nprb
            aco=one
            do k=1,ngroup-1
              aco=aco*tpowb(i,k)**resncb(j,k)
            enddo
            bm(i,j)=aco
          enddo
        enddo
        call lulapack (bm,ipvtb,nprb,info)
        if (info.eq.0) exit
      enddo
      if (allocated(w)) then
        deallocate (w,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate w()'
      endif
!
!.....run over all pairs of different beta configs
!
      if (.not.allocated(ordib)) then
        allocate (ordib(nmo),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ordib()'
      endif
      if (.not.allocated(ordjb)) then
        allocate (ordjb(nmo),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ordjb()'
      endif
      if (.not.allocated(indxb)) then
        allocate (indxb(mbev),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate indxb()'
      endif
      iper=0
      shundred=dble(npb)*(dble(npb)+1)/20D0
      xhundred=shundred
      indet=0d0

      if (.not.mem) then
        indprobb=61
        open (unit   = indprobb,
     .        file   = wfnfile(1:leng(wfnfile))//'.probalb.dat',
     .        access ='direct',
     .        recl   = RecLength*2*nprb,
     .        form   ='unformatted')
      endif

      if (.not.allocated(probb)) then
        allocate (probb(nprb),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate probb()'
      endif
      do m=1,ncore
        ordib(m)=m
        ordjb(m)=m
      enddo
      do i=1,npb
        ishift = 0
        k=0
        do ibin = 1,nmobin
          label = det(ibin,2,i)
          do while (label /= 0_8)
            k=k+1
            ordib(k)=trailz(label) + ishift + 1
            label = iand(label,label-1_8)
          enddo
          ishift = ishift + 63
        enddo
        do j=1,i
          ishift = 0
          k=0
          do ibin = 1,nmobin
            label = det(ibin,2,j)
            do while (label /= 0_8)
              k=k+1
              ordjb(k)=trailz(label) + ishift + 1
              label = iand(label,label-1_8)
            enddo
            ishift = ishift + 63
          enddo
          indet=indet+1d0
!
!..........Testing how the most costly part of the run is going on.
!
           if (indet.ge.xhundred) then
             iper=iper+10
             if (npb.gt.10) write (stderr,470) '(beta) ',iper
             call flush (stderr)
             xhundred=xhundred+shundred
           endif
!
!..........Overlap matrices in BETA  block.
!
           do m=1,mbev
             ivm=ordib(ival(m))
             do k=1,mbev
               ivk=ordjb(ival(k))
               oveb(1:ngroup,m,k)=sg(1:ngroup,ivm,ivk)
             enddo
           enddo

           if (mbev.gt.izero) then
             do n=1,nprb
               ovebb(1:mbev,1:mbev)=oveb(ngroup,1:mbev,1:mbev)
               do igr=1,ngroup-1
                 ovebb(1:mbev,1:mbev)=ovebb(1:mbev,1:mbev) + 
     .              tpowb(n,igr) * oveb(igr,1:mbev,1:mbev)
               enddo
               call detlapack (ovebb,indxb,mbev,info,probb(n))
             enddo
             call dgetrs ('N',nprb,1,bm,nprb,ipvtb,probb,nprb,info)
           else
             nprb=1
             probb(1)=one
           endif
           irec=i*(i-1)/2+j

           if (.not.mem) then
             write (indprobb,rec=irec)(probb(n),n=1,nprb)
           else
             if (irec.gt.mpab) then
               write (stderr,468)
               stop 
             endif
             pbcum(irec,1:nprb)=probb(1:nprb)
           endif
        enddo
      enddo
      if (allocated(tpowb)) then
        deallocate (tpowb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate tpowb()'
      endif
      if (allocated(bm)) then
        deallocate (bm,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate bm()'
      endif
      if (allocated(ipvtb)) then
        deallocate (ipvtb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ipvtb()'
      endif
      if (allocated(ordib)) then
        deallocate (ordib,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ordib()'
      endif
      if (allocated(ordjb)) then
        deallocate (ordjb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ordjb()'
      endif
      if (allocated(oveb)) then
        deallocate (oveb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate oveb()'
      endif
      if (allocated(oveaa)) then
        deallocate (oveaa,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate oveaa()'
      endif
      if (allocated(indxb)) then
        deallocate (indxb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate indxb()'
      endif
      if (allocated(ww)) then
        deallocate (ww,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ww()'
      endif
 
      iper=0
      shundred=dble(ndets)*dble((ndets)+1)/200D0
      xhundred=shundred
      indet=0d0
      call timer (4,ilinears,'_linearsys',-1)
!
!.....Compute probabilities
!
      if (.not.mem) then
        if (.not.allocated(probx)) then
          allocate (probx(nprb),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate probx()'
        endif
      endif
      if (.not.allocated(nconfa)) then
        allocate (nconfa(npa,2),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate nconfa()'
      endif

      call timer (2,iedfsplit,'_splitedf ',-1)
      iper=0
      shundred=dble(npa)*(dble(npa)+1)/20D0
      xhundred=shundred
      indet=0d0

      if (.not.allocated(probt)) then
        allocate (probt(npra*nprb),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate probt()'
      endif
      probt=zero
      write (stdout,44) '# Spin-Splitted probabilities'
      ncumconfa1=0
      do m1=1,npa
        nconf1=configa(m1)
        forall (i=1:nconf1) nconfa(i,1)=ncumconfa1+i
        ncumconfa2=0
        do m2=1,m1
          nconf2=configa(m2)
          forall (i=1:nconf2) nconfa(i,2)=ncumconfa2+i
          indet=indet+1d0
          if (indet.ge.xhundred) then
            iper=iper+10
            if (npa.gt.10) write (stderr,471) iper
            call flush (stderr)
            xhundred=xhundred+shundred
          endif
          ireca=m1*(m1-1)/2+m2

          if (.not.mem) then
            read (indproba,rec=ireca)(proba(k),k=1,npra)
          else
            proba(1:npra)=pacum(ireca,1:npra)
          endif
          if (maxval(proba).ge.epsproba) then
            probb(1:nprb)=0D0
            do n1=1,nconf1
              i=nconfa(n1,1)
              do n2=1,nconf2
                j=nconfa(n2,2)
                cd=codet(i)*codet(j)
                if (m1.ne.m2) cd=cd+cd
                if (abs(cd).gt.abs(epsdet)) then
                  iminb=min(kwb(i),kwb(j))
                  imaxb=max(kwb(i),kwb(j))
                  irecb=imaxb*(imaxb-1)/2+iminb
  
                  if (.not.mem) then
                    read (indprobb,rec=irecb)(probx(k),k=1,nprb)
                    probb(1:nprb)=probb(1:nprb)+cd*probx(1:nprb)
                  else
                    probb(1:nprb)= 
     .              probb(1:nprb)+cd*pbcum(irecb,1:nprb)
                  endif
                endif
              enddo
            enddo
            ij=0
            do ia=1,npra
              do ib=1,nprb
                ij=ij+1
                probt(ij)=probt(ij)+proba(ia)*probb(ib)
              enddo
            enddo
          endif
          ncumconfa2=ncumconfa2+configa(m2)
        enddo
        ncumconfa1=ncumconfa1+configa(m1)
      enddo
      call timer (4,iedfsplit,'_splitedf ',-1)

      if (.not.mem) then
        if (allocated(probx)) then
          deallocate (probx,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate probx()'
        endif
      endif

      if (allocated(codet)) then
        deallocate (codet,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate codet()'
      endif
      if (mem) then
        if (allocated(pacum)) then
          deallocate (pacum,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate pacum()'
        endif
        if (allocated(pbcum)) then
          deallocate (pbcum,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate pbcum()'
        endif
      endif
!
!-----Deallocate arrays related to the ALPHA and BETA configs
!
      if (allocated(nconfa)) then
        deallocate (nconfa,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate nconfa()'
      endif
!
!.....Write probabilities.
!
      if (.not.allocated(edfa)) then
        allocate (edfa(ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate edfa()'
      endif
      if (.not.allocated(edfb)) then
        allocate (edfb(ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate edfb()'
      endif

      write (stdout,200) ngroup,npra*nprb
      nprobg=0
      suma=zero
      sumtot=zero
      ij=0
      do i=1,npra
        edfa(1:ngroup)=resnca(i,1:ngroup)+mocogrp(1:ngroup)
        do j=1,nprb
          edfb(1:ngroup)=resncb(j,1:ngroup)+mocogrp(1:ngroup)
          ij=ij+1
          sumtot=sumtot+probt(ij)
          if (probt(ij).ge.probcut) then
            nprobg=nprobg+1
            suma=suma+probt(ij)
            write (stdout,100) probt(ij),(edfa(k),edfb(k),k=1,ngroup)
          endif
        enddo
      enddo
      write (stdout,106) suma, nprobg, probcut, sumtot

      if (allocated(edfa)) then
        deallocate (edfa,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate edfa()'
      endif
      if (allocated(edfb)) then
        deallocate (edfb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate edfb()'
      endif
!
!.....Average population of alpha and beta electrons in each group.
!     Store also in proba() and probb() arrays the alpha and beta 
!     probabilities.
!     
      proba=zero
      probb=zero
      ik=0
      do i=1,npra
        do k=1,nprb
          ik=ik+1
          proba(i)=proba(i)+probt(ik)
          probb(k)=probb(k)+probt(ik)
        enddo
      enddo
!
!.....Write EDF for ALPHA and BETA electrons independently
!
      write (stdout,203) 'ALPHA','BETA ',ngroup,'ALPHA',npra
      nprobg=0
      suma=zero
      sumtot=zero
      do i=1,npra
        sumtot=sumtot+proba(i)
        if (proba(i).ge.probcut) then
          nprobg=nprobg+1
          suma=suma+proba(i)
          write (stdout,100) proba(i),
     .      (resnca(i,igr)+mocogrp(igr),igr=1,ngroup)
        endif
      enddo
      if (allocated(proba)) then
        deallocate (proba,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate proba()'
      endif
      write (stdout,106) suma, nprobg, probcut, sumtot
      write (stdout,203) 'BETA ','ALPHA',ngroup,'BETA ',nprb
      nprobg=0
      suma=zero
      sumtot=zero
      do i=1,nprb
        sumtot=sumtot+probb(i)
        if (probb(i).ge.probcut) then
          nprobg=nprobg+1
          suma=suma+probb(i)
          write (stdout,100) probb(i),
     .      (resncb(i,igr)+mocogrp(igr),igr=1,ngroup)
        endif
      enddo
      if (allocated(probb)) then
        deallocate (probb,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate probb()'
      endif
      write (stdout,106) suma, nprobg, probcut, sumtot
c
c.....Compute spin-splitted localization and delocatization indices.
c
      if (.not.allocated(p1a)) then
        allocate (p1a(ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate p1a()'
      endif
      if (.not.allocated(p1b)) then
        allocate (p1b(ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate p1b()'
      endif
      if (.not.allocated(xlis)) then
        allocate (xlis(ngroup,2),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate xlis()'
      endif
      if (ngroup.gt.1) then
        ngpair=ngroup*(ngroup-1)/2
        if (.not.allocated(p2aa)) then
          allocate (p2aa(ngpair),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate p2aa()'
        endif
        if (.not.allocated(p2bb)) then
          allocate (p2bb(ngpair),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate p2bb()'
        endif
        if (.not.allocated(p2ab)) then
          allocate (p2ab(ngpair),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate p2ab()'
        endif
        if (.not.allocated(p2ba)) then
          allocate (p2ba(ngpair),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate p2ba()'
        endif
        if (.not.allocated(d1aa)) then
          allocate (d1aa(ngpair),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate d1aa()'
        endif
        if (.not.allocated(d1bb)) then
          allocate (d1bb(ngpair),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate d1bb()'
        endif
        if (.not.allocated(d1ab)) then
          allocate (d1ab(ngpair),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate d1ab()'
        endif
        if (.not.allocated(d1ba)) then
          allocate (d1ba(ngpair),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate d1ba()'
        endif
      endif
      p1a  = zero
      p1b  = zero
      p2aa = zero
      p2ab = zero
      p2ba = zero
      p2bb = zero
      d1aa = zero
      d1ab = zero
      d1ba = zero
      d1bb = zero
      xlis = zero
      np=0
      do ia=1,npra
        do ib=1,nprb
          np=np+1
          p=probt(np)
          do i=1,ngroup
            riai=resnca(ia,i)+mocogrp(i)
            ribi=resncb(ib,i)+mocogrp(i)
            p1a(i)=p1a(i)+p*riai
            p1b(i)=p1b(i)+p*ribi
            xlis(i,1)=xlis(i,1)+p*riai*riai
            xlis(i,2)=xlis(i,2)+p*ribi*ribi
          enddo
          if (ngroup.gt.1) then
            ipair=0
            do i=1,ngroup
              riai=resnca(ia,i)+mocogrp(i)
              ribi=resncb(ib,i)+mocogrp(i)
              do j=1,i-1
                riaj=resnca(ia,j)+mocogrp(j)
                ribj=resncb(ib,j)+mocogrp(j)
                ipair=ipair+1
                p2aa(ipair)=p2aa(ipair)+riai*riaj*p
                p2ab(ipair)=p2ab(ipair)+riai*ribj*p
                p2ba(ipair)=p2ba(ipair)+ribi*riaj*p
                p2bb(ipair)=p2bb(ipair)+ribi*ribj*p
              enddo
            enddo
          endif
        enddo
      enddo
 
      write (stdout,119)
      do i=1,ngroup
        mcoi=mocogrp(i)
        write (stdout,19) i,p1a(i)+mcoi,i,p1b(i)+mcoi
      enddo
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (stdout,29) 
     &       i,j,p2aa(ipair),i,j,p2ab(ipair),
     &       i,j,p2ba(ipair),i,j,p2bb(ipair)
          enddo
        enddo
      endif
  
      np=0
      do ia=1,npra
        do ib=1,nprb
          np=np+1
          p=probt(np)
          twop=p+p
          if (ngroup.gt.1) then
            ipair=0
            do i=1,ngroup
              do j=1,i-1
                ipair=ipair+1
                nai=resnca(ia,i)+mocogrp(i)
                naj=resnca(ia,j)+mocogrp(j)
                nbi=resncb(ib,i)+mocogrp(i)
                nbj=resncb(ib,j)+mocogrp(j)
                d1aa(ipair)=d1aa(ipair)-twop*(p1a(i)-nai)*(p1a(j)-naj)
                d1ab(ipair)=d1ab(ipair)-twop*(p1a(i)-nai)*(p1b(j)-nbj)
                d1ba(ipair)=d1ba(ipair)-twop*(p1b(i)-nbi)*(p1a(j)-naj)
                d1bb(ipair)=d1bb(ipair)-twop*(p1b(i)-nbi)*(p1b(j)-nbj)
              enddo
            enddo
          endif
        enddo
      enddo
 
      do i=1,ngroup
        xlis(i,1)=-(xlis(i,1)-p1a(i)*(p1a(i)+one))
        xlis(i,2)=-(xlis(i,2)-p1b(i)*(p1b(i)+one))
        terma=zero
        termb=zero
        if (p1a(i).gt.epsq8) terma=hundred*xlis(i,1)/p1a(i)
        if (p1b(i).gt.epsq8) termb=hundred*xlis(i,2)/p1b(i)
        write (stdout,69) i,i,xlis(i,1),terma,i,i,xlis(i,2),termb
      enddo
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (stdout,49) i,j,d1aa(ipair),i,j,d1ab(ipair),
     .                        i,j,d1ba(ipair),i,j,d1bb(ipair)
          enddo
        enddo
      endif
 
      if (ngroup.gt.1) then
        if (allocated(p2aa)) then
          deallocate (p2aa,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p2aa()'
        endif
        if (allocated(p2bb)) then
          deallocate (p2bb,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p2bb()'
        endif
        if (allocated(p2ab)) then
          deallocate (p2ab,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p2ab()'
        endif
        if (allocated(p2ba)) then
          deallocate (p2ba,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p2ba()'
        endif
        if (allocated(d1aa)) then
          deallocate (d1aa,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate d1aa()'
        endif
        if (allocated(d1bb)) then
          deallocate (d1bb,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate d1bb()'
        endif
        if (allocated(d1ab)) then
          deallocate (d1ab,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate d1ab()'
        endif
        if (allocated(d1ba)) then
          deallocate (d1ba,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate d1ba()'
        endif
      endif
      if (allocated(xlis)) then
        deallocate (xlis,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate xlis()'
      endif
!
!.....Obtain spinless multiple-fragment electron population covariances.
!
c
c-----Disccount core electrons from the population as required
c     by the ndelta routine.
c
      do i=1,ngroup
        p1a(i)=p1a(i)-mocogrp(i)
        p1b(i)=p1b(i)-mocogrp(i)
      enddo
c
c-----Undo this disccount
c
      do i=1,ngroup
        p1a(i)=p1a(i)+mocogrp(i)
        p1b(i)=p1b(i)+mocogrp(i)
      enddo
 
!
!.....Computes spinless probabilities from spin resolved probabilities.
!.....Obtain the number of spinless real space resonant structures.
!
      combi=one
      do i=0,nel-1
        combi=combi*dble(nel+ngroup-1-i)/dble(nel-i)
      enddo
      nprobt=int(combi+epsq10)
 
      if (.not.allocated(pnew)) then
        allocate (pnew(nprobt),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate pnew()'
      endif
      if (.not.allocated(resnc)) then
        allocate (resnc(nprobt,ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate resnc()'
      endif
      if (.not.allocated(ind)) then
        allocate (ind(ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ind()'
      endif
      pnew=zero
      n=0
      do ia=1,npra
        cycleib: do ib=1,nprb
          n=n+1
          if (n.eq.1) then
            np=1
            resnc(np,:)=resnca(ia,:)+resncb(ib,:)
            pnew(np)=probt(n)
          else
            ind(:)=resnca(ia,:)+resncb(ib,:)
            do j=1,np
              inlist=.true.
              do k=1,ngroup
                inlist=inlist.and.(ind(k).eq.resnc(j,k))
              enddo
              if (inlist) then
                pnew(j)=pnew(j)+probt(n)
                cycle cycleib
              endif
            enddo
            np=np+1
            pnew(np)=probt(n)
            resnc(np,:)=resnca(ia,:)+resncb(ib,:)
          endif
        enddo cycleib
      enddo
      if (allocated(ind)) then
        deallocate (ind,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ind()'
      endif
!
!.....Determine if probabilities are ordered or not.
!
      write (stdout,202) ngroup,np
      npnew=0
      sumtot=zero
      suma=zero
      do n=1,np
        sumtot=sumtot+pnew(n)
        if (pnew(n).ge.probcut) then
          npnew=npnew+1
          suma=suma+pnew(n)
        endif
      enddo
 
      if (npnew.le.nstack.and.orderp) then
!
!.......Order by increasing value the highest probabilities.
!
        allocate (ioprob(npnew),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate ioprob()'
        allocate (resncord(npnew,ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate resncord()'
        allocate (probord(npnew),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate probord()'
        npronew=0
        do n=1,np
          if (pnew(n).ge.probcut) then
            npronew=npronew+1
            probord(npronew)=pnew(n) 
            ioprob(npronew)=npronew
            resncord(npronew,:) = resnc(n,:)+2*mocogrp(:)
          endif
        enddo
        call qqsort (probord,ioprob,1,npnew,nstack)
        do n=npnew,1,-1
          m=ioprob(n)
          write (stdout,100) probord(m),(resncord(m,igr),igr=1,ngroup)
          nn=nn+1
        enddo
        deallocate (ioprob,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate ioprob()'
        deallocate (resncord,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate resncord()'
        deallocate (probord,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate probord()'
      else
!
!.......Probabilities are not ordered.
!
        do n=1,np
          if (pnew(n).ge.probcut) write (stdout,100) pnew(n),
     .        (resnc(n,igr)+2*mocogrp(igr),igr=1,ngroup)
        enddo
      endif
      write (stdout,106) suma, npnew, probcut, sumtot
c
c-----Computation of Mutual Entropy Information.
c
      if (doentropy) call 
     &  mutent (pnew,resnc,mocogrp,nprobt,np,ngroup,nel,largwr,stdout) 
c
      call ndelta (probt,p1a,p1b,resnca,resncb,
     &    ngroup,npra,nprb,stdout)
      if (allocated(p1a)) then
        deallocate (p1a,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p1a()'
      endif
      if (allocated(p1b)) then
        deallocate (p1b,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p1b()'
      endif
      if (allocated(probt)) then
        deallocate (probt,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate probt()'
      endif
!
!.....Compute spinless localization and delocalization indices
!
      if (.not.allocated(p1)) then
        allocate (p1(ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate p1()'
      endif
      if (.not.allocated(xli)) then
        allocate (xli(ngroup),stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot allocate xli()'
      endif
      p1=zero
      xli=zero
      allocate (npop(ngroup))
      if (ngroup.gt.1) then
        if (.not.allocated(p2)) then
          allocate (p2(ngroup*(ngroup-1)/2),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate p2()'
        endif
        if (.not.allocated(d1)) then
          allocate (d1(ngroup*(ngroup-1)/2),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate d1()'
        endif
        p2=zero
        d1=zero
      endif
      if (ngroup.gt.2) then
        if (.not.allocated(p3)) then
          allocate (p3(ngroup*(ngroup-1)*(ngroup-2)/6),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate p3()'
        endif
        if (.not.allocated(d2)) then
          allocate (d2(ngroup*(ngroup-1)*(ngroup-2)/6),stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot allocate d2()'
        endif
        p3=zero
        d2=zero
      endif
 
      nprob=np
      do np=1,nprob
        p=pnew(np)
c
c-------Add core electrons to each resonant structure
c
        npop(:)=resnc(np,:)+2*mocogrp(:)
        do i=1,ngroup
          p1(i)=p1(i)+npop(i)*p
          xli(i)=xli(i)+p*npop(i)*npop(i)
        enddo
        if (ngroup.gt.1) then
          ipair=0
          do i=1,ngroup
            do j=1,i-1
              ipair=ipair+1
              p2(ipair)=p2(ipair)+npop(i)*npop(j)*p
            enddo
          enddo
        endif
        if (ngroup.gt.2) then
          iter=0
          do i=1,ngroup
            do j=1,i-1
              do k=1,j-1
              iter=iter+1
              p3(iter)=p3(iter)+npop(i)*npop(j)*npop(k)*p
              enddo
            enddo
          enddo
        endif
      enddo
 
      write (stdout,11)
      do i=1,ngroup
        mcoi=mocogrp(i)
        write (stdout,15) i,p1(i)+2*mcoi
      enddo
c
c-----Compute variances
c
      call variance (ngroup,nprob,stdout,p1,pnew,resnc,mocogrp)

      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (stdout,2) i,j,p2(ipair)
          enddo
        enddo
      endif
 
      if (ngroup.gt.2) then
        iter=0
        do i=1,ngroup
          do j=1,i-1
            do k=1,j-1
            iter=iter+1
            write (stdout,3) i,j,k,p3(iter)
            enddo
          enddo
        enddo
      endif
 
      do np=1,nprob
        p=pnew(np)
c
c-------Add core electrons to each resonant structure
c
        npop(:)=resnc(np,:)+2*mocogrp(:)
        if (ngroup.gt.1) then
          ipair=0
          do i=1,ngroup
            do j=1,i-1
              ipair=ipair+1
              d1(ipair)=d1(ipair)-2*p*(p1(i)-npop(i))*(p1(j)-npop(j))
            enddo
          enddo
        endif
        if (ngroup.gt.2) then
          iter=0
          do i=1,ngroup
            do j=1,i-1
              do k=1,j-1
                iter=iter+1
                d2(iter)=d2(iter)-
     .          2*p*(p1(i)-npop(i))*(p1(j)-npop(j))*(p1(k)-npop(k))
              enddo
            enddo
          enddo
        endif
      enddo  
 
      do i=1,ngroup
        xli(i)=-(xli(i)-p1(i)*(p1(i)+one))
        write (stdout,6) i,i,xli(i),hundred*xli(i)/p1(i)
      enddo
      write (stdout,33)
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (stdout,4) i,j,d1(ipair)
          enddo
        enddo
      endif
      if (ngroup.gt.2) then
        iter=0
        do i=1,ngroup
          do j=1,i-1
            do k=1,j-1
            iter=iter+1
               write (stdout,5) i,j,k,d2(iter)
            enddo
          enddo
        enddo
      endif
c
c-----Diagonalize the covariance matrix
c
      call diagcov (ngroup,nprob,stdout,p1,xli,d1,pnew,resnc,mocogrp)

      if (allocated(p1)) then
        deallocate (p1,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p1()'
      endif
      if (allocated(xli)) then
        deallocate (xli,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate xli()'
      endif
      if (allocated(npop)) then
        deallocate (npop,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate npop()'
      endif
      if (allocated(resnc)) then
        deallocate (resnc,stat=ier)
        if (ier.ne.0) stop ' # binedf.f: Cannot deallocate resnc()'
      endif
      if (ngroup.gt.1) then
        if (allocated(p2)) then
          deallocate (p2,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p2()'
        endif
      endif
      if (ngroup.gt.2) then
        if (allocated(p3)) then
          deallocate (p3,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate p3()'
        endif
      endif
      if (ngroup.gt.1) then
        if (allocated(d1)) then
          deallocate (d1,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate d1()'
        endif
      endif
      if (ngroup.gt.2) then
        if (allocated(d2)) then
          deallocate (d2,stat=ier)
          if (ier.ne.0) stop ' # binedf.f: Cannot deallocate d2()'
        endif
      endif

      if (.not.mem) then
        close (indproba,status='delete')
        close (indprobb,status='delete')
      endif

      call timer (4,ibinedf,'_binedf   ',-1)
      return
!
!.....Formats.
!
 200  format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     . ' INCLUDING SPIN', /,1x, '# ',72('-'),/,
     . ' # NUMBER OF GROUPS               = ',I18,/,
     . ' # TOTAL NUMBER OF PROBABILITIES  = ',I18,/,
     . ' # Gi(a) Gi(b) ARE THE NUMBER OF ALPHA AND BETA ELECTRONS ',
     . ' IN GROUP i',/,1x,'# ',72('-'),/,
     . ' #     Probability',11x,'G1(a) G1(b) G2(a) G2(b) G3(a) G3(b)')
 203  format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION ',
     . 'FOR ',a,' ELECTRONS',/,
     . ' # FOR EACH VALUE, A SUM OVER ALL ',a,
     . ' RESONANT STRUCTURES HAS BEEN DONE',/,' #',72('-'),/,
     .  ' # NUMBER OF GROUPS',33x,' = ',I8,/,
     .  ' # TOTAL NUMBER OF PROBABILITIES FOR ',a,' ELECTRONS = ',I8,/,
     .   1x,'#',72('-'),/,
     . ' #     Probability            n1    n2    n3 ...')
 44   format (/,1x,a,/)
 444  format (1x,'# RECIPROCAL OF RCOND = ',E18.10,' FOR ',a,
     .        1x,'ELECTRONS')
 100  format (' # ',F22.16,1x,30I6)
 106  format (1x,'#',72('-'),/,' # ',F22.16,2x,'<-- SUM,',I8,
     . ' PROBABILITIES > ',E16.9,/,' # ',F22.16,2x,'<--- TOTAL SUM',/,
     . 1x,'#',72('-'))
 112  format (1x,'#',/,1x,'# multiple-group delocalization indices',/,
     .        1x,'#')
 113  format (1x,'# DELTA = ',1x,F16.9,5x,20(I3))
 202  format (/,' # M-BASINS SPINLESS ELECTRON DISTRIBUTION FUNCTION',
     . /,1x,'#',72('-'),/,
     .  ' # NUMBER OF GROUPS               = ',I8,/,
     .  ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     .   1x,'#',72('-'),/,
     . ' #     Probability            n1    n2    n3 ...')
 11   format (//1x,'Average populations and localization indices',
     &          1x,'(CORE ELECTRONS ARE INCLUDED)')
 33   format (//1x,'Delocalization indices,',
     .             ' Eq. (28) J. Chem. Phys.  126, 094102 (2007)')
 15   format (1x,'# <n(',I3,')>               = ',F16.9)
 2    format (1x,'# <n(',I3,') n(',I3,')>        = ',F16.9)
 3    format (1x,'# <n(',I3,') n(',I3,') n(',I3,')> = ',F16.9)
 4    format (1x,'# delta_(',I3,I3,')         = ',F16.9)
 5    format (1x,'# delta_(',I3,I3,I3,')      = ',F16.9)
 6    format (1x,'# delta_(',I3,I3,')         = ',F16.9,
     .        2x,'% Localization = ',F8.4)
 119  format (//1x,'Average populations and delocalization indices',
     &          1x,'(CORE ELECTRONS ARE INCLUDED)')
 19   format (1x,'<n(',I3,')_alpha>                  = ',F16.9,/,
     .        1x,'<n(',I3,')_beta>                   = ',F16.9)
 29   format (/1x,'<n(',I3,')_alpha n(',I3,')_alpha>     = ',F16.9,/,
     .         1x,'<n(',I3,')_alpha n(',I3,')_beta>      = ',F16.9,/,
     .         1x,'<n(',I3,')_beta  n(',I3,')_alpha>     = ',F16.9,/,
     .         1x,'<n(',I3,')_beta  n(',I3,')_beta>      = ',F16.9)
 49   format (/1x,'delta_(',I3,I3,')_{alpha,alpha)    = ',F16.9,/,
     .         1x,'delta_(',I3,I3,')_{alpha,beta )    = ',F16.9,/
     .         1x,'delta_(',I3,I3,')_{beta ,alpha)    = ',F16.9,/
     .         1x,'delta_(',I3,I3,')_{beta ,beta )    = ',F16.9)
 69   format (/1x,'delta_(',I3,I3,')_alpha            = ',F16.9,
     .         2x,'% Localization = ',F8.4,/,
     .         1x,'delta_(',I3,I3,')_beta             = ',F16.9,
     .         2x,'% Localization = ',F8.4)
 1120 format (3(1x,'#',/),1x,'# ',80('+'),/,
     .        1x,'# EXACT CALCULATION OF PROBABILITIES',/,1x,'#')
 1616 format (' # config ',a,I6,' : ',I6,' dets')
 1617 format (8I10)
 468  format (///" # binedf.f",/,
     . " # Fatal error: Too low value of MPAB, the first dimension",/,
     . " # of PACUM() and PBCUM() arrays. Increase the MPAB value",/,
     . " # in 'binedf.f' or use the NOMEM option in the input file")
 469  format (///" # binedf.f",/,
     . " # ========================================================="/,
     . " #  Warning: Too low value of MPAB, the first dimension of  "/,
     . " #  pacum() and pbcum() arrays. Intermediate files are used "/,
     . " #  instead in 'binedf.f'. If you prefer, please abort, put "/,
     . " #  MPAB = ",I8," in 'binedf.f and run again the program"/,
     . " # =========================================================")
 470  format (' # Solving linear system ',a,1x,i4,'% done')
 471  format (' # prob (alpha) x prob (beta)    ',I4,'% done')
 473  format (' # Codifying determinants  ',I4,'% done')
 474  format (' # Classifying configs  ',I4,'% done')
 475  format (' #  Ordering determinants  ',I4,'% done')
 2323 format (' # Config ',I6,' = ',64a1)
 1002 format (' #',/,' # ALPHA spin configs: Excitation Table',/,' #')
 1001 format (1000(' #',30I3,/))
 677  format (/,' # ',41('-'),
     &  /,' # Mutual information entropy between groups',/,' # ',
     & 41('-'))
 678  format (' # ',a,2x,20I4)
 680  format (' # Sum = ',F7.4)
 681  format (' # Pair (',2I3,' ) ',3x,'Prob(',2I3,' ) = ',F7.4,
     &  ' ?= ',F7.4,' = ',F7.4,' * ', F7.4)
 682  format (' # Entropy (',2I3,' ) = ',E15.8)
 683  format (1x,'# p(',I3,') = ',F7.4)
      end
