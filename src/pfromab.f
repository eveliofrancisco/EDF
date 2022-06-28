c
c-----This routine takes the ALPHA and BETA probabilities, stored
c     in 'filea' and 'fileb' respectively, of a molecule divided into
c     'ngroup' fragments and computes the SPINLESS probabilities.
c
c-----------------------------------------------------------------------
c
      subroutine pfromab (pcut,ngroup,lw,minp,maxp,filea,fileb)
      implicit none

      integer(kind=4), parameter :: nstack = 5000
      real   (kind=8) :: pcut,sumpcut
      real   (kind=8), allocatable,dimension (:)   :: pord
      integer(kind=4), allocatable,dimension (:)   :: iord
      integer(kind=4), allocatable,dimension (:,:) :: resncord
      integer(kind=4)    :: maxp(ngroup)
      integer(kind=4)    :: minp(ngroup)
      character(len=*)   :: filea,fileb
      integer(kind=4)    :: ngroup,ngroupx
      integer(kind=4)    :: lw
      character(len=200) :: line,word

      real   (kind=8), allocatable,dimension (:)  :: pnew
      integer(kind=4), allocatable,dimension(:,:) :: resnc
      integer(kind=4), allocatable,dimension(:,:) :: na,nb
      integer(kind=4), allocatable,dimension(:)   :: ind
      real   (kind=8), allocatable,dimension(:)   :: pa
      real   (kind=8), allocatable,dimension(:)   :: pb

      integer(kind=4)    :: i,npa,npb,np,m,j,ia,ib,k,n,nprob,ls,ios
      integer(kind=4)    :: napnb,ngreater
      logical            :: inlist,exfil

      ls  =  20
      exfil = .true.
      do while (exfil)
         inquire (unit=ls,opened=exfil)
         if (.not.exfil) then
           open (unit=ls,file=trim(filea),iostat=ios,status='old')
           if (ios.ne.0) then
             write (0,*) ' # pfromab.f: Error openning '//trim(filea)
             stop
           endif
         else
           ls=ls+1
         endif
      enddo
c
c-----Read ALPHA probabilities
c
      read (ls,*) ngroup,npa
      allocate (pa(npa))
      allocate (na(npa,ngroup))
      do j=1,npa
        read (ls,*) pa(j),(na(j,k),k=1,ngroup)
      enddo
      close (ls)

      exfil = .true.
      do while (exfil)
         inquire (unit=ls,opened=exfil)
         if (.not.exfil) then
           open (unit=ls,file=trim(fileb),iostat=ios,status='old')
           if (ios.ne.0) then
             write (0,*) ' # pfromab.f: Error openning '//trim(fileb)
             stop
           endif
         else
           ls=ls+1
         endif
      enddo
c
c-----Read BETA  probabilities
c
      read (ls,*) ngroupx,npb
      if (ngroup.ne.ngroupx) then
        stop '# pfromab.f: ALPHA/BETA probs different NGROUP values'
      endif
      allocate (pb(npb))
      allocate (nb(npb,ngroup))
      do j=1,npb
        read (ls,*) pb(j),(nb(j,k),k=1,ngroup)
      enddo
      close (ls)

      allocate (pnew(npa*npb))
      allocate (resnc(npa*npb,ngroup))
      allocate (ind(ngroup))
      pnew=0d0
      n=0
      do ia=1,npa
        cycleib: do ib=1,npb
c
c---------Skip this probability if minp(i) > elec(i) or maxp(i) < elec(i)
c
          do k=1,ngroup
            napnb=na(ia,k)+nb(ib,k)
            if (napnb.gt.maxp(k).or.napnb.lt.minp(k)) cycle cycleib
          enddo
          n=n+1
          if (n.eq.1) then
            np=1
            resnc(np,1:ngroup)=na(ia,1:ngroup)+nb(ib,1:ngroup)
            pnew(np)=pa(ia)*pb(ib)
          else
            ind(1:ngroup)=na(ia,1:ngroup)+nb(ib,1:ngroup)
            do j=1,np
              inlist=.true.
              do k=1,ngroup
                inlist=inlist.and.(ind(k).eq.resnc(j,k))
              enddo
              if (inlist) then
                pnew(j)=pnew(j)+pa(ia)*pb(ib)
                cycle cycleib
              endif
            enddo
            np=np+1
            pnew(np)=pa(ia)*pb(ib)
            resnc(np,1:ngroup)=na(ia,1:ngroup)+nb(ib,1:ngroup)
          endif
        enddo cycleib
      enddo
*     do i=1,np
*       write (6,'(F22.15,30I4)') pnew(i),(resnc(i,j),j=1,ngroup)
*     enddo

      allocate (pord(np))
      allocate (iord(np))
      allocate (resncord(np,ngroup))
      pord(1:np)=pnew(1:np)
      forall (i=1:np) iord(i)=i
      nprob=min(np,nstack)
      call qqsort (pord,iord,1,nprob,nprob)
      write (lw,20) nprob,pcut
      sumpcut=0d0
      ngreater=0
      do n=nprob,1,-1
        m=iord(n)
        if (pord(m).ge.pcut) then
          sumpcut=sumpcut+pord(m)
          ngreater=ngreater+1
          write (lw,10) pord(m),(resnc(m,j),j=1,ngroup)
        endif
      enddo
      write (lw,3) sum(pord(1:nprob)),sumpcut,pcut
      write (lw,4) ngreater
 10   format (1x,'# ',F22.15,3x,20I6)
 20   format (//,' #',6x,'THERE ARE ',I8,' Spinless probabilities',/,
     & ' #',6x,'Printing probabilities greater than ',F15.8,/,
     & ' #',6x,'probability             RSRS --->')
 3    format (1x,'# ',F22.15,2x,'<--- SUM',/,
     &        1x,'# ',F22.15,2x,'<--- SUM ( p >',1x,F15.8,' )')
 4    format (1x,'# ',I10,' probabilities printed')
      end
