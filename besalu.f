
      integer(kind=4), allocatable,dimension(:) :: ilo(ngroup)
      integer(kind=4), allocatable,dimension(:) :: jlo(ngroup)
      integer(kind=4), allocatable,dimension(:) :: flo(ngroup)
      integer(kind=4), allocatable,dimension(:) :: slo(ngroup)




      allocate (ilo(nelectrons))
      allocate (jlo(nelectrons))
      allocate (flo(nelectrons))
      allocate (slo(nelectrons))
      allocate (rsrs(nregions))


      allocate (region(10))
      region( 1)='A'
      region( 2)='B'
      region( 3)='C'
      region( 4)='D'
      region( 5)='E'
      region( 6)='F'
      region( 7)='G'
      region( 9)='H'
      region(10)='I'
      region(10)='J'

      do k=1,nelectrons
        ilo(k)=1        ! Initial value
        flo(k)=nregions ! Final value
        slo(k)=1        ! Step increment
        jlo(k)=ilo(k)   ! Initialization of NSS variable
      enddo
      k=nelectrons
      ic = 0
      do while (k.gt.0)
         if ((jlo(k)-flo(k))*slo(k).gt.0) then
            jlo(k)=ilo(k)
            k=k-1
         else
            ic = ic + 1
            write (6,1,advance='no') (m,region(jlo(m)),m=1,nelectrons)
            rsrs(1:nregions)=0
            do m=1,nelectrons
              mm=jlo(m)
              rsrs(mm)=rsrs(mm)+1
            enddo
            write (6,2) (rsrs(m),m=1,nregions)
            k=nelectrons
         end if
         if (k.gt.0) jlo(k)=jlo(k)+slo(k)
      end do
      write (6,*) 'terms =',ic
 1    format (1x,1000(' (El',I2,' in ',A1,')',3x))
 2    format (5x,'RSRS = ',20I3)
      end





c
c-----------------------------------------------------------------------
c
      subroutine sieteprsrs (epsdet,ndets,okord,ngroup,nta,ntb,nmo,
     &  ncore,nel,lw,lr,maxp,minp,probal,line)
      USE space_for_sgarea
      include     'implicit.inc'
      include     'fact.inc'
      include     'constants.inc'
      parameter (maxg=10)
      integer occa(ngroup),occb(ngroup),occup(ngroup)
      integer popul(ngroup),maxp(ngroup),minp(ngroup)
      real(kind=8)    epsdet,probal,pab
      character*(*) line
      logical setint,ok,okord,docalc,calcres
c
      okord=.true.
c
      probal=zero
      if (nta+ntb.ne.nel) return
      if (ngroup.gt.maxg) then
        write (lr,100) maxg
        okord=.false.
        return
      endif
c
c.....Built in the set of integers that define the RSRS.
c
      lp=1
      npop=0
      id=0
      nugr=1
      do while (nugr.le.ngroup)
        ok = setint (ipol,line,lp)
        if (ok) then
          id=id+1
          if (ipol.le.nel.and.ipol.ge.0) then
            npop=npop+ipol
            if (npop.gt.nel) then
              write (lr,101)
              okord=.false.
              return
            else
              popul(id)=ipol
            endif
            nugr=nugr+1
          elseif (ipol.gt.nel) then
            write (lr,102) id
            okord=.false.
            return
          else
            write (lr,103) id
            okord=.false.
            return
          endif
        else
          write (lr,104)
          okord=.false.
          return
        endif
      enddo
      occup(1:ngroup)=popul(1:ngroup)

      write (lw,10)
      write (lw,11) ngroup,nta,ntb,(occup(i),i=1,ngroup)
      write (lw,111) (minp(i),i=1,ngroup)
      write (lw,112) (maxp(i),i=1,ngroup)
c
      ntis=0
      do i=1,ngroup
        ntis=ntis+popul(i)
      enddo
      if (ntis.ne.nel) then
         write (lr,105) 
         okord=.false.
         return
      endif

      probal=zero
      nspinres=0
      if (ngroup.eq.1) then
        occa(1)=nta
        docalc=.true.
        do i=1,ngroup
           if (occa(i).gt.maxp(i).or.
     &         occa(i).lt.minp(i)) docalc=.false.
        enddo
        if (docalc) then
          call wrsrs (epsdet,ndets,nta,ntb,nmo,
     &         ncore,ngroup,lw,pab,occup,occa,occb)
          probal=probal+pab
          if (ndets.gt.1) then
            call wriprob (pab,occa,occup,ngroup,lw)
          endif
        else
          write (lw,6,advance='no') (occa(i),i=1,ngroup)
          write (lw,7) (occup(i)-occa(i),i=1,ngroup)
        endif
        nspinres=nspinres+1
      elseif (ngroup.eq.2) then
        ia1=0
        do while (ia1.le.occup(1))
          occa(1)=ia1
          ia2=0
          ia12=ia1+ia2
          do while (ia2.le.occup(2).and.ia12.le.nta)
            occa(2)=ia2
            if (ia12.eq.nta) then
              docalc=.true.
              do i=1,ngroup
                if (occa(i).gt.maxp(i).or.
     &              occa(i).lt.minp(i)) docalc=.false.
              enddo
              if (docalc) then
                call wrsrs (epsdet,ndets,nta,ntb,nmo,
     &               ncore,ngroup,lw,pab,occup,occa,occb)
                probal=probal+pab
                if (ndets.gt.1) then
                  call wriprob (pab,occa,occup,ngroup,lw)
                endif
              else
                write (lw,6,advance='no') (occa(i),i=1,ngroup)
                write (lw,7) (occup(i)-occa(i),i=1,ngroup)
              endif
              nspinres=nspinres+1
            endif
            ia2=ia2+1
            ia12=ia1+ia2
          enddo  
          ia1=ia1+1
        enddo
      elseif (ngroup.eq.3) then
        ia1=0
        do while (ia1.le.occup(1))
          occa(1)=ia1
          ia2=0
          ia12=ia1+ia2
          do while (ia2.le.occup(2).and.ia12.le.nta)
            occa(2)=ia2
            ia3=0
            ia123=ia12+ia3
            do while (ia3.le.occup(3).and.ia123.le.nta)
              occa(3)=ia3
              if (ia123.eq.nta) then
                docalc=.true.
                do i=1,ngroup
                  if (occa(i).gt.maxp(i).or.
     &                occa(i).lt.minp(i)) docalc=.false.
                enddo
                if (docalc) then
                  call wrsrs (epsdet,ndets,nta,ntb,nmo,
     &                 ncore,ngroup,lw,pab,occup,occa,occb)
                  probal=probal+pab
                  if (ndets.gt.1) then
                    call wriprob (pab,occa,occup,ngroup,lw)
                  endif
                else
                  write (lw,6,advance='no') (occa(i),i=1,ngroup)
                  write (lw,7) (occup(i)-occa(i),i=1,ngroup)
                endif
                nspinres=nspinres+1
              endif
              ia3=ia3+1
              ia123=ia12+ia3
            enddo
            ia2=ia2+1
            ia12=ia1+ia2
          enddo  
          ia1=ia1+1
        enddo
      elseif (ngroup.eq.4) then
        do ia1=0,occup(1)
          occa(1)=ia1
          do ia2=0,occup(2)
            occa(2)=ia2
            do ia3=0,occup(3)
              occa(3)=ia3
              do ia4=0,occup(4)
                occa(4)=ia4
c-----------------------------------------------------------------------
                if (sum(occa(1:ngroup)).ne.nta) cycle
                docalc=calcres(ngroup,occa,minp,maxp)
                if (docalc) then
                  call wrsrs (epsdet,ndets,nta,ntb,nmo,
     &                ncore,ngroup,lw,pab,occup,occa,occb)
                  probal=probal+pab
                  if (ndets.gt.1) then
                    call wriprob (pab,occa,occup,ngroup,lw)
                  endif
                else
                  write (lw,6,advance='no') (occa(i),i=1,ngroup)
                  write (lw,7) (occup(i)-occa(i),i=1,ngroup)
                endif
c-----------------------------------------------------------------------
              enddo
            enddo
          enddo
        enddo
      elseif (ngroup.eq.5) then
        ia1=0
        do while (ia1.le.occup(1))
          occa(1)=ia1
          ia2=0
          ia12=ia1+ia2
          do while (ia2.le.occup(2).and.ia12.le.nta)
            occa(2)=ia2
            ia3=0
            ia123=ia12+ia3
            do while (ia3.le.occup(3).and.ia123.le.nta)
              occa(3)=ia3
              ia4=0
              ia1234=ia123+ia4
              do while (ia4.le.occup(4).and.ia1234.le.nta)
                occa(4)=ia4
                ia5=0
                ia12345=ia1234+ia5
                do while (ia5.le.occup(5).and.ia12345.le.nta)
                   occa(5)=ia5
                   if (ia12345.eq.nta) then
                     docalc=.true.
                     do i=1,ngroup
                       if (occa(i).gt.maxp(i).or.
     &                     occa(i).lt.minp(i)) docalc=.false.
                     enddo
                     if (docalc) then
                       call wrsrs (epsdet,ndets,nta,ntb,nmo,
     &                      ncore,ngroup,lw,pab,occup,occa,occb)
                       probal=probal+pab
                       if (ndets.gt.1) then
                         call wriprob (pab,occa,occup,ngroup,lw)
                       endif
                     else
                       write (lw,6,advance='no') (occa(i),i=1,ngroup)
                       write (lw,7) (occup(i)-occa(i),i=1,ngroup)
                     endif
                     nspinres=nspinres+1
                   endif
                   ia5=ia5+1
                   ia12345=ia1234+ia5
                enddo
                ia4=ia4+1
                ia1234=ia123+ia4
              enddo
              ia3=ia3+1
              ia123=ia12+ia3
            enddo
            ia2=ia2+1
            ia12=ia1+ia2
          enddo  
          ia1=ia1+1
        enddo
      elseif (ngroup.eq.6) then
        ia1=0
        do while (ia1.le.occup(1))
          occa(1)=ia1
          ia2=0
          ia12=ia1+ia2
          do while (ia2.le.occup(2).and.ia12.le.nta)
            occa(2)=ia2
            ia3=0
            ia123=ia12+ia3
            do while (ia3.le.occup(3).and.ia123.le.nta)
              occa(3)=ia3
              ia4=0
              ia1234=ia123+ia4
              do while (ia4.le.occup(4).and.ia1234.le.nta)
                occa(4)=ia4
                ia5=0
                ia12345=ia1234+ia5
                do while (ia5.le.occup(5).and.ia12345.le.nta)
                   occa(5)=ia5
                   ia6=0
                   ia123456=ia12345+ia6
                   do while (ia6.le.occup(6).and.ia123456.le.nta)
                     occa(6)=ia6
                     if (ia123456.eq.nta) then
                       docalc=.true.
                       do i=1,ngroup
                         if (occa(i).gt.maxp(i).or.
     &                       occa(i).lt.minp(i)) docalc=.false.
                       enddo
                       if (docalc) then
                         call wrsrs (epsdet,ndets,nta,ntb,nmo,
     &                        ncore,ngroup,lw,pab,occup,occa,occb)
                         probal=probal+pab
                         if (ndets.gt.1) then
                           call wriprob (pab,occa,occup,ngroup,lw)
                         endif
                       else
                         write (lw,6,advance='no') (occa(i),i=1,ngroup)
                         write (lw,7) (occup(i)-occa(i),i=1,ngroup)
                       endif
                       nspinres=nspinres+1
                     endif
                     ia6=ia6+1
                     ia123456=ia12345+ia6
                   enddo
                   ia5=ia5+1
                   ia12345=ia1234+ia5
                enddo
                ia4=ia4+1
                ia1234=ia123+ia4
              enddo
              ia3=ia3+1
              ia123=ia12+ia3
            enddo
            ia2=ia2+1
            ia12=ia1+ia2
          enddo  
          ia1=ia1+1
        enddo
      elseif (ngroup.eq.7) then
        ia1=0
        do while (ia1.le.occup(1))
          occa(1)=ia1
          ia2=0
          ia12=ia1+ia2
          do while (ia2.le.occup(2).and.ia12.le.nta)
            occa(2)=ia2
            ia3=0
            ia123=ia12+ia3
            do while (ia3.le.occup(3).and.ia123.le.nta)
              occa(3)=ia3
              ia4=0
              ia1234=ia123+ia4
              do while (ia4.le.occup(4).and.ia1234.le.nta)
                occa(4)=ia4
                ia5=0
                ia12345=ia1234+ia5
                do while (ia5.le.occup(5).and.ia12345.le.nta)
                   occa(5)=ia5
                   ia6=0
                   ia123456=ia12345+ia6
                   do while (ia6.le.occup(6).and.ia123456.le.nta)
                     occa(6)=ia6
                     ia7=0
                     ia1234567=ia123456+ia7
                     do while (ia7.le.occup(7).and.ia1234567.le.nta)
                       occa(7)=ia7
                       if (ia1234567.eq.nta) then
                         docalc=.true.
                         do i=1,ngroup
                           if (occa(i).gt.maxp(i).or.
     &                         occa(i).lt.minp(i)) docalc=.false.
                         enddo
                         if (docalc) then
                           call wrsrs (epsdet,ndets,nta,ntb,nmo,
     &                          ncore,ngroup,lw,pab,occup,occa,occb)
                           probal=probal+pab
                           if (ndets.gt.1) then
                             call wriprob (pab,occa,occup,ngroup,lw)
                           endif
                         else
                           write (lw,6,advance='no')(occa(i),i=1,ngroup)
                           write (lw,7) (occup(i)-occa(i),i=1,ngroup)
                         endif
                         nspinres=nspinres+1
                       endif
                       ia7=ia7+1
                       ia1234567=ia123456+ia7
                     enddo
                     ia6=ia6+1
                     ia123456=ia12345+ia6
                   enddo
                   ia5=ia5+1
                   ia12345=ia1234+ia5
                enddo
                ia4=ia4+1
                ia1234=ia123+ia4
              enddo
              ia3=ia3+1
              ia123=ia12+ia3
            enddo
            ia2=ia2+1
            ia12=ia1+ia2
          enddo  
          ia1=ia1+1
        enddo
      else
        write (lr,100) 6
        okord=.false.
        return
      endif
      write (lw,*) '# ',nspinres,' spin-resolved probabilities'
      return
c
 6    format (' #  Ignoring config (alpha',100I3)
 7    format ('   --- beta',100I3)
 10   format (//1x,'# ',88('-'),/,' #',10x,
     & 'Spin  resolved components of a spin unresolved probability',
     & /1x,'# ',88('-'))
 11   format (1x,'# Number of domains                  = ',I3,/,
     &        1x,'# Number of alpha and beta electrons = ',2I3,/,
     &        1x,'# Real space resonance structure     = ',20I3,
     &        /,' #')
 111  format (1x,'# Mininum population allowed         = ',20I3)
 112  format (1x,'# Maximum population allowed         = ',20I3)
 100  format (' #',/,' # prsrs.f: Too many groups, Maximum value is',I4)
 101  format (' #',/,' # prsrs.f: Wrong number of electrons,'
     &        ' N_alpha + N_beta .ne. total number of electrons')
 102  format (' #',/,' # prsrs.f: Too high population in domain ',I2)
 103  format (' #',/,' # prsrs.f: Negative population in domain ',I1)
 104  format (' #',/,' # prsrs.f: Erroneus PRSRS order in input file')
 105  format (' #',/,' # prsrs.f: Improper number of electrons in RSRS')
      end
c
c-----------------------------------------------------------------------
c
      subroutine wriprob (pab,occa,occup,ngroup,lw)
      real(kind=8) pab
      integer occa,occup,ngroup,lw
      dimension occa(ngroup),occup(ngroup)
c
      write (lw,400,advance='no') 
 400  format (' # Computing config (alpha')
      write (lw,401,advance='no') (occa(i),i=1,ngroup)
      write (lw,403,advance='no') (occup(i)-occa(i),i=1,ngroup)
 403  format ('   --- beta',100I3)
 401  format (100I3)
      write (lw,402) pab
 402  format (' ---> PROB = ',F16.10,1x)
      return

      write (lw,40,advance='no') pab
      write (lw,20,advance='no')(occa(i),i=1,ngroup)
      write (lw,'(a)',advance='no') '  --- '
      write (lw,30,advance='no')(occup(i)-occa(i),i=1,ngroup)
      write (lw,31) 
      return
c
 20   format (' for alpha ',20I3) 
 30   format ('beta ',20I3,' resonant structure')
 31   format (' resonant structure')
 40   format (1x,'# PROB = ',F16.10,1x)
      end
c
c-----------------------------------------------------------------------
c
      subroutine wrsrs 
     &  (epsdet,ndets,nta,ntb,nmo,ncore,ngroup,lw,prob,n,na,nb)
      USE space_for_sgarea
      include     'fact.inc'
      include     'constants.inc'
      parameter (maxg=10)
      real(kind=8) epsdet,prob,tottimes

      integer, allocatable,dimension (:)    :: iacta,iactb,iactt
      integer, allocatable,dimension (:)    :: in1,in2,in3,in4,in5,in6
      integer, allocatable,dimension (:)    :: ou1,ou2,ou3,ou4,ou5,ou6
      integer, allocatable,dimension (:)    :: io1,io2,io3,io4,io5,io6
      integer, allocatable,dimension (:)    :: wc
      integer, allocatable,dimension (:,:)  :: wca,wcb


      integer ngroup,lw,n(ngroup),na(ngroup),nb(ngroup)
      integer da(maxg),db(maxg),dd(maxg)
      integer*8 ntimes(maxg),nttim(2)
      logical first(maxg)
 
      if (ndets.eq.1) then
        call sdwprsrs (nta,ntb,nmo,ngroup,lw,prob,n,na,nb)
        return
      endif
      nel=nta+ntb
c
      allocate (iacta(nel))
      allocate (iactb(nel))
      allocate (iactt(nel))
      allocate (wc(nmo))
c
      allocate (in1(nel),ou1(nel),io1(nel))
      if (ngroup.gt.1) allocate (in2(nel),ou2(nel),io2(nel))
      if (ngroup.gt.2) allocate (in3(nel),ou3(nel),io3(nel))
      if (ngroup.gt.3) allocate (in4(nel),ou4(nel),io4(nel))
      if (ngroup.gt.4) allocate (in5(nel),ou5(nel),io5(nel))
      if (ngroup.gt.5) allocate (in6(nel),ou6(nel),io6(nel))

      nacta=0
      nactb=0
      do i=1,ngroup
        nb(i)=n(i)-na(i)
        if (na(i).gt.0) then
          nacta=nacta+1
          iacta(nacta)=na(i)
          da(nacta)=i
        endif
        if (nb(i).gt.0) then
          nactb=nactb+1
          iactb(nactb)=nb(i)
          db(nactb)=i
        endif
      enddo
      do isp=1,2
        if (isp.eq.1) then
          nactt=nacta
          ntab=nta
          iactt(1:nactt)=iacta(1:nacta)
          dd(1:nactt)=da(1:nacta)
        else
          nactt=nactb
          ntab=ntb
          iactt(1:nactt)=iactb(1:nactb)
          dd(1:nactt)=db(1:nactb)
        endif
        nrest=0
        tottimes=one
        if (ntab.gt.maxfac) then
          write (0,*) 'NTAB,MAXFAC = ',ntab,maxfac
          stop 'prsrs.f: !!! NTAB > MAXFAC. Increase MAXFAC in fact.inc'
        endif
        do i=1,nactt-1
          num=ntab-nrest
          n1=iactt(i)
          ntimes(i)=nint(fact(num)/fact(n1)/fact(num-n1))
          tottimes=tottimes*ntimes(i)
          nrest=nrest+n1
        enddo

        if (tottimes.gt.huge(l)) then
          write (0,123) tottimes
          write (lw,123) tottimes
          stop 
        endif
 123    format (3(' #',/),
     &    ' # prsrs.f: I am ver sorry.',/,
     &    ' # prsrs.f: A very HUGE integer number should be needed'/,
     &    ' # prsrs.f: nttim() array should need ',F22.0,' elements'/
     &    ' # --- STOP --- STOP --- STOP --- STOP --- STOP --- STOP',/,
     &   3(' #',/))
        nttim(isp)=idnint(tottimes)
        if (isp.eq.1) then
          allocate (wca(nttim(isp),nmo))
        else
          allocate (wcb(nttim(isp),nmo))
        endif

        first(1:ngroup)=.true.
        if (nactt.eq.1) then
          ia1=iactt(1)
          do k=1,ia1
            wc(k)=dd(1)
          enddo
          ncom=1
          if (isp.eq.1) then
            wca(ncom,1:ntab)=wc(1:ntab)
          else
            wcb(ncom,1:ntab)=wc(1:ntab)
          endif
        elseif (nactt.eq.2) then
          first(1)=.true.
          ntt=ntab
          do i=1,ntab
            io1(i)=i
          enddo
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            ncom=ncom+1
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            do k=1,ia1
              wc(in1(k))=dd(1)
            enddo
            ia2=iactt(2)
            do k=1,ia2
              wc(ou1(k))=dd(2)
            enddo
            if (isp.eq.1) then
              wca(ncom,1:ntab)=wc(1:ntab)
            else
              wcb(ncom,1:ntab)=wc(1:ntab)
            endif
          enddo
        elseif (nactt.eq.3) then
          first(1)=.true.
          ntt=ntab
          do i=1,ntab
            io1(i)=i
          enddo
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            do k=1,ia1
              wc(in1(k))=dd(1)
            enddo
            first(2)=.true.
            ntt1=ntab-ia1
            io2(1:ntt1)=ou1(1:ntt1)
            ia2=iactt(2)
            do i2=1,ntimes(2)
              ncom=ncom+1
              call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
              do k=1,ia2
                wc(in2(k))=dd(2)
              enddo
              ia3=iactt(3)
              do k=1,ia3
                wc(ou2(k))=dd(3)
              enddo
              if (isp.eq.1) then
                wca(ncom,1:ntab)=wc(1:ntab)
              else
                wcb(ncom,1:ntab)=wc(1:ntab)
              endif
            enddo
          enddo
        elseif (nactt.eq.4) then
          first(1)=.true.
          ntt=ntab
          do i=1,ntab
            io1(i)=i
          enddo
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            do k=1,ia1
              wc(in1(k))=dd(1)
            enddo
            first(2)=.true.
            ntt1=ntab-ia1
            io2(1:ntt1)=ou1(1:ntt1)
            ia2=iactt(2)
            do i2=1,ntimes(2)
              call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
              do k=1,ia2
                wc(in2(k))=dd(2)
              enddo
              first(3)=.true.
              ntt2=ntab-ia1-ia2
              io3(1:ntt2)=ou2(1:ntt2)
              ia3=iactt(3)
              do i3=1,ntimes(3)
                ncom=ncom+1
                call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
                do k=1,ia3
                  wc(in3(k))=dd(3)
                enddo
                ia4=iactt(4)
                do k=1,ia4
                  wc(ou3(k))=dd(4)
                enddo
                if (isp.eq.1) then
                  wca(ncom,1:ntab)=wc(1:ntab)
                else
                  wcb(ncom,1:ntab)=wc(1:ntab)
                endif
              enddo
            enddo
          enddo
        elseif (nactt.eq.5) then
          first(1)=.true.
          ntt=ntab
          do i=1,ntab
            io1(i)=i
          enddo
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            do k=1,ia1
              wc(in1(k))=dd(1)
            enddo
            first(2)=.true.
            ntt1=ntab-ia1
            io2(1:ntt1)=ou1(1:ntt1)
            ia2=iactt(2)
            do i2=1,ntimes(2)
              call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
              do k=1,ia2
                wc(in2(k))=dd(2)
              enddo
              first(3)=.true.
              ntt2=ntab-ia1-ia2
              io3(1:ntt2)=ou2(1:ntt2)
              ia3=iactt(3)
              do i3=1,ntimes(3)
                call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
                do k=1,ia3
                  wc(in3(k))=dd(3)
                enddo
                first(4)=.true.
                ntt3=ntab-ia1-ia2-ia3
                io4(1:ntt3)=ou3(1:ntt3)
                ia4=iactt(4)
                do i4=1,ntimes(4)
                  ncom=ncom+1
                  call gengen(4,ia4,ntt3,maxg,nel,io4,in4,ou4,first)
                  do k=1,ia4
                    wc(in4(k))=dd(4)
                  enddo
                  ia5=iactt(5)
                  do k=1,ia5
                    wc(ou4(k))=dd(5)
                  enddo
                  if (isp.eq.1) then
                    wca(ncom,1:ntab)=wc(1:ntab)
                  else
                    wcb(ncom,1:ntab)=wc(1:ntab)
                  endif
                enddo
              enddo
            enddo
          enddo
        elseif (nactt.eq.6) then
          first(1)=.true.
          ntt=ntab
          do i=1,ntab
            io1(i)=i
          enddo
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            do k=1,ia1
              wc(in1(k))=dd(1)
            enddo
            first(2)=.true.
            ntt1=ntab-ia1
            io2(1:ntt1)=ou1(1:ntt1)
            ia2=iactt(2)
            do i2=1,ntimes(2)
              call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
              do k=1,ia2
                wc(in2(k))=dd(2)
              enddo
              first(3)=.true.
              ntt2=ntab-ia1-ia2
              io3(1:ntt2)=ou2(1:ntt2)
              ia3=iactt(3)
              do i3=1,ntimes(3)
                call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
                do k=1,ia3
                  wc(in3(k))=dd(3)
                enddo
                first(4)=.true.
                ntt3=ntab-ia1-ia2-ia3
                io4(1:ntt3)=ou3(1:ntt3)
                ia4=iactt(4)
                do i4=1,ntimes(4)
                  call gengen(4,ia4,ntt3,maxg,nel,io4,in4,ou4,first)
                  do k=1,ia4
                    wc(in4(k))=dd(4)
                  enddo
                  first(5)=.true.
                  ntt4=ntab-ia1-ia2-ia3-ia4
                  io5(1:ntt4)=ou4(1:ntt4)
                  ia5=iactt(5)
                  do i5=1,ntimes(5)
                    ncom=ncom+1
                    call gengen(5,ia5,ntt4,maxg,nel,io5,in5,ou5,first)
                    do k=1,ia5
                      wc(in5(k))=dd(5)
                    enddo
                    ia6=iactt(6)
                    do k=1,ia6
                      wc(ou5(k))=dd(6)
                    enddo
                    if (isp.eq.1) then
                      wca(ncom,1:ntab)=wc(1:ntab)
                    else
                      wcb(ncom,1:ntab)=wc(1:ntab)
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        else
          stop 'prsrs.f: Too many domains'
        endif
      enddo
      call sumdetp (epsdet,prob,nttim(1),nttim(2),
     &   nta,ntb,ncore,wca,wcb,ngroup,ndets,nmo,.false.)

      deallocate (iacta,iactb,iactt,wca,wcb,wc)
      deallocate (in1,ou1,io1)
      if (ngroup.gt.1) deallocate (in2,ou2,io2)
      if (ngroup.gt.2) deallocate (in3,ou3,io3)
      if (ngroup.gt.3) deallocate (in4,ou4,io4)
      if (ngroup.gt.4) deallocate (in5,ou5,io5)
      if (ngroup.gt.5) deallocate (in6,ou6,io6)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine gengen (ig,n,m,maxg,nel,iord,inside,outside,first)
c
c.....Extracts a subset of N elements from a set of M>N natural numbers
c     iord(1), ..., iord(M), non necessarily consecutive storing them in 
c     the set inside[1..N]. The remaining M-N elements are stored in the 
c     set outside[1..M-N]. Each time the routine is called, it produces 
c     new subsets of 'N' and 'M-N' elements. Before 'gengen' is called
c     for the first time the logical variable FIRST(IG) must be set to 
c     .TRUE. This first time the returned inside[1..N] and outside[1..M-N] 
c     sets will contain the first N and the last M-N numbers of the set. 
c     This first time, 'gengen' makes FIRST(IG) equal to .FALSE. One 
c     must be careful that the routine is called as much (M over N) times.
c
      parameter (maxinout=100)
      integer combin(maxinout,maxg),combout(maxinout,maxg)
      integer n,m,i,j,isin,isout,ig
      integer inside(nel),outside(nel)
      integer iord(nel)
      logical first(maxg)
c
      if (n.gt.maxinout.or.m-n.gt.maxinout) then
        stop 'gengen.f: Increase MAXINOUT parameter in gengen.f'
      endif
      if (first(ig)) then
        first(ig)=.false.
        do i=1,n
          combin(i,ig)=i
          inside(i)=iord(i)
        enddo
        do i=n+1,m
          in=i-n
          combout(in,ig)=i
          outside(in)=iord(i)
        enddo
      else
        i=n
 1      do while (i.ge.1)
          if (combin(i,ig).lt.m) then
            combin(i,ig)=combin(i,ig)+1
            do j=i+1,n
              combin(j,ig)=combin(j-1,ig)+1
              if (combin(j,ig).gt.m) then
                 i=i-1
                 goto 1
              endif
            enddo
            isout=0
            do j=1,m
              do i=1,n
                if (combin(i,ig).eq.j) goto 2
              enddo
              isout=isout+1
              combout(isout,ig)=j
 2            continue
            enddo
            inside (1:n)  =iord(combin(1:n,ig))
            outside(1:m-n)=iord(combout(1:m-n,ig))
            return
          else
            i=i-1
          endif
        enddo
      endif
      return
      end
