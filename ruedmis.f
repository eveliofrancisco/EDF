c
c-----------------------------------------------------------------------
c
      subroutine ruedmis (sg,c,natp,nmo,ngroup,lw,largwr,okw)
c
c.....ruedmis - Performs the Cioslowski isopynic localization of a set 
c     of natural molecular orbitals into ngroup space exhausting domains
c     by maximizing 
c                            
c     L = sum [i=1,NMO] nu_i^2 L_i
c                                   
c     L_i = SUM [klmn] A_klmn
c
c     A_klmn = U^*(ik) U(il) U^*(im) U(in) N(klmn) T(klmn)
c
c     N(klmn) = (n_k n_l n_m n_n)^(1/2)
c
c     T(klmn)=SUM_A <NMO_k|NMO_l>_A <NMO_m|NMO_n>_A
c
c.....where
c
c     NMO ....... Number of natural molecular orbitals (NMO)
c     U(rs) ....  The elements of an unitary matrix
c     A ......... Index that runs over all the fragments
c     NMO_r ..... The r^th molecular orbital (MO)
c     n_r ....... Natural population of r^th MO.
c     nu_i ...... Natural population of the i^th isopycnic MO
c                            
c-----PARAMETERS--------------------------------------------------------
c
c.....sg()                           INPUT/OUTPUT
c
c     Group overlap matrix between input MO's when entering the 
c     routine or between output MO's at the end of the routine.
c
c.....c()                            OUTPUT 
c
c     Rotation matrix giving the transformed MO's from the input MO's. 
c
c.....natp()                         INPUT/OUTPUT
c
c     Natural orbitals populations/Transformed MOs populations
c
c.....nmo                            INPUT
c
c     Number of MO's.
c
c.....ngroup                         INPUT
c
c     Number of fragments 
c
c     Maximum number of groups
c
c.....lw                             INPUT
c
c     Output logical unit
c
c.....largwr                         INPUT
c
c     .true.  ====>  Large output
c     .false. ====>  Short output
c
c.....okw                            INPUT
c
c     .false. ====> The routine writes nothing.
c     .true.  ====> Output according to the value of the LARGWR value
c-----------------------------------------------------------------------
c                         
      implicit real(kind=8) (a-h,o-z)
      include 'constants.inc'
c
      parameter (maxit=200)
      parameter (xmoloc=0.05d0)

      real(kind=8)   c(nmo,nmo)
      real(kind=8)   pi,pi4
      real(kind=8)   natp(nmo)
      real(kind=8)   sg(ngroup,nmo,nmo)

      real(kind=8),allocatable,dimension (:)     :: xloci,xion,degloc
      real(kind=8),allocatable,dimension (:,:)   :: over
      real(kind=8),allocatable,dimension (:,:,:) :: sgtmp
      real(kind=8),allocatable,dimension (:,:)   :: ctmp
      real(kind=8),allocatable,dimension (:)     :: delta,deltad,rotinc
      integer, allocatable,dimension (:)      :: indi,ixion

      logical   warn,largwr,okw,goon
      logical   moreiter
      integer   ipair(2)
c
      call timer (2,ipid,'_ruedmis  ',-1)
c
c.....Allocate arrays
c
      n = nmo
      if (.not.allocated(over)) then
        allocate (over(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate over()'
      endif

      if (.not.allocated(sgtmp)) then
        allocate (sgtmp(ngroup,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate sgtmp()'
      endif

      if (.not.allocated(xloci)) then
        allocate (xloci(nmo),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate xloci()'
      endif

      if (.not.allocated(xion)) then
        allocate (xion(nmo),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate xion()'
      endif
      if (.not.allocated(degloc)) then
        allocate (degloc(nmo),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate degloc()'
      endif

      if (.not.allocated(ixion)) then
        allocate (ixion(nmo),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate ixion()'
      endif

      if (.not.allocated(ctmp)) then
        allocate (ctmp(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate ctmp()'
      endif

      if (.not.allocated(delta)) then
        allocate (delta(ngroup*(ngroup-1)/2),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate delta()'
      endif

      if (.not.allocated(deltad)) then
        allocate (deltad(ngroup*(ngroup-1)/2),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate deltad()'
      endif

      if (.not.allocated(rotinc)) then
        allocate (rotinc(nmo*(nmo-1)/2),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate rotinc()'
      endif

      if (.not.allocated(indi)) then
        allocate (indi(nmo),stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot allocate indi()'
      endif

      if (okw) write (lw,321)
 321  format (1x,'# ',60('+'),/,' # RUEDMIS routine',/,1x,'# ',60('+'))

      pi4  = atan(one)
      four = two+two
      pi   = four*pi4
      warn =.false.
c
c.....Print input group overlap integrals.
c
      if (okw.and.largwr) then
        write (lw,*) '# INPUT ATOMIC OVERLAP MATRIX BETWEEN NATURAL MOs'
      endif
      over(1:n,1:n)=zero
      do k=1,ngroup
        over(1:n,1:n)=over(1:n,1:n)+sg(k,1:n,1:n)
        if (okw.and.largwr) then
          write (lw,'(a,I2)') ' # FRAGMENT ',k
          do ip=1,n
             write (lw,30)(sg(k,ip,iq),iq=1,ip)
          enddo
        endif
      enddo
      if (okw) then
        if (largwr) then
          write (lw,*) '# SUM OVER FRAGMENTS'
          do ip=1,n
             write (lw,30)(over(ip,iq),iq=1,ip)
          enddo
        endif
        write (lw,233)
        write (lw,30) (natp(ip),ip=1,n)
      endif
c
c.....Initial rotation matrix is the unit matrix.
c
      c(1:n,1:n)=zero
      do i = 1, n
        c(i,i)=one
      enddo
c
c.....Start process.
c
      moreiter = .true.
      maxitn=maxit*n*(n-1)/2
      nrot=1
      goon=.true.
      do while (nrot.le.maxitn.and.goon)
        ipq=0
c
c.......Compute and print (DE)Localization funtion before rotating.
c
        xloctot=zero
        do ip=1,n
          xloci(ip)=zero
          do k = 1, ngroup
            spp = sg(k,ip,ip)
            xloci(ip) = xloci(ip) + spp * spp
          enddo
          xloctot=xloctot+xloci(ip)*natp(ip)*natp(ip)
        enddo
        if ((largwr.or.nrot.eq.1).and.okw) then
          write (lw,51) nrot,xloctot
          write (lw,54) 
          write (lw,30) (xloci(ip),ip=1,n)
          write (lw,520)
          write (lw,30)  (one/xloci(ip),ip=1,n)
        endif
c
c.......Find the pair of MOs increasing L the most.
c
        ipq=0
        xincmax = - 1q20
        do ip=1,n-1
          do iq=ip+1,n
            ipq=ipq+1
            apq = zero
            bpq = zero
            do k = 1, ngroup
              spp  = sg(k,ip,ip) * natp(ip)
              sqq  = sg(k,iq,iq) * natp(iq)
              spq  = sg(k,ip,iq) * sqrt(natp(ip)*natp(iq))
              spp2 = spp * spp
              sqq2 = sqq * sqq
              spq2 = spq * spq
              apq  = apq + spq2 - (spp2+sqq2-two*spp*sqq)/four
              bpq  = bpq + spq * (spp-sqq)
            enddo
            raiz = sqrt (apq*apq + bpq*bpq)
            xinc = apq + raiz
            rotinc(ipq)=xinc
            if (xinc.gt.xincmax) then
              xincmax=xinc
              ipair(1)=ip
              ipair(2)=iq
              apqmax=apq
              bpqmax=bpq
              raizmax=raiz
            endif
          enddo
        enddo
        if (okw.and.largwr) then
          write (lw,612)'Increase in L for the ',
     &                  'rotation of all pairs of MOs'
*       endif
*       if ((largwr.or.nrot.eq.1).and.okw) then
          imin=1
          do k=2,n
            imax=k*(k-1)/2
            write (lw,611) (rotinc(i),i=imin,imax)
            imin=imax+1
          enddo
          write (lw,'(a)') ' #'
        endif
c
        if (abs(xincmax).lt.epsq12) moreiter=.false.
        if (moreiter) then
c
c.........Sin(4*ALPHA) and Cos(4*ALPHA)
c
          sina=+bpqmax/raizmax
          cosa=-apqmax/raizmax
c
c.........Obtain the ALPHA angle in the range 0 le ALPHA le Pi/2
c
          if (bpqmax.ge.zero) then
            if (apqmax.lt.zero) then
              gam=asin(sina)/four
            else
              gam=(pi-asin(sina))/four
            endif
          else
            if (apqmax.lt.zero) then
              gam=(two*pi+asin(sina))/four
            else
              gam=(pi-asin(sina))/four
            endif
          endif
c
c.........Print rotation angle.
c
          ip=ipair(1)
          iq=ipair(2)
          gamdeg=gam*180d0/pi
          if (largwr.and.okw) write (lw,56) gamdeg,ip,iq
          cc=cos(gam)
          ss=sin(gam)
          cc2=cc*cc
          ss2=ss*ss
          cs=ss*cc
c
c.........Update populations
c
          xnuip=cc2*natp(ip)+ss2*natp(iq)
          xnuiq=ss2*natp(ip)+cc2*natp(iq)
c
c.........Compute the 2 x 2 non-unitary transformation matrix
c
          cpp=+cc*sqrt(natp(ip)/xnuip)
          cpq=+ss*sqrt(natp(iq)/xnuip)
          cqp=-ss*sqrt(natp(ip)/xnuiq)
          cqq=+cc*sqrt(natp(iq)/xnuiq)
c
c.........Modify cpp, cpq, cqp, cqq such that the transformed
c         orbitals are normalized.
c
          spqr3=zero
          do k=1,ngroup
            spqr3=spqr3+sg(k,ip,iq)
          enddo
          xnormp=cpp**2+cpq**2+two*cpp*cpq*spqr3
          xnormq=cqp**2+cqq**2+two*cqq*cqp*spqr3
          cpp=cpp/sqrt(xnormp)
          cpq=cpq/sqrt(xnormp)
          cqp=cqp/sqrt(xnormq)
          cqq=cqq/sqrt(xnormq)
c
          natp(ip)=xnuip
          natp(iq)=xnuiq
          if (largwr.and.okw) then
            write (lw,32)
            write (lw,31) cpp,cpq,cqp,cqq
          endif
c
c.........Update rotation matrix
c  
          do j=1,n
            ctmp(ip,j)=+cpp*c(ip,j)+cpq*c(iq,j)
            ctmp(iq,j)=+cqp*c(ip,j)+cqq*c(iq,j)
            c(ip,j)=ctmp(ip,j)
            c(iq,j)=ctmp(iq,j)
          enddo
c
c.........Write orbital electronic populations
c
          if (largwr.and.okw) then
            write (lw,*) '# New orbital populations'
            write (lw,30) (natp(j),j=1,n)
            write (lw,71) nrot,xloctot+xincmax
          endif
c
c.........Update the group overlap matrix.
c
          m=0
          do i=1,n
            if (i.ne.ip.and.i.ne.iq) then
              m=m+1
              indi(m)=i
            endif
          enddo
          do k=1,ngroup
            sgtmp(k,1:n,1:n)=sg(k,1:n,1:n)
            do ii=1,n-2
              i=indi(ii)
              tmp1=sg(k,i,ip)
              tmp2=sg(k,i,iq)
              sgtmp(k,i,ip)=+cc*tmp1+ss*tmp2
              sgtmp(k,i,iq)=-ss*tmp1+cc*tmp2
            enddo
            do jj=1,n-2
              j=indi(jj)
              tmp1=sg(k,ip,j)
              tmp2=sg(k,iq,j)
              sgtmp(k,ip,j)=+cc*tmp1+ss*tmp2
              sgtmp(k,iq,j)=-ss*tmp1+cc*tmp2
            enddo
            opp=sg(k,ip,ip)
            oqq=sg(k,iq,iq)
            opq=sg(k,ip,iq)
            oqp=sg(k,iq,ip)
            tmp1=+cs*opq+cs*oqp
            tmp2=-cs*opp+cs*oqq
            sgtmp(k,ip,ip)=cc2*opp+ss2*oqq+tmp1
            sgtmp(k,iq,iq)=ss2*opp+cc2*oqq-tmp1
            sgtmp(k,ip,iq)=tmp2+cc2*opq-ss2*oqp
            sgtmp(k,iq,ip)=tmp2+cc2*oqp-ss2*opq
            sg(k,1:n,1:n)=sgtmp(k,1:n,1:n)
          enddo
        endif
c
c.......Compute the property that should be zero at convergence
c
        if (nrot.eq.maxitn) warn = .true.
        if ((.not.moreiter).or.warn) then
c
c.........Compute and print (DE)Localization funtion after convergence
c
CC        xloctot=zero
CC        do ip=1,n
CC          xloci(ip)=zero
CC          do k = 1, ngroup
CC            spp = sg(k,ip,ip)
CC            xloci(ip) = xloci(ip) + spp * spp
CC          enddo
CC          xloctot=xloctot+xloci(ip)*natp(ip)*natp(ip)
CC        enddo
CC        if (okw) then
CC          write (lw,53) nrot,xloctot
CC          write (lw,54) 
CC          write (lw,30) (xloci(ip),ip=1,n)
CC          write (lw,520)
CC          write (lw,30)  (one/xloci(ip),ip=1,n)
CC          if (largwr) then
CC            write (lw,'(a)') ' #'
CC            write (lw,200)
CC            do ip=1,n
CC              write (lw,33) ip
CC              write (lw,30) (c(ip,iq),iq=1,n)
CC            enddo
CC          endif
CC          write (lw,23)
CC          write (lw,30) (natp(ip),ip=1,n)
CC        endif
c
c.........Test convergence. Print final results if converged.
c
          if (warn) then
            call error 
     &      ('_ruedmis_','Too many iterations (Change maxit)',warning)
            write (lw,1120) 
          else
            if (okw) then
              write (lw,'(a)') ' #'
              write (lw,300)
            else
              write (lw,400)
            endif
          endif
          if (okw) write (lw,'(a)') ' #'
c
c---------If a diagonal element is > 1.0, this element is set to 1.0,
c         the diagonal elements of other groups are set to 0.0, and
c         the non-diagonal elements of all the groups are set to 0.0.
c
          do i=1,nmo
            do k=1,ngroup
              if (sg(k,i,i).gt.one) then
                sg(1:ngroup,1:nmo,i)=zero
                sg(1:ngroup,i,1:nmo)=zero
                sg(k,i,i)=one
              endif
            enddo
          enddo
c
c.........Compute and print (DE)Localization funtion after convergence
c
          xloctot=zero
          do ip=1,n
            xloci(ip)=zero
            do k = 1, ngroup
              spp = sg(k,ip,ip)
              xloci(ip) = xloci(ip) + spp * spp
            enddo
            xloctot=xloctot+xloci(ip)*natp(ip)*natp(ip)
          enddo
          if (okw) then
            write (lw,510) xloctot
            if (largwr) then
              write (lw,610)
              imin=1
              do k=2,n
                imax=k*(k-1)/2
                write (lw,611) (rotinc(i),i=imin,imax)
                imin=imax+1
              enddo
            endif
            write (lw,'(a)') ' #'
            write (lw,54) 
            write (lw,30)  (xloci(ip),ip=1,n)
            write (lw,520) 
            write (lw,30)  (one/xloci(ip),ip=1,n)
            write (lw,'(a)') ' #'
            if (largwr) then
              write (lw,20)
              do ip=1,n
                write (lw,33) ip
                write (lw,30) (c(ip,iq),iq=1,n)
              enddo
              write (lw,23)
              write (lw,30) (natp(ip),ip=1,n)
              write (lw,66) 
              do k=1,ngroup
                write (lw,67) k
                write (lw,30) (sg(k,ip,ip)*natp(ip),ip=1,n)
              enddo
            endif
            if (largwr) then
              write (lw,'(a)') ' #'
              write (lw,'(a)') ' # FINAL ATOMIC OVERLAP MATRIX'
            endif
          endif
          over(1:n,1:n)=zero
          do k=1,ngroup
            over(1:n,1:n)=over(1:n,1:n)+sg(k,1:n,1:n)
            if (okw) then
              if (largwr) then
                write (lw,*) '# FRAGMENT ',k
                trace=zero
                do ip=1,n
                   trace=trace+sg(k,ip,ip)
                   write (lw,30)(sg(k,ip,iq),iq=1,ip)
                enddo
                write (lw,301) trace
              endif
            endif
          enddo
          if (okw) then
            if (largwr) then
              write (lw,*) '# SUM OVER FRAGMENTS'
              trace=zero
              do ip=1,n
                 trace=trace+over(ip,ip)
                 write (lw,30)(over(ip,iq),iq=1,ip)
              enddo
              write (lw,301) trace
            endif
            write (lw,'(a)') ' #'
            write (lw,'(a)') ' # Cioslowski-Mixon bond orders'
          endif
c
c.........Compute and write Cioslowski-Mixon bond orders and ionicity
c         indices
c
          npair=0
          do k=2,ngroup
            npair0=npair
            do l=1,k-1
              npair=npair+1
              del=zero
              delx=zero
              do ip=1,n
                del=del+natp(ip)**2*sg(k,ip,ip)*sg(l,ip,ip)
                do iq=1,n
                  delx=delx+natp(ip)*natp(iq)*sg(k,ip,iq)*sg(l,ip,iq)
                enddo
              enddo
              deltad(npair)=del
              delta (npair)=delx
            enddo
            if (okw) write (lw,38) k,k,k-1,(deltad(i),i=npair0+1,npair)
          enddo
          if (okw) then
            write (lw,'(a)') ' #'
            write (lw,'(a)') ' # Corrected Cioslowski-Mixon bond orders'
          endif
          npair=0
          do k=2,ngroup
            npair0=npair
            npair=npair+k-1
            if (okw) write (lw,38) k,k,k-1,(delta(i),i=npair0+1,npair)
          enddo
          if (okw) write (lw,1011) xmoloc*100
          do k=2,ngroup
            do l=1,k-1
              nb=0
              do ip=1,n
                deno=abs(sg(k,ip,ip))+abs(sg(l,ip,ip))
                if (abs(deno).gt.xmoloc) then
                  nb=nb+1
                  degloc(nb)=min(100d0*deno,100D0)
                  ixion(nb)=ip
                  xion(nb)=abs(sg(k,ip,ip))-abs(sg(l,ip,ip))
                  xion(nb)=100d0*abs(xion(nb)/deno)
                endif
              enddo
              if (okw) write (lw,69) k,l,
     &           (ixion(ip),xion(ip),degloc(ip),ip=1,nb)
            enddo
          enddo
          goon=.false.
        endif
        nrot=nrot+1
      enddo
c
      if (allocated (over)) then
        deallocate (over,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate over()'
      endif
      if (allocated (sgtmp)) then
        deallocate (sgtmp,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate sgtmp()'
      endif
      if (allocated (xloci)) then
        deallocate (xloci,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate xloci()'
      endif
      if (allocated (xion)) then
        deallocate (xion,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate xion()'
      endif
      if (allocated (ixion)) then
        deallocate (ixion,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate ixion()'
      endif
      if (allocated (degloc)) then
        deallocate (degloc,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate degloc()'
      endif
      if (allocated (ctmp)) then
        deallocate (ctmp,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate ctmp()'
      endif
      if (allocated (delta)) then
        deallocate (delta,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate delta()'
      endif
      if (allocated (deltad)) then
        deallocate (deltad,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate deltad ()'
      endif
      if (allocated (rotinc)) then
        deallocate (rotinc,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate rotinc()'
      endif
      if (allocated (indi)) then
        deallocate (indi,stat=ier)
        if (ier.ne.0) stop 'ruedmis.f: Cannot deallocate indi()'
      endif
c
      if (okw) then
        write (lw,*) '#'
        write (lw,*) '# Done _ruedmis'
      endif
      call timer (4,ipid,'_ruedmis  ',-1)
      return
c
c.....Formats
c
 20   format (' # Rotation Matrix (READ BY ROWS)',
     &        '    locMO_i =  SUM_j C_ij natMO_j')
 200  format (' # ROTATION MATRIX (READ BY ROWS)',
     &        '    locMO_i =  SUM_j C_ij natMO_j')
 30   format (10(1x,F12.6))
 301  format (5x,'Trace =  ',F12.6)
 31   format (2(/,2(5x,E14.8)))
 32   format (1x,'# Transformation matrix')
 300  format (' #',/,' #',80('-'),/,' # CONVERGED LOCALIZATION',/,' #')
 400  format (1x,'# CONVERGED LOCALIZATION')
 51   format (' #',80('-'),/,' #',/,
     &        ' # L value before rotation ' ,I4,' = ',E16.10)
 71   format (' # L value after  rotation ' ,I4,' = ',E16.10)
 53   format (' # L value after completing ',I4,' rotations = ',E16.10)
 510  format (' # L value             = ',E16.10)
 54   format (' # L_i (i=1..N)')
 520  format (' # Fragments expanded by MOs, L_i^(-1) (i=1..N)')
 56   format (1x,'# ROTATION ANGLE = ',F12.6,
     &  ' degree for the pair of atoms ',2I3,/,' #')
 23   format (/1x,'# Rotated orbitals natural populations')
 233  format (/1x,'# INPUT orbitals natural populations')
 1120 format (1x,'#',/,1x,'# ',80('+'),/,
     &  ' # UNCONVERGED RESULTS',/,' #')
 66   format (' # Fragment occupation numbers of the localized MOs')
 67   format (' # FRAGMENT ',I2)
 38   format (' # delta(',I4,' 1) to delta(',2I4,') = ',/,
     &        10(1x,F8.5))
*69   format (' # Ionicity indices for the pair of atoms ',2I4,/,
*    &   ' # The associated MO is indicated previous to the index',/,
*    &   6(2x,I5,1x,F5.1,'%'))
 69   format (' # Pair of atoms ',2I4,/,
     &        5(1x,I3,1x,F5.1,'%',' (',F6.2'%)'))
 610  format (' #',/,
     &   ' # Final potential increase in L for the rotation of all',
     &   ' pairs of MOs')
 612  format (' # ',a,a)
 611  format (10(1x,E12.6))
 33   format (' # Localized MO number ',I3)
 1011 format (' #',/,
     &' # Cioslowski-Mixon ionicity indices: ',/,
     &' # ---------------------------------- ',/,
     &' # Only MOs localized more than ',F5.2,'% between the two',1x,
     &'atoms of the pair are considered',/,
     &' # The associated MO is indicated previous to the index',/,
     &' # and the degree of localization between the two atoms in',1x,
     &'parenthesis')
      end
