c 
c.......................................................................
c
      subroutine coremos 
     &   (sg,ng,nmo,ovc,mocore,core,icore,moval,ival,lw,ok,largwr,line)
      include     'implicit.inc'
      real(kind=8) sg(ng,nmo,nmo)
      integer(kind=4) icore(nmo),ival(nmo)
      integer(kind=4) core(ng)
      integer(kind=4) icogrp(ng,nmo)
      logical iscore,ok,largwr
      parameter (lline  = 2000)  
      integer(kind=4) totcore
      character(len=4) fourchar
      character*(lline) line
      integer(kind=4) leng
c
c-----------------------------------------------------------------------
c
      write (lw,65) 
      mocore = 0
      core   = 0
      do ig=1,ng
        do i=1,nmo
          if (largwr) then
            write (lw,66) ig,i,sg(ig,i,i)
          endif
          if (sg(ig,i,i).gt.ovc) then
            mocore=mocore+1
            core(ig)=core(ig)+1
            icore(mocore)=i
            icogrp(ig,core(ig))=i
          endif
        enddo
      enddo
c
c.....Re-compute the valence orbitals
c
      moaux = nmo-mocore
      moval = 0
      do i=1,nmo
        j=1
        iscore=.false.
        do while (j.le.mocore.and.(.not.iscore))
          if (i.eq.icore(j)) iscore=.true.
          j=j+1
        enddo
        if (.not.iscore) then
          moval=moval+1
          ival(moval)=i
        endif
      enddo
      if (moval.ne.moaux) then
        write (0,*)
        write (0,*) '# coremos.f: Wrong COREMO order !!!'
        write (0,*)
        stop
      endif
c
      write (lw,223) ovc
      totcore=0
      do ig=1,ng
        write (lw,224) ig,core(ig),(icogrp(ig,k),k=1,core(ig))
        totcore=totcore+core(ig)
        do k=1,core(ig)
           ick=icogrp(ig,k)
           write (lw,226) ick,ick,ig,sg(ig,ick,ick)
        enddo
      enddo
      ok=.false.
      if (totcore.le.lline/4-1) then
        ok=.true.
        line(1:4) = '    '
        iw=5
        do ig=1,ng
          if (core(ig).gt.0) then
            do k=1,core(ig)
             line(iw:iw+4) = fourchar(icogrp(ig,k))
             do j=iw,iw+4
               if (line(j:j).eq.'0') then
                 line(j:j)=' '
               else
                 exit
               endif
             enddo
             iw=iw+4
            enddo
          endif
        enddo
*       write (6,'(a)') line(1:leng(line))
      endif
c
c.....Write list of CORE and VALENCE orbitals.
c
      if (mocore.gt.0) then
        write (lw,29) mocore,(icore(i),i=1,mocore)
      else
        write (lw,29) mocore
      endif
      if (moval.gt.0) then
         write (lw,31) moval,(ival(i),i=1,moval)
      else
         write (lw,31) moval
      endif
      return
      stop
 223  format (1x,'#',/,
     & 1x,'#',/,1x,'# Automatic analysis of CORE MOs',/,
     & 1x,'# Overlap criterium of Localization = ',F16.10)
 224  format (1x,'# FRAGMENT ',I3,' HAS ',I3,
     & ' MOs ALMOST FULLY LOCALIZED ON IT: ',/,100(1x,'# ---> ',20I3,/))
 226  format (1x,'# <MO_',I3,'|MO_',I3,'>_',I2,' = ',F16.10)
 29   format (//,' # THERE ARE ',I3,' CORE    ORBITALS :',200I3)
 31   format (   ' # THERE ARE ',I3,' VALENCE ORBITALS :',200I3)
 66   format (' # Group ',I2,'  MO ',I4,'  Self-Overlap = ',F15.8)
 65   format (//,'  Automatic analysis of CORE MOs')
      end
