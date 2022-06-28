c
c-----------------------------------------------------------------------
c
      subroutine otherfno (aominp,inside,rhf,rohf,uhf,cciqa,icorr,
     & mal,mbe,natoms,nprims,nmo,ncent,epseigen,udat,wfnfile,lw,lerr,
     & verbose)
c
c-----PARAMETERS--------------------------------------------------------
c
c-----aominp()                        INPUT
c
c     Atomic Overlap Matrix (AOM) of Canonical MOs in all the centers.
c
c-----INSIDE (1..NATOMS)              INPUT 
c
c     Indices of the NATOMS that define the fragment.
c
c-----rhf,rohf,uhf,cciqa,icorr        INPUT
c
c     Each of them is .TRUE. or .FALSE. depending on the type of WFN used
c
c-----mal,mbe                         INPUT
c
c     Number of ALPHA and BETA electrons     
c
c-----natoms                          INPUT
c
c     Number of atoms of the fragment (fragment is called A from now on)
c
c-----nprims                          INPUT
c
c     Number of primitives of the WFN
c
c-----nmo                             INPUT 
c
c     Number of canonical MOs in the WFN
c
c-----ncent                           INPUT
c
c     Number of atoms of the molecule
c
c
c-----epseigen                        INPUT
c
c     In the diagonalization of [S^A]^(1/2) * 1RDM * [S^A]^(1/2), the
c     eigenvalues smaller than 'epseigen', as well as the associated
c     eigenvectors, correspond to FNOs mostly localized outside the 
c     fragment A and will ge ignored after the above diagonalizaition 
c     is performed.
c
c.....udat                            INPUT
c
c     Unit number of original WFN file
c
c-----wfnfile                         INPUT
c
c     Name of the WFN file
c
c-----lw                              INPUT
c
c     Output unit
c
c-----lerr
c
c     Output unit for errors
c
c-----verbose                         INPUT
c
c     Logical variable: if verbose=.true. a large output is requested,
c     otherwise a short output is used.
c
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfncoef
      USE      space_for_wfnbasis
      USE      space_for_rdm1
      include 'implicit.inc'
      include 'corr.inc'
      include 'mline.inc'
      real   (kind=8) aominp(ncent,nmo,nmo)
      integer(kind=4) inside(ncent)
      integer(kind=4) insidex(ncent)
      real   (kind=8) rho(nmo,nmo)
      real   (kind=8) occwfn(nmo)
      integer(kind=4) natoms
      integer(kind=4) mal,mbe
      integer(kind=4) iordatoms(ncent)
      integer(kind=4) typewfn
      logical rhf,rohf,uhf,cciqa,icorr
      character(len=4) fourchar
      logical          verbose
      integer(kind=4)  udat,udatnw
      character(len=*) wfnfile
      character (len=mline) wfnloc

      real   (kind=8),allocatable, dimension (:,:,:) :: aom
      real   (kind=8),allocatable, dimension (:,:,:) :: aomab
      real   (kind=8),allocatable, dimension (:,:)   :: sg
      real   (kind=8),allocatable, dimension (:,:)   :: sgx
      real   (kind=8),allocatable, dimension (:,:)   :: sginv
      real   (kind=8),allocatable, dimension (:,:)   :: cmat
      real   (kind=8),allocatable, dimension (:,:)   :: cmatord
      real   (kind=8),allocatable, dimension (:,:)   :: wvec
      real   (kind=8),allocatable, dimension (:,:)   :: uvec
      real   (kind=8),allocatable, dimension (:,:)   :: ueu
      real   (kind=8),allocatable, dimension (:,:)   :: rhosig
      real   (kind=8),allocatable, dimension (:,:)   :: srhos
      real   (kind=8),allocatable, dimension (:,:)   :: ufin
      real   (kind=8),allocatable, dimension (:,:)   :: rhoa
      real   (kind=8),allocatable, dimension (:,:)   :: dcoef
      real   (kind=8),allocatable, dimension (:,:)   :: cc
      real   (kind=8),allocatable, dimension (:)     :: epsi
      real   (kind=8),allocatable, dimension (:)     :: lambda
      real   (kind=8),allocatable, dimension (:)     :: eigx
      real   (kind=8),allocatable, dimension (:)     :: xloci
      real   (kind=8),allocatable, dimension (:)     :: diagord
      real   (kind=8),allocatable, dimension (:)     :: diagaom
      real   (kind=8),allocatable, dimension (:)     :: ocup
      integer(kind=4),allocatable, dimension (:)     :: iaom
      integer(kind=4),allocatable, dimension (:)     :: iloc
      integer(kind=4) isigma(nmo)
      character (len=1) spin

c
c-----------------------------------------------------------------------
c
      call timer (2,iotherfno,'_otherfno ',-1)
!
!-----CCIQA wavefunctions with Nalpha # Nbeta can not be used yet.
!
      if (cciqa.and.mal.ne.mbe) then
        write (lerr,*) '# otherfno.f: Open-shell CCWFN not allowed yet'
        return
      endif 

      ncode = 1  ! Nalpha = Nbeta and Ialpha():Ibeta()
      if (uhf.or.rohf.or.(icorr.and.mal.ne.mbe)) ncode=2
      xfac = 1d0
      if (ncode.eq.1) xfac = 2d0
c
c-----Ordering the atoms of the fragment.
c
      insidex = inside
      forall (i=1:natoms) iordatoms(i)=i
      call iqcksort (insidex, iordatoms, ncent, 1, natoms)
      forall (i=1:natoms) inside(i)=insidex(iordatoms(i))
      do icode=1,ncode
        if (icode.eq.1) then
          spin = 'a'
          nab = nalpha
          if (rhf.or.(icorr.and.mal.eq.mbe)) then
            rho = 0.5d0 * c1et
          else
            rho = c1ea
          endif
          isigma(1:nab)=ialpha(1:nab)
        else
          spin = 'b'
          nab = nbeta
          rho = c1eb
          isigma(1:nab)=ibeta(1:nab)
        endif
        allocate (aomab(1:ncent,1:nab,1:nab))
        allocate (sg(1:nab,1:nab))
        allocate (sgx(1:nab,1:nab))
        allocate (wvec(1:nab,1:nab))
        allocate (rhosig(1:nab,1:nab))
        allocate (epsi(1:nab))
        allocate (iloc(1:nab))
        allocate (srhos(nab,nab))
c
c-------Extract from the current AOM matrix the part associated
c       to electrons with spin ALPHA or BETA
c
        aomab(:,1:nab,1:nab)=aominp(:,isigma(1:nab),isigma(1:nab))
c
c-------Extract from the current 1RDM matrix the part associated
c       to electrons with spin ALPHA or BETA
c
        rhosig(1:nab,1:nab)=rho(isigma(1:nab),isigma(1:nab))
c
c-------Group overlap matrix (S^A) in the fragment
c
        sg = 0d0
        do i=1,natoms
          sg = sg+aomab(inside(i),:,:)
        enddo
        sgx=sg
c
c-------Diagonalize S^A
c
        call jacobi (sgx,nab,nab,epsi,wvec,nrot)
c
c-------Determine the number of non-zero eigenvalues (nloc) and 
c       their indices.
c
        nloc=0
        do i=1,nab
          if (abs(epsi(i)).gt.epseigen) then
            nloc=nloc+1
            iloc(nloc)=i
          endif
        enddo
c
        allocate (eigx(nloc))
        allocate (ufin(nab,nloc))
        allocate (rhoa(nloc,nloc))
        allocate (lambda(1:nloc))
        allocate (uvec(1:nloc,1:nloc))
        allocate (cmat(1:nab,1:nloc))
        allocate (ueu(nloc,nloc))
        allocate (aom(1:ncent,1:nloc,1:nloc))
        allocate (xloci(1:nloc))
        allocate (dcoef(1:nab,1:nprims))
        allocate (cc(1:nloc,1:nprims))
c
c-------Store in eigx() the non-zero eigenvalues and associated 
c       eigenvectors. Actually, we store in ufin(:,i) the ith-eigenvector 
c       divided by SQRT(eigx(i)).
c
        do i=1,nloc
          eigx(i)=epsi(iloc(i))
        enddo
c
c-------Write the relevant eigenvalues and eigenvectors (=FNOs)
c
        
        write (lw,987) (inside(i),i=1,natoms)
        if (verbose) then
          write (lw,*) '# Fragment Overlap matrix of canonical MOs'
          do k=1,nab
            write (lw,'(5(1x,F15.8))') (sg(k,l),l=1,k)
          enddo
        endif
        if (ncode.eq.1) then
          write (lw,*) '# ALPHA and BETA FNOs are equivalent'
          write (lw,*) '# ALPHA and BETA FNOs are determined'
        else
          if (icode.eq.1) write (lw,*) '# ALPHA FNOs are determined'
          if (icode.eq.2) write (lw,*) '# BETA  FNOs are determined'
        endif
        write (lw,43) epseigen
        write (lw,42) epseigen,nloc
        write (lw,988) (eigx(i),i=1,nloc)
        do i=1,nloc
          if (verbose) write (lw,990) i,(wvec(j,iloc(i)),j=1,nab)
          ufin(:,i)=wvec(:,iloc(i))/sqrt(eigx(i))
        enddo
c
c-------Obtain the matrix representation of rho in the basis of the
c       above eigenvectors
c
        srhos=matmul(sg,matmul(rhosig,sg))
        rhoa=matmul(transpose(ufin),matmul(srhos,ufin))
        call jacobi (rhoa,nloc,nloc,lambda,uvec,nrot)
c
c-------AOM in R^3 of FNOs
c
        do i=1,nloc
          do j=1,nloc
            value=0d0
            do k=1,nloc
              value=value+uvec(k,i)*uvec(k,j)/eigx(k)
            enddo
            ueu(i,j)=value
          enddo
        enddo

        if (verbose) then
          write (lw,*) ' # AOM in R^3 of FNOs'
          do k=1,nloc
            write (lw,'(5(1x,F15.8))') (ueu(k,l),l=1,k)
          enddo
        endif
c
c-------Matrix transformation from Canonical MOs to FNOs
c
        cmat=matmul(ufin,uvec)
        if (verbose) then
          write (lw,*) ' # Matrix from CANMOs to FNOs'
          do k=1,nloc
            write (lw,*) ' # FNO ',k
            write (lw,'(5(1x,F15.8))') (cmat(l,k),l=1,nab)
          enddo
        endif
c
c-------Normalized FNOs in R^3
c-------cmat() contains now the Normalized FNOs in R^3
c
        forall (i=1:nloc) cmat(:,i)=cmat(:,i)/sqrt(ueu(i,i))
c-------Overlaps in R^3 of normalized in R^3 FNOs
c
        do i=1,nloc
          do j=1,nloc
            ueu(i,j)=ueu(i,j)/sqrt(ueu(i,i)*ueu(j,j))
          enddo
        enddo
        if (verbose) then
          write (lw,*) ' # AOM in R^3 of the normalized in R^3 FNOs'//
     &                 ', .i.e. AOM in R^3 of NFNOs'
          do k=1,nloc
            write (lw,'(5(1x,F15.8))') (ueu(k,l),l=1,k)
          enddo
        endif
c
c-------Compute AOM between normalized in R^3 FNOs in all the atoms
c
        do i=1,ncent
          aom(i,:,:) = matmul(transpose(cmat),matmul(aomab(i,:,:),cmat))
        enddo
c
c-------Ordering FNOs by decreasing value of the diagonal overlap in the
c       group. 
c
        allocate (diagaom(1:nloc))
        allocate (diagord(1:nloc))
        allocate (iaom(1:nloc))
        diagaom=0d0
        do ip=1,nloc
          iaom(ip)=ip
          do i=1,natoms
            diagaom(ip) = diagaom(ip) + aom(inside(i),ip,ip)
          enddo
        enddo
        call qcksort (diagaom, iaom, 1, nloc)
        forall (i=1:nloc) diagord(i)=diagaom(iaom(nloc-i+1))
        write (lw,332) (diagord(i),i=1,nloc)

        write (lw,68) 'normalized in R^3 FNO'
        write (lw,*) '#    MO\ATOM --->   1       2       3       4  ..'
        do ip=1,nloc
          i=iaom(nloc-ip+1)
          xloci(ip) = 0d0  ! Effective number of centers of each FNO
          do k = 1, ncent
            spp = aom(k,i,i)
            xloci(ip) = xloci(ip) + spp * spp
          enddo
          write (lw,410,advance='no') ip
          if (ncent.gt.10) then
            write (lw,41 ) (100d0*aom(k,i,i),k=1,10)
            write (lw,411) (100d0*aom(k,i,i),k=11,ncent)
            write (lw,*)
          else
            write (lw,41) (100d0*aom(k,i,i),k=1,ncent)
          endif
        enddo
        write (lw,991) (xfac*lambda(iaom(nloc-ip+1)),ip=1,nloc)
        allocate (ocup(nloc),stat=ier)
        if (ier.ne.0) stop '# otherfno.f: Cannot allocate ocup()'
        forall (ip=1:nloc) ocup(ip)=xfac*lambda(iaom(nloc-ip+1))
        write (lw,994) xfac*sum(lambda(1:nloc))
        write (lw,520) (1d0/xloci(ip),ip=1,nloc)
c
c-------Find FNOs in the primitive basis, cc(). 
c
        allocate (cmatord(1:nab,1:nloc))
        do i=1,nloc
          cmatord(:,i) = cmat(:,iaom(nloc-i+1))
        enddo
        if (verbose) then
          write (lw,*) '# Normalized in R^3 FNOs from canonical MOs'
          do i=1,nloc
            write (lw,*) ' Normalized in R^3 FNO number ',i
            write (lw,'(5(1x,F15.8))') (cmatord(j,i),j=1,nab)
          enddo 
        endif
        do i=1,nab
          dcoef(i,:) = coef(nmo+isigma(i),:)
        enddo
        cc = matmul(transpose(cmatord),dcoef)
c
c-------write a WFN file with FNOs
c
        inpr = nprims
        wfnloc = trim(wfnfile)//"-fno"
        if (ncode.eq.2) wfnloc = trim(wfnfile)//"-fno"//spin(1:1)
        do i=1,natoms
          wfnloc = trim(wfnloc)//"-"//fourchar(inside(i))
        enddo
        udatnw = udat+10
        open (unit=udatnw,file=trim(wfnloc),status='unknown')
        write (lw,*)
        write (lw,351) trim(wfnloc)
        call cdafhmos (udat,udatnw,0,ocup,cc,nloc,inpr)
        deallocate (ocup,stat=ier)
        if (ier.ne.0) stop '# otherfno.f: Cannot deallocate ocup()'
c
        deallocate (aom)
        deallocate (aomab)
        deallocate (sg)
        deallocate (sgx)
        deallocate (cmat)
        deallocate (cmatord)
        deallocate (wvec)
        deallocate (uvec)
        deallocate (rhosig)
        deallocate (srhos)
        deallocate (rhoa)
        deallocate (ufin)
        deallocate (ueu)
        deallocate (dcoef)
        deallocate (cc)
        deallocate (epsi)
        deallocate (lambda)
        deallocate (eigx)
        deallocate (xloci)
        deallocate (iloc)
        deallocate (diagaom)
        deallocate (diagord)
        deallocate (iaom)
      enddo
      call timer (4,iotherfno,'_otherfno   ',-1)
      return
c
c.....Formats
c
 410  format (' # ',I4)
 41   format (1000(10x,10(1x,F6.2,'%')))
 411  format (1000(17x,10(1x,F6.2,'%'),/))
 68   format (/' # Localization of each ',a,' in each atom')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 711  format (' # SUM  ',17x,'= ',F18.8)
 351  format (1x,"# Writing file '",a,"'",1x,"with NFNOs",/,' # ',
     &        '(Normalized in R^3 Orbitals are written)')
 353  format (1x,"# Writing file '",a,"'",1x,"with AOM for FNOs")
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 28   format (' # The value of EPSNEG is ',E15.8)
 222  format (' # Negative eigenvalue ',E15.8,
     &  ' set to 0.0 for fragment ',100I4)
 101  format (1000(' # ',5(1x,F15.8),/))
 80   format (6(1x,e16.10))
 3340 format (/,' # Second & fourth moment orbital spread of ',a)
 520  format (' # Atoms expanded by each FNO',/,1000(5(1x,F15.8),/))
 987  format (//' # ',77('+'),/,
     & ' # FNO ANALYSIS - FNO ANALYSIS - FNO ANALYSIS - FNO ANALYSIS',/
     & ' # FNOs are not normalized in R^3',/,
     & ' # NFNOs ARE normalized in R^3',/,
     & ' # ',77('+'),/,' # Fragment formed by atoms:',/,
     & 1000(' # ',20(1x,I3),/))
 988  format (' # Eigenvalues of fragment overlap matrix',
     & /,1000(5(1x,F15.8),/))
 994    format (' # Full set of FNOs TOTAL POPULATION = ',F15.8,/,
     & ' # ',77('-'))
 991  format (' # ',77('-'),/' # FNO occupation numbers',/,
     &  1000(5(1x,F15.8),/))
 989  format (' # Eigenvectors (canonical basis)')
 990  format (' # Eigenvector ',I3,/,1000(5(1x,F15.8),/))
 42   format (' # Eigenvalues    > ',1PE15.8,' are ',I4)
 43   format (' # EPSEIGEN value = ',1PE15.8)
 332  format (' # Diagonal Overlaps in the fragment'
     & ' of normalized in R^3 FNOs',/,1000(5(1x,F15.8),/))
c
      end
