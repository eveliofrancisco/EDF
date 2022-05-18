c
c-----------------------------------------------------------------------
c
      subroutine fno (aominp,inside,rhf,rohf,uhf,cciqa,icorr,mal,mbe,
     & natoms,nprims,nmo,ncent,epsneg,udat,wfnfile,lw,lerr,largwr)
c
c-----PARAMETERS--------------------------------------------------------
c
c     aominp()                        INPUT/OUTPUT
c
c     Atomic Overlap Matrix (AOM) of Canonical MOs in all the centers.
c
c-----AOMO()                          OUTPUT
c
c     AOM between orthonormalized Localized MOs in all the centers in the
c     output
c
c.....c()                            OUTPUT 
c
c     Rotation matrix giving the nonorthonormal localized MOs in terms of 
c     the input MOs: Loc_i = SUM_j C_ij CANMO_j
c
c.....mxcent                         INPUT
c
c     Maximum number of N in the search of (Nc,2e) bonds.
c
c.....covx                           INPUT
c
c     Factor multiplying the sum of covalent radii of two atoms A and B  
c     used in the connect routine to determined whether A and B are 
c     linked or not from a purely geomtrical point of view.
c
c-----epsneg                         INPUT
c
c     Diagonalizations in this routine involve definite positive 
c     matrices that should have all of their eigvenvalues positive.
c     Numerical errors, however, make that this not necessarily happens 
c     and some of the eigenvalues can be marginally negative. The value 
c     of 'epsneg' is the largest value (in absolute value) allowed in order
c     that the method ends nicely. For instance, if epsneg=-1D-6
c     and all the eigenvalues are greater that this number the 
c     process continues up to the end. However, if an eigvenvalue
c     is smaller that -1D-6, for instance -1D-5, this routine stops.
c
c-----maxb                          INPUT
c
c     Maximum number of links between two atoms A and B for which
c     a (2c,2e) localized MOs between these atoms is searched. Atoms
c     directly connected in a geometrical sence have a 'geometrical
c     distance' or link of 1, if A is connected to C which, in turn,
c     is connected to B the 'geometrical distance' or link between
c     A and B is 2, and so on. 
c
c.....nloc                           OUTPUT
c
c     Number of localized MOs found. It should be equal to nloctarget
c     (see below) if the routine ends nicely
c
c.....nloctarget                      INPUT
c
c     Half the number of electrons of the molecule (assumed to have
c     a even number of electrons, i.e. nalpha = nbeta)
c
c.....udat                           INPUT
c
c     Unit number of original WFN file
c
c-----wfnfile                        INPUT
c
c     Name of the WFN file
c
c-----maxcritov                      INPUT
c
c     Maximum order 'n' in the search of (nc,2e) localized MOs
c
c.....critov(1)...critov(mxcent)     INPUT
c
c     Tolerance used in localization: If an eigvenvalue of
c     S^(1/2) rho S^{1/2} is greater that critov(i) it means that the
c     associated eigenvector corresponds to an MO assumed to be loca-
c     lized in the 'i' centers involved in the diagonalization.
c
c.....ixcent(1)...ixcent(mxcent)     INPUT
c
c     Logical variable: If ixcent(i)=.true. (i-center,2e) bonds are
c     searched, otherwise (i-center,2e) bonds are skipped.
c
c-----largwr                         INPUT
c
c     Logical variable: if largwr=.true. a large output is requested,
c     otherwise a short output is used.
c
c-----warn                           OUTPUT
c
c     Returned  .false. if the routine ends nicely. Otherwise warn is
c     returned .true.
c
c-----skiph                          INPUT
c
c     Logical variable: If .true. hydrogen atoms are skipped in one-
c     center localizations. Otherwise hydrogens are not skipped.
c
c-----lw                             INPUT
c
c     Logical output unit
c
c-----lerr                           INPUT
c
c     Logical output unit of error messages
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
      real    (kind=8) epsneg
      logical          largwr
      integer(kind=4)  udat,udatnw
      character(len=*) wfnfile
      character (len=mline) wfnloc

      real   (kind=8),allocatable, dimension (:,:,:) :: aom,aomn
      real   (kind=8),allocatable, dimension (:,:,:) :: aomx
      real   (kind=8),allocatable, dimension (:,:,:) :: aomab
      real   (kind=8),allocatable, dimension (:,:)   :: sg
      real   (kind=8),allocatable, dimension (:,:)   :: uvec
      real   (kind=8),allocatable, dimension (:,:)   :: cmat
      real   (kind=8),allocatable, dimension (:,:)   :: rhosig
      real   (kind=8),allocatable, dimension (:,:)   :: v
      real   (kind=8),allocatable, dimension (:,:)   :: vx
      real   (kind=8),allocatable, dimension (:,:)   :: dcoef
      real   (kind=8),allocatable, dimension (:,:)   :: cc
      real   (kind=8),allocatable, dimension (:)     :: eigen
      real   (kind=8),allocatable, dimension (:)     :: eigx
      real   (kind=8),allocatable, dimension (:)     :: eigs
      real   (kind=8),allocatable, dimension (:)     :: xloci
      real   (kind=8),allocatable, dimension (:)     :: uveceig
      real   (kind=8),allocatable, dimension (:)     :: ocup
      real   (kind=8),allocatable, dimension (:)     :: work
      integer(kind=4),allocatable, dimension (:)     :: iord
      integer(kind=4) isigma(nmo)
      character (len=1) spin

c
c-----------------------------------------------------------------------
c
      call timer (2,ifno,'_fno      ',-1)
!
!-----CCIQA wavefunctions with Nalpha # Nbeta can not be used yet.
!
      if (cciqa.and.mal.ne.mbe) then
        write (lerr,*) '# fno.f: Open-shell CCWFN not allowed yet'
        return
      endif 

      ncode = 1  ! Nalpha = Nbeta and Ialpha():Ibeta()
      if (uhf.or.rohf.or.(icorr.and.mal.ne.mbe)) ncode=2
c
c-----Store temporarily the Occupation number of WFN file in occwfn()
c
      occwfn(1:nmo) = occ(1:nmo)
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
        n3 = 3*nab-1
        allocate (aom(1:ncent,1:nab,1:nab))
        allocate (aomn(1:ncent,1:nab,1:nab))
        allocate (aomx(1:ncent,1:nab,1:nab))
        allocate (aomab(1:ncent,1:nab,1:nab))
        allocate (sg(1:nab,1:nab))
        allocate (uvec(1:nab,1:nab))
        allocate (cmat(1:nab,1:nab))
        allocate (rhosig(1:nab,1:nab))
        allocate (v(1:nab,1:nab))
        allocate (vx(1:nab,1:nab))
        allocate (dcoef(1:nab,1:nprims))
        allocate (cc(1:nab,1:nprims))
        allocate (eigen(1:nab))
        allocate (eigx(1:nab))
        allocate (eigs(1:nab))
        allocate (xloci(1:nab))
        allocate (uveceig(1:nab))
        allocate (iord(1:nab))
        allocate (work(1:n3))
CCC     aomab(:,1:nab,1:nab)=aominp(:,isigma(1:nab),isigma(1:nab))
CCC     rhosig(1:nab,1:nab)=rho(isigma(1:nab),isigma(1:nab))
c
c-------aomab = Block of the original AOM matrix that corresponds
c       exclusively to sigma electrons (sigma=alpha or beta)
c
        sg = 0d0
        do i=1,nab
          do j=1,nab
            rhosig(i,j)=rho(isigma(i),isigma(j))
            do k=1,ncent
              aomab(k,i,j)=aominp(k,isigma(i),isigma(j))     
            enddo
            do k=1,natoms
              sg(i,j) = sg(i,j)+aomab(inside(k),i,j)
            enddo
          enddo
        enddo
c
c-------If the fragment is denoted as A, sg() = S^A. Diagonalize S^A
c
        call jacobi (sg,nab,nab,eigen,uvec,nrot)
         write (6,*) eigen
          STOP
c
c-------Test that sg() is definite positive
c
        do m=1,nab
          if (eigen(m).gt.0d0) then
          elseif (eigen(m).gt.epsneg) then
            write (lerr,222) eigen(m),inside(1:natoms)
            eigen(m)=0d0
          else
            write (lerr,1124) m,eigen(m)
            return
          endif
        enddo
c
c-------Compute S^{A,1/2} storing it in the sg() array.
c
        sg = 0d0
        do k=1,nab
          do l=1,k
            forall (m=1:nab) uveceig(m)=uvec(l,m)*sqrt(eigen(m))
            sg(k,l) = dot_product(uvec(k,:),uveceig(:))
            sg(l,k) = sg(k,l)
          enddo
        enddo
c
c-------Now, sg() is the array [S^A]^(1/2)
c
        vx = matmul(sg,matmul(rhosig,sg))
c
c-------Diagonalize rho^A
c
        call jacobi (vx,nab,nab,eigs,uvec,nrot)
        write (lw,987) (inside(j),j=1,natoms)
        write (lw,988) (eigs(j),j=1,nab)
        do j=1,nab
          write (lw,990) j,(uvec(i,j),i=1,nab)
        enddo
c
c-------Compute AOM between the FNOs, that are orthonormal in A, in 
c       all the atoms of the system
c
        do i=1,ncent
          aom(i,:,:) = matmul(transpose(uvec),matmul(aomab(i,:,:),uvec))
        enddo
c
c-------At this point, the new MOs, say |psi>, are not orthonormal in R^3,
c       but <psi|psi>_{R^3} = S # I.  However, they are orthonormal in A
c       <psi|psi>_{A} = I. To obtain FNOs, |phi>, that are orthonormal in 
c       R^3, we will diagonalize S and apply |phi> = |psi> S^{-1/2}
c       
c-------Compute S and store it in the array sg()
c
        sg = 0d0
        do i=1,ncent
          sg = sg + matmul(transpose(uvec),matmul(aomab(i,:,:),uvec))
        enddo
c
c-------Diagonalize sg()
c
        call jacobi (sg,nab,nab,eigen,v,nrot)
c
c-------Compute S^{-1/2} storing it in the sg() array.
c
        sg = 0d0
        do k=1,nab
          do l=1,k
            if (eigen(m).lt.0d0) then
              write (lerr,1125) m,eigen(m)
              stop
            endif
            if (abs(eigen(m)).lt.1d-6) then
              write (lerr,1126) m,eigen(m)
              stop
            endif
            forall (m=1:nab) uveceig(m)=v(l,m)/sqrt(eigen(m))
            sg(k,l) = dot_product(v(k,:),uveceig(:))
            sg(l,k) = sg(k,l)
          enddo
        enddo
c
c-------Finally, Compute AOM between the latter FNOs (that now are 
c       orthonormal in R^3) in all the atoms of the system.
c
        do i=1,ncent
          aomn(i,:,:) = matmul(sg,matmul(aom(i,:,:),sg))
        enddo
c
c-------Store aomn() in aom()
c
        aom = aomn
c
c-------Write diagonal AOM elements in each atom or fragment. 
c
        write (lw,68) 'FNO'
        write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
        do ip=1,nab
          xloci(ip) = 0d0  ! Effective number of centers of each FNO
          do k = 1, ncent
            spp = aom(k,ip,ip)
            xloci(ip) = xloci(ip) + spp * spp
          enddo
          write (lw,410,advance='no') ip
          if (ncent.gt.10) then
            write (lw,41 ) (100d0*aom(k,ip,ip),k=1,10)
            write (lw,411) (100d0*aom(k,ip,ip),k=11,ncent)
            write (lw,*)
          else
            write (lw,41) (100d0*aom(k,ip,ip),k=1,ncent)
          endif
        enddo
        write (lw,520) (1d0/xloci(ip),ip=1,nab)
c
c-------Find FNOs in the primitive basis
c
        forall (i=1:nab) dcoef(i,:) = coef(nmo+isigma(i),:)
*       cc = matmul(v,dcoef)

        cc = matmul(transpose(v),dcoef)
c
c-------write a WFN file with FNOs
c
        occ(1:nab) = eigs(1:nab)
        if (ncode.eq.1) occ(1:nab) = 2d0 * eigs(1:nab)
        inpr = nprims
        nloc = nab
        wfnloc = trim(wfnfile)//"-fno"
        if (ncode.eq.2) wfnloc = trim(wfnfile)//"-fno"//spin(1:1)
        do i=1,natoms
          wfnloc = trim(wfnloc)//"-"//fourchar(inside(i))
        enddo
        udatnw = udat+10
        open (unit=udatnw,file=trim(wfnloc),status='unknown')
        write (lw,*)
        write (lw,351) trim(wfnloc)
        allocate (ocup(nab),stat=ier)
        if (ier.ne.0) stop '# newfno.f: Cannot allocate ocup()'
        ocup(1:nab)=occ(1:nab)
        call cdafhmos (udat,udatnw,0,ocup,cc,nab,inpr)
        deallocate (ocup,stat=ier)
        if (ier.ne.0) stop '# newfno.f: Cannot deallocate ocup()'
c
c-------Write file with AOMs between FNOs
c
        lu18 = udat+11
        wfnloc = trim(wfnfile)//"-fnoaom"
        if (ncode.eq.2) wfnloc = trim(wfnfile)//"-fno"//spin(1:1)
        do i=1,natoms
          wfnloc = trim(wfnloc)//"-"//fourchar(inside(i))
        enddo
        open (unit=lu18,file=trim(wfnloc),status='unknown')
        write (lw,353) trim(wfnloc)
        do k=1,ncent
          write (lu18,'(I4,a)') k,' <=== AOM within this center'
          write (lu18,80) ((aom(k,m,j),m=1,j),j=1,nab)
        enddo
        close (unit=lu18)
c
c-------Analyze localization in all the centers of FNOs
c
        critics = 0.02d0
        aomx = aom
*       call aomtest (aomx,nab,ncent,critics,lw)
        write (lw,3340) 'FNOs'
        call newdrvmult (nab,nprims,ncent,cc,lw)

        deallocate (aom)
        deallocate (aomn)
        deallocate (aomx)
        deallocate (aomab)
        deallocate (sg)
        deallocate (uvec)
        deallocate (cmat)
        deallocate (rhosig)
        deallocate (v)
        deallocate (vx)
        deallocate (dcoef)
        deallocate (cc)
        deallocate (eigen)
        deallocate (eigx)
        deallocate (eigs)
        deallocate (xloci)
        deallocate (uveceig)
        deallocate (iord)
        deallocate (work)
      enddo
c
c-----Recover the original occupation numbers
c
      occ(1:nmo) = occwfn(1:nmo)
      call timer (4,ifno,'_fno      ',-1)
      return
c
c.....Formats
c
 410  format (' # ',I4)
 41   format (1000(10x,10(1x,F5.1,'%')))
 411  format (1000(17x,10(1x,F5.1,'%'),/))
 68   format (/' # Localization of each ',a,' OQS MO in each atom')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 711  format (' # SUM  ',17x,'= ',F18.8)
 351  format (1x,"# Writing file '",a,"'",1x,"with FNOs")
 353  format (1x,"# Writing file '",a,"'",1x,"with AOM for FNOs")
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 1125 format (' # Eigenvalue ',I4,' of S_[R^3] is negative, ',E17.10)
 1126 format (' # Eigenvalue ',I4,' of S_[R^3] is too low, ',E17.10)
 28   format (' # The value of EPSNEG is ',E15.8)
 222  format (' # Negative eigenvalue ',E15.8,
     &  ' set to 0.0 for fragment ',100I4)
 101  format (1000(' # ',5(1x,F15.8),/))
 80   format (6(1x,e16.10))
 3340 format (/,' # Second & fourth moment orbital spread of ',a)
 520  format (' # Atoms expanded by each FNO in each FNO',/,
     &    1000(5(1x,F15.8),/))
 987  format (' # Fragment formed by atoms:',/,1000(20(1x,I3),/))
 988  format (' # Eigenvalues ',/,1000(5(1x,F15.8),/))
 989  format (' # Eigenvectors (canonical basis)')
 990  format (' # Eigenvector ',I3,/,1000(5(1x,F15.8),/))
c
      end
