c
c-----------------------------------------------------------------------
c
      subroutine fragno (aominp,inside,rhf,rohf,uhf,cciqa,icorr,mal,mbe,
     & natoms,nprims,nmo,ncent,epsneg,epseigen,smalleigen,udat,wfnfile,
     & lw,lerr,largwr)
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
c-----epsneg                          INPUT
c
c     The AOM matrix in fragment A, S^A, is positive definite and should 
c     have all of its eigvenvalues positive. Numerical errors, however, 
c     make that this not necessarily happens and some of the eigenvalues 
c     of S^A can be marginally negative. The value of 'epsneg' is the 
c     largest value (in absolute value) allowed in order that the method 
c     ends nicely. For instance, if epsneg=-1D-6 and all the eigenvalues 
c     of S^A are greater that this number the process continues up to the 
c     end. However, if an eigvenvalue is smaller that -1D-6, for instance 
c     -1D-5, this routine stops. When an eigenvalue of S^A is negative
c     but greater than epsneg, such eigenvalue is set to 0.0
c
c-----epseigen                        INPUT
c
c     In the diagonalization of [S^A]^(1/2) * 1RDM * [S^A]^(1/2), the
c     eigenvalues smaller than 'epseigen', as well as the associated
c     eigenvectors, correspond to FNOs mostly localized outside the 
c     fragment A and will ge ignored after the above diagonalizaition 
c     is performed.
c
c-----smalleigen
c
c-----In the computation of [S^A]^(-1/2) a division by the square root
c     of an eigenvalue of S^A is involved. In many cases, such eigenvalue
c     is zero or very small. To avoid a division by 0.0 the eigenvalue
c     of S^A is replaced by 'smalleigen' if it is smaller than this number.
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
c-----largwr                         INPUT
c
c     Logical variable: if largwr=.true. a large output is requested,
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
      real    (kind=8) epsneg
      logical          largwr
      integer(kind=4)  udat,udatnw
      character(len=*) wfnfile
      character (len=mline) wfnloc

      real   (kind=8),allocatable, dimension (:,:,:) :: aom
      real   (kind=8),allocatable, dimension (:,:,:) :: aomab
      real   (kind=8),allocatable, dimension (:,:)   :: sg
      real   (kind=8),allocatable, dimension (:,:)   :: sginv12
      real   (kind=8),allocatable, dimension (:,:)   :: sginv
      real   (kind=8),allocatable, dimension (:,:)   :: cmat
      real   (kind=8),allocatable, dimension (:,:)   :: aomina
      real   (kind=8),allocatable, dimension (:)     :: w2
      real   (kind=8),allocatable, dimension (:,:)   :: uvec
      real   (kind=8),allocatable, dimension (:,:)   :: rhosig
      real   (kind=8),allocatable, dimension (:,:)   :: vv
      real   (kind=8),allocatable, dimension (:,:)   :: ufin
      real   (kind=8),allocatable, dimension (:,:)   :: rhoa
      real   (kind=8),allocatable, dimension (:,:)   :: dcoef
      real   (kind=8),allocatable, dimension (:,:)   :: cc
      real   (kind=8),allocatable, dimension (:)     :: eigen
      real   (kind=8),allocatable, dimension (:)     :: eigs
      real   (kind=8),allocatable, dimension (:)     :: eigx
      real   (kind=8),allocatable, dimension (:)     :: xloci
      real   (kind=8),allocatable, dimension (:)     :: uveceig
      real   (kind=8),allocatable, dimension (:)     :: uveceiginv
      real   (kind=8),allocatable, dimension (:)     :: ocup
      integer(kind=4),allocatable, dimension (:)     :: ipiv
      integer(kind=4),allocatable, dimension (:)     :: iloc
      integer(kind=4) isigma(nmo)
      character (len=1) spin

c
c-----------------------------------------------------------------------
c
      call timer (2,ifragno,'_fragno   ',-1)
!
!-----CCIQA wavefunctions with Nalpha # Nbeta can not be used yet.
!
      if (cciqa.and.mal.ne.mbe) then
        write (lerr,*) '# fragno.f: Open-shell CCWFN not allowed yet'
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
        allocate (aomab(1:ncent,1:nab,1:nab))
        allocate (sg(1:nab,1:nab))
        allocate (sginv12(1:nab,1:nab))
        allocate (uvec(1:nab,1:nab))
        allocate (rhosig(1:nab,1:nab))
        allocate (eigen(1:nab))
        allocate (eigs(1:nab))
        allocate (uveceig(1:nab))
        allocate (uveceiginv(1:nab))
        allocate (iloc(1:nab))
c
c-------Extract from the current AOM matrix the part associated
c       to electrons with spin sigma = alpha or beta
c
        aomab(:,1:nab,1:nab)=aominp(:,isigma(1:nab),isigma(1:nab))
c
c-------Extract from the current 1RDM matrix the part associated
c       to electrons with spin sigma = alpha or beta
c
        rhosig(1:nab,1:nab)=rho(isigma(1:nab),isigma(1:nab))
c
c-------Group overlap matrix (S^A) in the fragment
c
        sg = 0d0
        do i=1,natoms
          sg = sg+aomab(inside(i),:,:)
        enddo
c
c-------Inverse of S^A
c
        allocate (sginv(1:nab,1:nab))
        allocate (w2(1:nab))
        allocate (ipiv(1:nab))
        sginv=sg
        call dgetrf (nab,nab,sginv,nab,ipiv,info)
        if (info.ne.0) then
          stop ' # fragno.f: Error in the dgetrf.f subroutine'
        endif
        call dgetri (nab,sginv,nab,ipiv,w2,nab,info)
        if (info.ne.0) then
          stop ' # fragno.f: Error in the dgetri.f subroutine'
        endif
c
c-------Diagonalize S^A
c
        call jacobi (sg,nab,nab,eigen,uvec,nrot)
c
c-------Test that sg()=S^A is definite positive
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
c-------Compute [S^A]^{1/2} storing it in the sg() array.
c-------Compute also [S^A]^{-1/2} storing it in the sginv12() array.
c
*       SMALLEIGEN = 1D-10
        sg = 0d0
        sginv12 = 0d0
        do k=1,nab
          do l=1,nab
            forall (m=1:nab) uveceig(m)=uvec(l,m)*sqrt(eigen(m))
c
c-----------Here avoid division by zero. If an eigvenvalue of S^A,
c           say eigen(m), is smaller than SMALLEIGEN, eigen(m) is 
c           replaced by SMALLEIGEN in the following division.
c
            forall (m=1:nab) 
              uveceiginv(m)=uvec(l,m)/sqrt(max(eigen(m),SMALLEIGEN))
            endforall
            sg(k,l) = dot_product(uvec(k,:),uveceig(:))
            sginv12(k,l) = dot_product(uvec(k,:),uveceiginv(:))
          enddo
        enddo
c
c-------Compute rho^A = [S^A]^{1/2} * 1RDM * [S^A]^{1/2} 
c
        allocate (rhoa(1:nab,1:nab))
        rhoa = matmul(sg,matmul(rhosig,transpose(sg)))
        call jacobi (rhoa,nab,nab,eigs,uvec,nrot)
c
c-------Select eigevalues larger than EPSEIGEN
c
*       EPSEIGEN = 1.0d-04
        nloc=0
        do i=1,nab
          if (abs(eigs(i)).gt.EPSEIGEN) then
            nloc=nloc+1
            iloc(nloc)=i
          endif
        enddo
        allocate (eigx(nloc))
        allocate (ufin(nab,nloc))
c
c-------Store in eigx() the NLOC eigenvalues eigs > EPSEIGEN
c
        forall (j=1:nloc) eigx(j)=eigs(iloc(j))
c
c-------Store in ufin() the associated NLOC eigenvectors
c
        forall (j=1:nloc) ufin(:,j)=uvec(:,iloc(j))
c
c-------Write the relevant eigenvalues and eigenvectors (=FNOs)
c
        write (lw,987) (inside(j),j=1,natoms)
        write (lw,43) epsneg,epseigen,smalleigen
        write (lw,42) EPSEIGEN,nloc
        write (lw,988) (eigx(j),j=1,nloc)
        do j=1,nloc
          write (lw,990) j,(ufin(i,j),i=1,nab)
        enddo
c
c-------Compute the overlap matrix of relevant FNOs in R^3 and store
c       it in vv()
c
        allocate (vv(1:nloc,1:nloc))
        vv=matmul(transpose(ufin),matmul(sginv,ufin))
c
c-------Write the diagonal elements and their inverses
c
        write (lw,*) ' # Diagonal <psi|psi>_{R^3} elements and inverses'
        do i=1,nloc
          write (lw,'("  #",1x,2E15.8)') vv(i,i),1d0/vv(i,i)
        enddo
c
c-------Compute the matrix C = S^{-1/2} * U
c
        allocate (cmat(1:nab,1:nloc))
        cmat=matmul(sginv12,ufin)
        write (lw,*) ' # CMAT = Matrix that give the DNOS from CANMOs'
        do l=1,nloc
          write (lw,*) ' # FNO ',l
          write (lw,'(5(1x,F15.8))') (cmat(k,l),k=1,nab)
          write (lw,*) 
        enddo
c
c-------Compute AOM between FNOs in all the atoms
c
        allocate (aom(1:ncent,1:nloc,1:nloc))
        do i=1,ncent
          aom(i,:,:) = matmul(transpose(cmat),matmul(aomab(i,:,:),cmat))
        enddo
c
c-------AOM between FNOs in fragment A
c
        allocate (aomina(nloc,nloc))
        write (lw,*) 'AOM in A (It should be the unit matrix)'
        do k=1,nloc
          do l=1,nloc
            aomina(k,l)=sum(aom(inside(1:natoms),k,l))
          enddo
        enddo
        do k=1,nloc
          write (lw,'(5(1x,F15.8))') (aomina(k,l),l=1,k)
        enddo
c
c-------Compute AOM between FNOs in R^3
c
        aomina = matmul(transpose(cmat),cmat)
        write (lw,*) 'AOM in R^3'
        do k=1,nloc
          write (lw,'(5(1x,F15.8))') (aomina(k,l),l=1,k)
        enddo
c
c-------Normalize in R^3 the FNOs
c
        forall (i=1:nloc) cmat(:,i)=cmat(:,i)/sqrt(aomina(i,i))
c
c-------Test that the diagonal overlap in R^3 of the normalized FNOs
c       are computed Ok!
c
        do i=1,ncent
          aomina = matmul(transpose(cmat),cmat)
        enddo
        write (lw,*) ' # AOM in R^3 of FNOs normalized in R^3'
        do k=1,nloc
          write (lw,'(5(1x,F15.8))') (aomina(k,l),l=1,k)
        enddo
c
c-------Compute AOM of the FNOs normalized in R^3 in all the atoms
c
        do i=1,ncent
          aom(i,:,:) = matmul(transpose(cmat),matmul(aomab(i,:,:),cmat))
        enddo
c
c-------Write diagonal AOM elements in each atom or fragment. 
c
        allocate (xloci(1:nloc))
        write (lw,68) 'FNO'
        write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
        do ip=1,nloc
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
        write (lw,520) (1d0/xloci(ip),ip=1,nloc)
c
c-------Find normalized in R^3 FNOs in the primitive basis, cc().
c
        allocate (dcoef(1:nab,1:nprims))
        do i=1,nab
          dcoef(i,:) = coef(nmo+isigma(i),:)
        enddo
        allocate (cc(1:nloc,1:nprims))
        cc = matmul(transpose(cmat),dcoef)
c
c-------write a WFN file with FNOs
c
        do i=1,nloc
          occ(i)=eigx(i)
        enddo
        if (ncode.eq.1) occ(1:nloc) = 2d0 * occ(1:nloc)
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
        allocate (ocup(nloc),stat=ier)
        if (ier.ne.0) stop '# fragno.f: Cannot allocate ocup()'
        ocup(1:nloc)=occ(1:nloc)
        call cdafhmos (udat,udatnw,0,ocup,cc,nloc,inpr)
        deallocate (ocup,stat=ier)
        if (ier.ne.0) stop '# fragno.f: Cannot deallocate ocup()'
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
          write (lu18,80) ((aom(k,m,j),m=1,j),j=1,nloc)
        enddo
        close (unit=lu18)
c
c-------Analyze localization in all the centers of FNOs
c
CCCCC   critics = 0.02d0
CCCCC   allocate (aomx(1:ncent,1:nloc,1:nloc))
CCCCC   aomx = aom
CCCCC   call aomtest (aomx,nab,ncent,critics,lw)
        write (lw,3340) 'FNOs'
        call newdrvmult (nloc,nprims,ncent,cc,lw)

        deallocate (aom)
        deallocate (aomab)
        deallocate (aomina)
        deallocate (sg)
        deallocate (sginv)
        deallocate (sginv12)
        deallocate (cmat)
        deallocate (w2)
        deallocate (ipiv)
        deallocate (uvec)
        deallocate (rhosig)
        deallocate (vv)
        deallocate (rhoa)
        deallocate (ufin)
        deallocate (dcoef)
        deallocate (cc)
        deallocate (eigen)
        deallocate (eigs)
        deallocate (eigx)
        deallocate (xloci)
        deallocate (uveceig)
        deallocate (uveceiginv)
        deallocate (iloc)
      enddo
c
c-----Recover the original occupation numbers
c
      occ(1:nmo) = occwfn(1:nmo)
      call timer (4,ifragno,'_fragno   ',-1)
      return
c
c.....Formats
c
 410  format (' # ',I4)
 41   format (1000(10x,10(1x,F6.2,'%')))
 411  format (1000(17x,10(1x,F6.2,'%'),/))
 68   format (/' # Localization of each ',a,' OQS MO in each atom')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 711  format (' # SUM  ',17x,'= ',F18.8)
 351  format (1x,"# Writing file '",a,"'",1x,"with FNOs")
 353  format (1x,"# Writing file '",a,"'",1x,"with AOM for FNOs")
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 28   format (' # The value of EPSNEG is ',E15.8)
 222  format (' # Negative eigenvalue ',E15.8,
     &  ' set to 0.0 for fragment ',100I4)
 101  format (1000(' # ',5(1x,F15.8),/))
 80   format (6(1x,e16.10))
 3340 format (/,' # Second & fourth moment orbital spread of ',a)
 520  format (' # Atoms expanded by each FNO in each FNO',/,
     &    1000(5(1x,F15.8),/))
 987  format (' # ',77('+'),/,
     & ' # FNO ANALYSIS - FNO ANALYSIS - FNO ANALYSIS - FNO ANALYSIS',/
     & ,' # ',77('+'),/,' # Fragment formed by atoms:',/,
     & 1000(20(1x,I3),/))
 988  format (' # Eigenvalues ',/,1000(5(1x,F15.8),/))
 989  format (' # Eigenvectors (canonical basis)')
 990  format (' # Eigenvector ',I3,/,1000(5(1x,F15.8),/))
 42   format (' # Eigenvalues > ',F15.8,' are ',I4)
 43   format (' # EPSNEG     value = ',E15.8,/,
     &        ' # EPSEIGEN   value = ',E15.8,/,
     &        ' # SMALLEIGEN value = ',E15.8)
c
      end
