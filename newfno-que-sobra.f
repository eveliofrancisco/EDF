
        rhoa=matmul(transpose(wvec),matmul(rhosig,wvec))
        do i=1,nab
          do j=1,nab
            rhoa(i,j)=rhoa(i,j) * sqrt(epsilon(i) * epsilon(j))
          enddo
        enddo
c
c-------Diagonalize 
c
        call jacobi (rhoa,nab,nab,lambda,uvec,nrot)
        allocate (cmat(1:nab,1:nab))
        cmat=matmul(wvec,uvec)
        nloc=0
        do i=1,nab
*         if (abs(lambda(i)).gt.epseigen) then
            nloc=nloc+1
            iloc(nloc)=i
*         endif
        enddo
        allocate (eigx(nloc))
        allocate (ufin(nab,nloc))
c
c-------Store in eigx() the NLOC eigenvalues eigs > EPSEIGEN
c
        forall (j=1:nloc) eigx(j)=lambda(iloc(j))
c
c-------Store in ufin() the associated NLOC eigenvectors
c
        forall (j=1:nloc) ufin(:,j)=cmat(:,iloc(j))
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
*       allocate (vv(1:nloc,1:nloc))
*       vv=matmul(transpose(ufin),matmul(sginv,ufin))
c
c-------Write the diagonal elements and their inverses
c
*       write (lw,*) ' # Diagonal <psi|psi>_{R^3} elements and inverses'
*       do i=1,nloc
*         write (lw,'("  #",1x,2E15.8)') vv(i,i),1d0/vv(i,i)
*       enddo
c
c-------Compute the matrix C = S^{-1/2} * U
c
*       allocate (cmat(1:nab,1:nloc))
*       cmat=matmul(sginv12,ufin)
*       write (lw,*) ' # CMAT = Matrix that give the DNOS from CANMOs'
*       do l=1,nloc
*         write (lw,*) ' # FNO ',l
*         write (lw,'(5(1x,F15.8))') (cmat(k,l),k=1,nab)
*         write (lw,*) 
*       enddo
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
        write (lw,*) 'AOM in A'
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
        write (lw,*) 'AOM in R^3 (It should be the unit matrix)'
        do k=1,nloc
          write (lw,'(5(1x,F15.8))') (aomina(k,l),l=1,k)
        enddo
c
c-------Normalize in R^3 the FNOs
c
CCCC    forall (i=1:nloc) cmat(:,i)=cmat(:,i)/sqrt(aomina(i,i))
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
        call cdafhmos (udat,udatnw,0,cc,nloc,inpr)
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
