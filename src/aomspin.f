c
c-----------------------------------------------------------------------
c
      module space_for_aomspin
      save
      real(kind=8), allocatable,dimension (:,:,:) :: aomalpha,aombeta
      end module space_for_aomspin
c
c-----------------------------------------------------------------------
c
      subroutine aomspin (ncent,nmo,aom,nalpha,nbeta)
      USE space_for_wfnbasis
      USE space_for_aomspin
      real   (kind=8) aom(ncent,nmo,nmo)
      integer(kind=4) ncent,nmo,nalpha,nbeta,i
c
      allocate (aomalpha(ncent,nalpha,nalpha))
      allocate (aombeta (ncent,nbeta ,nbeta ))
      do i=1,ncent
        aomalpha(i,:,:) = aom(i,ialpha(1:nalpha),ialpha(1:nalpha))
        aombeta (i,:,:) = aom(i,ibeta (1:nbeta ),ibeta (1:nbeta ))
      enddo

*     do i=1,ncent
*       do j=1,nalpha
*         do k=1,nalpha
*           aomalpha(i,j,k)=aom(i,ialpha(j),ialpha(k))
*         enddo
*       enddo
*       do j=1,nbeta
*         do k=1,nbeta
*           aombeta(i,j,k)=aom(i,ibeta(j),ibeta(k))
*         enddo
*       enddo
*     enddo

      end subroutine aomspin
