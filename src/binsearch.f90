      recursive subroutine binsearch (array,key,imin,imax,ifound,found)
      implicit none
      integer(kind=8) key
      integer imin,imax,imid,ifound
      logical found
      integer (kind=8) array(*)
!
!-----Test if arrau is empty
!
      if (imax < imin) then
!
!       set is empty, so return value showing not found
!
        found = .false.
        return
      else
!
!       calculate midpoint to cut set in half
!
        imid = (imin+imax)/2
!
!       three-way comparison
!
             if ( array(imid) > key ) then
!
!            key is in the lower subset
!
          call binsearch (array,key,imin,imid-1,ifound,found)
        else if ( array(imid) < key ) then
!
!            key is in the upper subset
!
          call binsearch (array,key,imid+1,imax,ifound,found)
        else
          found = .true.
          ifound = imid
          return
        endif
      endif
      end
