      recursive subroutine searchbin (arr,key,imi,ima,Nint,ifound,found)
      implicit none
      integer(kind=8), intent (inout) :: arr(Nint)
      integer(kind=8), intent (inout) :: key
      integer        , intent (in   ) :: Nint
      integer                         :: imi,ima,imid,ifound
      logical                         :: found
!
      if (ima < imi) then
        found = .false.
        return
      else
        imid = (imi+ima)/2
             if ( arr(imid) > key ) then
          call searchbin (arr,key,imi,imid-1,Nint,ifound,found)
        else if ( arr(imid) < key ) then
          call searchbin (arr,key,imid+1,ima,Nint,ifound,found)
        else
          found = .true.
          ifound = imid
          return
        endif
      endif
      end
