      subroutine parser(inpstr,npos,lc1,lc2,lcl,iflg)
c***********************************************************************
c     Purpose:  parse 210-character data
c
c     inpstr = input character string
c     npos   = sub-string position to find after column lc2-on-input
c     lc1    = first column of the sub-string (on output)
c     lc2    = last  column of the sub-string (on output)
c            = last  column of the previous sub-string (on input)
c     lcl    = last column of the last sub-string
c     iflg   = 1: do normal parsing: start from front and extract
c             -1: only find lcl
c
c     original routine provided by James O. Hager, Boeing Phantom Works,
c     Long Beach, CA
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      parameter (linp=210)
      character*210 inpstr
c
c     normal parsing operation

      if (iflg.eq.1) then
         lc1 = max(lc2+1,1)
c
c        loop over the sub-strings
c
         do ipos=1,npos
            if(ipos.gt.1) lc1 = lc2 +1
c
c           find the first non-blank position
c
            do lc=lc1,linp
               if(inpstr(lc:lc).ne.' ') go to 1100
            end do
 1100       continue
            lc1 = lc
c
c           find the last non-blank position
c
            do lc=lc1,linp
               if(inpstr(lc:lc).eq.' ') go to 1600
            end do
 1600       continue
            lc2 = lc -1
         end do
c
c     find lc2 of the last sub-string (lcl)
c
      else
         do lc=linp,1,-1
            if(inpstr(lc:lc).ne.' ') go to 2100
         end do
 2100    continue
         lcl = lc
      end if
c
      return
      end
