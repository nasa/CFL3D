      subroutine bc9999(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,iuns,nou,bou,nbuf,
     .                  ibufdim,nummem,x,y,z,smin)
c
c***********************************************************************
c     Purpose:  Set exact solution BCs
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension smin(jdim-1,kdim-1,idim-1)
c
      common /sklton/ isklton
      common /maxiv/ ivmx
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
c
      if (iexact_trunc .ne. 0 .or. iexact_disc .ne. 0) then
        call exact_flow_bc(jdim,kdim,idim,x,y,z,qi0,qj0,qk0,nface,
     +                       bcj,bck,bci,smin,tj0,tk0,ti0,
     +                       vj0,vk0,vi0,nummem,
     +                       iexact_trunc,iexact_disc)
      else
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'('' Error... bc9999 only usable with'',
     +   '' iexact_trunc .ne. 0 or iexact_disc .ne. 0'')')
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      return
      end
