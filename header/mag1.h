      parameter(mdim=14,jdim=300,kdim=300,ldim=2,mhldim=10,
     .          ibdim=2000,idim=8000,ndim=300,ipmax=28,
     .          iwrdim=idim+ibdim,mlen=jdim*kdim*ldim)
c
c     parameter definitions:
c
c     mdim................max number of input grids
c
c     jdim,kdim,ldim......max j, k, l dimensions of any grid 
c
c     mhldim..............max number of holes/outer boundaries (per grid)   
c
c     idim................max number of fringe points (total over all grids)
c                         (also max number of boundary points)
c 
c     ibdim...............max number of points allowed for storage of 
c                         the initial hole hole boundary (per hole)
c    
c     ndim................largest dimension of any coordinate surface used
c                         for hole cutting
c
c     ipmax...............max number of coordinate surfaces used for cutting
c                         holes (total over all grids) 
c
