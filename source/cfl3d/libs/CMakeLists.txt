cmake_minimum_required(VERSION 2.6)

set(
CM_FSRCS
add2x.F        resadd.F       blocki_d.F     blockj_d.F 
init_mast.F    blockk_d.F     mreal.F        barth3d.F 
outbuf.F       get_bvals.F    qface.F        xmukin.F 
rotateq2_d.F   abciz.F        diagk.F        wkstn.F 
abcjz.F        colldat.F      dlutr.F        abckz.F 
collv.F        dlutrp.F       init.F         q8vrev.F 
collx.F        dsmin.F        initvist.F     lamfix.F 
collxtb.F      expand.F       int2.F         bc9999.F 
extra.F        extrae.F       lead.F 
amafi.F        ffluxl.F       l2norm.F       bc_xmera.F 
amafj.F        dabciz.F       fhat.F         l2norm2.F 
amafk.F        dabcjz.F       fill.F         blkmax.F 
avgint.F       dabckz.F       loadgr.F       newfit.F 
blocki.F       dfbtr.F        getdhdr.F      ld_qc.F 
blockj.F       dfbtrp.F       gfluxl.F       project.F 
blockk.F       dfhat.F        grdmove.F      q8sdot.F 
bsub.F         dfluxpm.F      hfluxl.F       q8smax.F 
bsubp.F        diagi.F        q8smin.F       spalart.F 
ccf.F          diagj.F        hole.F         global0.F 
tau2x.F        tdq.F          tinvr.F        tmetric.F 
trace.F        trans.F        transmc.F      transp.F 
vlutrp.F       xe.F           xlim.F         unld_qc.F 
xtbatb.F       rie1d.F        rie1de.F       rotateq.F 
rotateq0.F     rotatmc.F      rotatp.F       twoeqn.F 
rp3d.F         setblk.F       setdqc0.F      triv.F 
setqc0.F       swafi.F        swafj.F        swafk.F 
af3f.F         resid.F        update.F       cellvol.F 
chkrap.F       chksym.F       gfluxv.F       hfluxv.F 
ffluxv.F       bcchk.F        blomax.F       dird.F 
direct.F       ctime1.F       dthole.F       diagnos.F 
fa2xi.F        fa2xj.F        fa2xk.F        ffluxr.F 
gfluxr.F       hfluxr.F       fluxm.F        fluxp.F 
force.F        coll2q.F       collq.F        collqc0.F 
collapse.F     collxt.F       rotate.F       rsmooth.F 
shear.F        topol2.F       xe2.F          vlutr.F 
rechk.F        metric.F       invert.F       ld_datj.F 
ld_datk.F      ld_dati.F      topol.F        arc.F 
wmag.F         chkrot.F       chkrotj_d.F    chkrotk_d.F 
chkroti_d.F    int2_d.F       i2x.F          i2xj_d.F 
i2xk_d.F       i2xi_d.F       i2xs.F         i2xsj_d.F 
i2xsk_d.F      i2xsi_d.F      readdat.F      cblki.F 
cblkj.F        cblkk.F        cblki_d.F      cblkj_d.F 
cblkk_d.F      bc2005.F       bc2005j_d.F    bc2005k_d.F 
bc2005i_d.F    pre_blockbc.F  pre_blocki.F   pre_blockj.F 
pre_blockk.F   bc_info.F      bc1000.F       bc1001.F 
bc1002.F       bc1003.F       bc1005.F       bc1008.F 
bc1011.F       bc1012.F       bc1013.F       bc2002.F 
bc2003.F       bc2004.F       bc2006.F       bc2007.F 
bc2008.F       bc2102.F       rotateqb.F     rpatch.F 
delv.F         xupdt.F        intrbc.F       fa.F 
tau.F          histout.F      csurf.F        csout.F 
mvdat.F        getibk.F       pre_embed.F    pre_period.F 
pre_patch.F    pre_cblki.F    pre_cblkj.F    pre_cblkk.F 
global.F       global2.F      bc.F           getibk0.F 
rpatch0.F      cntsurf.F      cctogp.F       bc2009.F 
modread.F      genforce.F     rotsurf.F      mms.F
trnsurf.F      tfiedge.F      tfiface.F      tfivol.F  
arclen.F       deform.F       bc_delt.F      getsurf.F 
setcorner.F    prolim.F       prolim2.F      readkey.F 
parser.F       ae_pred.F      collmod.F      init_ae.F 
aesurf.F       chkdef.F       moddefl.F      xyzintr.F 
delintr.F      rcfl.F         ccomplex.F     setseg.F 
cgnstools.F    fmaps.F        getdelt.F      rsurf.F 
avghole.F      augmntq.F      blnkfr.F       gradinfo.F 
rb_pred.F      rb_corr.F      init_rb.F      init_trim.F 
pltmode.F      bcnonin.F      initnonin.F    resnonin.F 
my_flush.F     bc2016.F       gfluxv1.F      hfluxv1.F 
ffluxv1.F      sijrate2d.F    threeeqn.F     lesdiag.F 
sijrate3d.F    foureqn.F      bc2026.F       bc2019.F 
)
set(
CM_ONLY_FSRCS
u_doubleprime.F bc2010.F
)
set(
CM_F90SRCS
module_profileout.F90 module_contour.F90  module_kwstm.F90 module_stm_2005.F90
)
set(
CM_CSRCS
bessel.c
)
set(
CM_FSRCS_SPEC
addx.F
)

set(CM_SRCS ${CM_F90SRCS} ${CM_FSRCS} ${CM_FSRCS_SPEC} ${CM_ONLY_FSRCS})
add_library(COBJ OBJECT ${CM_CSRCS})

add_library(common $<TARGET_OBJECTS:COBJ> ${CM_SRCS})
target_compile_options(common PUBLIC ${CFL3D_FLAGS})
if(USE_CGNS)
	target_link_libraries(common PUBLIC cgns)
endif(USE_CGNS)

if(BUILD_CMPLX)
	# complex version of libcommon do not contain F90 source files:  
	# module_profileout.F90 module_contour.F90  module_kwstm.F90 module_stm_2005.F90
	set(
		CM_C_ONLY_FSRCS
		histout_img.F
	)
	set(CM_C_SRCS ${CM_FSRCS} ${CM_FSRCS_SPEC} ${CM_ONLY_FSRCS} ${CM_C_ONLY_FSRCS})
	add_library(common_cmplx $<TARGET_OBJECTS:COBJ> ${CM_C_SRCS})
	
	target_compile_options(common_cmplx PUBLIC ${CFL3D_FLAGS})
	target_compile_options(common_cmplx PUBLIC -DCMPLX)
	if(USE_CGNS)
		target_link_libraries(common_cmplx PUBLIC cgns)
	endif(USE_CGNS)
endif(BUILD_CMPLX)
