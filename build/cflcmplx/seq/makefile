#     $Id$
#=============================================================================
#
#        creates the executable for the sequential version of cfl3d
#
#=============================================================================

# ***************************** CREATE LINKS *********************************

link:     lncode lnhead lnlibs

lncode:
	@ echo "        linking source code"
	ln -s $(CFLSRC_D)/*.F .
	ln -s $(CFLSRC_S)/ccomplex.F .

lnhead:

lnlibs:
	ln -s ../libs/$(COMMONLIB) .

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# **************************** CREATE LIBRARIES  *****************************

#Note: fsrc_dist dependancy list must not contain development.f

FSRC_DIST = \
	bc_blkint.F    findmin_new.F  plot3d.F       rrest.F \
	bc_patch.F     forceout.F     plot3t.F       rrestg.F \
	calyplus.F     pointers.F     setup.F        writ_buf.F \
	mgblk.F        qinter.F       prntcp.F       newalpha.F \
	cputim.F       patcher.F      qout.F         termn8.F \
	dynptch.F      plot3c.F       resp.F         usrint.F \
	wrest.F        wrestg.F       pre_bc.F       bc_embed.F \
	updateg.F      compg2n.F      resetg.F       bc_period.F \
	yplusout.F     sizer.F        cfl3d.F        trnsfr_vals.F \
	updatedg.F     ae_corr.F      mgbl.F         setslave.F \
	umalloc.F      reass.F        qoutavg.F      plot3davg.F \
        qout_coarse.F  qout_2d.F      plot3d_2d.F    plot3d_coarse.F

FSRC_SPEC =

FOBJ_DIST = $(FSRC_DIST:.F=.o)

FOBJ_SPEC = $(FSRC_SPEC:.F=.o)

DISTLIB = libdist.a

$(DISTLIB): $(FSRC_DIST) $(FOBJ_DIST) $(FSRC_SPEC) $(FOBJ_SPEC)
	ar $(AROPT) $(DISTLIB) $(FOBJ_DIST) $(FOBJ_SPEC)
	@$(RANLIB) $(DISTLIB)

HEAD_DIST = 

$(FOBJ_DIST): $(HEAD_DIST)
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(FOBJ_SPEC): $(HEAD_DIST)
	$(FTN) $(CPPOPT) $(FFLAG_SPEC) -c $*.F

# *************************** CREATE EXECUTABLE ******************************

#Note: for inlining on cray, ccomplex must appear first in fsrc_main

FSRC_MAIN = ccomplex.F development.F main.F

FOBJ_MAIN = $(FSRC_MAIN:.F=.o)

HEAD_MAIN = 

$(FOBJ_MAIN): $(HEAD_MAIN)
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

COMMONLIB   = libcommon.a

$(EXEC): $(FSRC_MAIN) $(FOBJ_MAIN) $(DISTLIB) $(COMMONLIB) 
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC) $(FOBJ_MAIN) \
	$(DISTLIB) $(COMMONLIB) $(LLIBS)
	@ echo "                                                              "
	@ echo "=============================================================="
	@ echo "                                                              "
	@ echo "                   DONE:  $(EXEC) created                     "
	@ echo "                                                              "
	@ echo "           the sequential executable can be found in:         " 
	@ echo "                                                              "
	@ echo "                    $(DIR)/$(EXEC)                            "
	@ echo "                                                              "
	@ echo "=============================================================="
	@ echo "                                                              "

# ******************************* CLEAN/SCRUB ********************************

# the @touch is used to (silently) create some temp files to prevent irksome
# warning messages are sometimes created if there are no *.whatever files and
# one tries to remove them

cleano:
	@touch temp.o
	-rm -f *.o

cleane:
	-rm -f $(EXEC) 

cleana:
	@touch temp.a
	-rm -f *.a 

cleanf:
	@touch temp.f
	-rm -f *.f

cleanh:
	@touch temp.h
	-rm -f *.h

cleang:
	@touch temp.F
	-rm -f *.F

scrub: cleana cleano cleane cleanf cleanh cleang
