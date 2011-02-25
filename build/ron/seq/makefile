#     $Id$
#=============================================================================
#
#                   builds the executable for ronnie
#
#=============================================================================

# ***************************** CREATE LINKS *********************************

link: lncode lnhead 

lncode:
	@ echo "        linking source code"
	ln -s  $(RONSRC)/*.F .
	ln -s  $(CFLSRC_D)/cputim.F .
	ln -s  $(CFLSRC_D)/writ_buf.F .
	ln -s  $(CFLSRC_D)/umalloc.F .
	ln -s  $(CFLSRC_D)/patcher.F .
	ln -s  $(CFLSRC_S)/lead.F .
	ln -s  $(CFLSRC_S)/pre_patch.F .
	ln -s  $(CFLSRC_S)/collx.F .
	ln -s  $(CFLSRC_S)/global2.F .
	ln -s  $(CFLSRC_S)/arc.F .
	ln -s  $(CFLSRC_S)/avgint.F .
	ln -s  $(CFLSRC_S)/collapse.F .
	ln -s  $(CFLSRC_S)/diagnos.F .
	ln -s  $(CFLSRC_S)/direct.F .
	ln -s  $(CFLSRC_S)/dsmin.F .
	ln -s  $(CFLSRC_S)/expand.F .
	ln -s  $(CFLSRC_S)/extra.F .
	ln -s  $(CFLSRC_S)/extrae.F .
	ln -s  $(CFLSRC_S)/invert.F .
	ln -s  $(CFLSRC_S)/loadgr.F .
	ln -s  $(CFLSRC_S)/newfit.F .
	ln -s  $(CFLSRC_S)/project.F .
	ln -s  $(CFLSRC_S)/rechk.F .
	ln -s  $(CFLSRC_S)/rp3d.F .
	ln -s  $(CFLSRC_S)/shear.F .
	ln -s  $(CFLSRC_S)/topol.F .
	ln -s  $(CFLSRC_S)/topol2.F .
	ln -s  $(CFLSRC_S)/trace.F .
	ln -s  $(CFLSRC_S)/xe.F .
	ln -s  $(CFLSRC_S)/xe2.F .
	ln -s  $(CFLSRC_S)/outbuf.F .
	ln -s  $(CFLSRC_S)/transp.F .
	ln -s  $(CFLSRC_S)/rotatp.F .
	ln -s  $(CFLSRC_S)/ccomplex.F .
	ln -s  $(CFLSRC_S)/my_flush.F .

lnhead:

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# *************************** CREATE EXECUTABLE ******************************

SOURCE = arc.F        diagnos.F    extrae.F     newfit.F     rp3d.F \
 	 avgint.F     direct.F     global2.F     patcher.F    setup.F \
	 collapse.F   dsmin.F      invert.F     project.F    shear.F \
	 expand.F     lead.F       rechk.F      topol.F      usrint.F \
	 collx.F      extra.F      loadgr.F     ronnie.F     topol2.F \
	 trace.F      xe.F         xe2.F        cputim.F     writ_buf.F \
	 termn8.F     outbuf.F     transp.F     rotatp.F     umalloc.F \
	 pre_patch.F  sizer.F      ccomplex.F   my_flush.F   main.F

FSRC_SPEC =

OBJECT = $(SOURCE:.F=.o)

FOBJ_SPEC = $(FSRC_SPEC:.F=.o)

RON_HEAD  =

$(OBJECT): $(RON_HEAD)
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(FOBJ_SPEC): $(RON_HEAD)
	$(FTN) $(CPPOPT) $(FFLAG_SPEC) -c $*.F

$(EXEC): $(SOURCE) $(OBJECT) $(FSRC_SPEC) $(FOBJ_SPEC)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC) $(OBJECT) $(FOBJ_SPEC) $(LLIBS)
	@ echo "                                                              "
	@ echo "=============================================================="
	@ echo "                                                              "
	@ echo "                  DONE:  $(EXEC) created                      "
	@ echo "                                                              "
	@ echo "          the sequential executable can be found in:          "
	@ echo "                                                              "
	@ echo "                      $(DIR)/$(EXEC)                          "
	@ echo "                                                              "
	@ echo "=============================================================="
	@ echo "                                                    "

# ****************************** CLEAN/SCRUB *********************************

# the @touch is used to (silently) create some temp files to prevent irksome
# warning messages are sometimes created if there are no *.whatever files and
# one tries to remove them

cleano:
	@touch temp.o
	-rm -f *.o

cleane:
	-rm -f $(EXEC)

cleanf:
	@touch temp.f
	-rm -f *.f

cleang:
	@touch temp.F
	-rm -f *.F

cleanh:
	@touch temp.h
	-rm -f *.h

scrub: cleano cleane cleanf cleang cleanh 
