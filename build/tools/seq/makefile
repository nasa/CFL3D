#     $Id$
#=============================================================================
#
#                 builds the executables for some cfl3d utilities
#
#=============================================================================

# ***************************** CREATE LINKS *********************************

link: lncode lnhead

lncode:
	@ echo "        linking source code"
	ln -s $(TOOLSSRC)/*.F .
	ln -s $(CFLSRC_D)/umalloc_r.F .
	ln -s $(CFLSRC_D)/umalloc_c.F .
	ln -s $(CFLSRC_S)/ufree_r.F .
	ln -s $(CFLSRC_S)/ufree_c.F .
	ln -s $(CFLSRC_S)/readkey.F .
	ln -s $(CFLSRC_S)/global0.F .
	ln -s $(CFLSRC_S)/parser.F .
	ln -s $(CFLSRC_S)/outbuf.F .
	ln -s $(CFLSRC_S)/cgnstools.F  .
	ln -s $(CFLSRC_D)/termn8.F .
	ln -s $(CFLSRC_D)/development.F_clean fake.F

lnhead:

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# **************************** CREATE EXECUTABLES ****************************

SOURCELIBS = umalloc_c.F ufree_c.F umalloc_r.F ufree_r.F \
	readkey.F global0.F cgnstools.F termn8.F parser.F outbuf.F

#place routines needing special compilation flags (e.g low-level optimization)
#in FSRC_SPEC...include a fake routine in case there are no others

FSRC_SPEC = fake.F

OBJECTLIBS = $(SOURCELIBS:.F=.o)

FOBJ_SPEC = $(FSRC_SPEC:.F=.o)

LIBS = libs.a

$(LIBS): $(SOURCELIBS) $(OBJECTLIBS) $(FSRC_SPEC) $(FOBJ_SPEC)
	ar $(AROPT) $(LIBS) $(OBJECTLIBS) $(FOBJ_SPEC)
	@$(RANLIB) $(LIBS)

$(OBJECTLIBS):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(FOBJ_SPEC):
	$(FTN) $(CPPOPT) $(FFLAG_SPEC) -c $*.F

SOURCE1  = grid_perturb_cmplx.F 
SOURCE2  = grid_perturb.F
SOURCE3  = Get_FD.F
SOURCE4  = moovmaker.F
SOURCE5  = plot3dg_to_cgns.F
SOURCE6  = XINTOUT_to_ovrlp.F
SOURCE7  = p3d_to_INGRID.F
SOURCE8  = INGRID_to_p3d.F
SOURCE9  = cfl3d_to_pegbc.F
SOURCE10 = everyother_xyz.F
SOURCE11 = initialize_field.F
SOURCE12 = p3d_to_cfl3drst.F
SOURCE13 = v6_restart_mod.F
SOURCE14 = v6inpdoubhalf.F

OBJECT1  = $(SOURCE1:.F=.o)
OBJECT2  = $(SOURCE2:.F=.o)
OBJECT3  = $(SOURCE3:.F=.o)
OBJECT4  = $(SOURCE4:.F=.o)
OBJECT5  = $(SOURCE5:.F=.o)
OBJECT6  = $(SOURCE6:.F=.o)
OBJECT7  = $(SOURCE7:.F=.o)
OBJECT8  = $(SOURCE8:.F=.o)
OBJECT9  = $(SOURCE9:.F=.o)
OBJECT10 = $(SOURCE10:.F=.o)
OBJECT11 = $(SOURCE11:.F=.o)
OBJECT12 = $(SOURCE12:.F=.o)
OBJECT13 = $(SOURCE13:.F=.o)
OBJECT14 = $(SOURCE14:.F=.o)

$(OBJECT1):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT2):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT3):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT4):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT5):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT6):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT7):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT8):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT9):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT10):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT11):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT12):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT13):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT14):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F


EXEC1  = grid_perturb_cmplx
EXEC2  = grid_perturb
EXEC3  = Get_FD
EXEC4  = moovmaker
EXEC5  = plot3dg_to_cgns
EXEC6  = XINTOUT_to_ovrlp
EXEC7  = p3d_to_INGRID
EXEC8  = INGRID_to_p3d
EXEC9  = cfl3d_to_pegbc
EXEC10 = everyother_xyz
EXEC11 = initialize_field
EXEC12 = p3d_to_cfl3drst
EXEC13 = v6_restart_mod
EXEC14 = v6inpdoubhalf

$(EXEC1): $(SOURCE1) $(OBJECT1) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC1) $(OBJECT1) $(LIBS) $(LLIBS)

$(EXEC2): $(SOURCE2) $(OBJECT2) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC2) $(OBJECT2) $(LIBS) $(LLIBS)

$(EXEC3): $(SOURCE3)  $(OBJECT3) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC3) $(OBJECT3) $(LIBS) $(LLIBS)

$(EXEC4): $(SOURCE4)  $(OBJECT4) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC4) $(OBJECT4) $(LIBS) $(LLIBS)

$(EXEC5): $(SOURCE5)  $(OBJECT5) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC5) $(OBJECT5) $(LIBS) $(LLIBS)

$(EXEC6): $(SOURCE6)  $(OBJECT6) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC6) $(OBJECT6) $(LIBS) $(LLIBS)

$(EXEC7): $(SOURCE7)  $(OBJECT7) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC7) $(OBJECT7) $(LIBS) $(LLIBS)

$(EXEC8): $(SOURCE8)  $(OBJECT8) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC8) $(OBJECT8) $(LIBS) $(LLIBS)

$(EXEC9): $(SOURCE9)  $(OBJECT9) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC9) $(OBJECT9) $(LIBS) $(LLIBS)

$(EXEC10): $(SOURCE10)  $(OBJECT10) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC10) $(OBJECT10) $(LIBS) $(LLIBS)

$(EXEC11): $(SOURCE11)  $(OBJECT11) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC11) $(OBJECT11) $(LIBS) $(LLIBS)

$(EXEC12): $(SOURCE12)  $(OBJECT12) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC12) $(OBJECT12) $(LIBS) $(LLIBS)

$(EXEC13): $(SOURCE13)  $(OBJECT13) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC13) $(OBJECT13) $(LIBS) $(LLIBS)

$(EXEC14): $(SOURCE14)  $(OBJECT14) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC14) $(OBJECT14) $(LIBS) $(LLIBS)

$(EXEC): $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) $(EXEC5) $(EXEC6) $(EXEC7) \
	 $(EXEC8) $(EXEC9) $(EXEC10) $(EXEC11) $(EXEC12) $(EXEC13) $(EXEC14)
	@ echo "                                                              "
	@ echo "=============================================================="
	@ echo "                                                              "
	@ echo "                  DONE:  $(EXEC) created                      "
	@ echo "                                                              "
	@ echo "          the sequential executables can be found in:         "
	@ echo "                                                              "
	@ echo "                      $(DIR)/$(EXEC1)                         "
	@ echo "                      $(DIR)/$(EXEC2)                         "
	@ echo "                      $(DIR)/$(EXEC3)                         "
	@ echo "                      $(DIR)/$(EXEC4)                         "
	@ echo "                      $(DIR)/$(EXEC5)                         "
	@ echo "                      $(DIR)/$(EXEC6)                         "
	@ echo "                      $(DIR)/$(EXEC7)                         "
	@ echo "                      $(DIR)/$(EXEC8)                         "
	@ echo "                      $(DIR)/$(EXEC9)                         "
	@ echo "                      $(DIR)/$(EXEC10)                        "
	@ echo "                      $(DIR)/$(EXEC11)                        "
	@ echo "                      $(DIR)/$(EXEC12)                        "
	@ echo "                      $(DIR)/$(EXEC13)                        "
	@ echo "                      $(DIR)/$(EXEC14)                        "
	@ echo "                                                              "
	@ echo "=============================================================="
	@ echo "                                                              "

# ****************************** CLEAN/SCRUB *********************************

# the @touch is used to (silently) create some temp files to prevent irksome
# warning messages are sometimes created if there are no *.whatever files and
# one tries to remove them

cleana:
	@touch temp.a
	-rm -f *.a
cleano:
	@touch temp.o
	-rm -f *.o

cleane:
	-rm -f $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) $(EXEC5) $(EXEC6) $(EXEC7) \
	 $(EXEC8) $(EXEC9) $(EXEC10) $(EXEC11) $(EXEC12) $(EXEC13) $(EXEC14)

cleanf:
	@touch temp.f
	-rm -f *.f

cleang:
	@touch temp.F
	-rm -f *.F

cleanh:
	@touch temp.h
	-rm -f *.h

scrub: cleana cleano cleane cleanf cleang cleanh 
