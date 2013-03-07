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
	ln -s $(TOOLSSRC)/*.F90 .
	ln -s $(CFLSRC_D)/umalloc_r.F .
	ln -s $(CFLSRC_D)/umalloc_c.F .
	ln -s $(CFLSRC_S)/readkey.F .
	ln -s $(CFLSRC_S)/global0.F .
	ln -s $(CFLSRC_S)/parser.F .
	ln -s $(CFLSRC_S)/outbuf.F .
	ln -s $(CFLSRC_S)/cgnstools.F  .
	ln -s  $(CFLSRC_S)/my_flush.F .
	ln -s $(CFLSRC_D)/termn8.F .

lnhead:

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# **************************** CREATE EXECUTABLES ****************************

SOURCELIBS = umalloc_c.F umalloc_r.F \
	readkey.F global0.F cgnstools.F termn8.F parser.F outbuf.F my_flush.F

FSRC_SPEC =

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
SOURCE15 = cgns_to_cfl3dinput.F
SOURCE16 = cgns_readhist.F
SOURCE17 = v6_ronnie_mod.F90
SOURCE18 = cfl3d_to_nmf.F
SOURCE19 = cfl3dinp_to_FVBND.F
SOURCE20 = nmf_to_cfl3dinput.F
SOURCE21 = gridswitchijk.F
SOURCE22 = v6inpswitchijk.F

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
OBJECT15 = $(SOURCE15:.F=.o)
OBJECT16 = $(SOURCE16:.F=.o)
OBJECT17 = $(SOURCE17:.F90=.o)
OBJECT18 = $(SOURCE18:.F=.o)
OBJECT19 = $(SOURCE19:.F=.o)
OBJECT20 = $(SOURCE20:.F=.o)
OBJECT21 = $(SOURCE21:.F=.o)
OBJECT22 = $(SOURCE22:.F=.o)

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
$(OBJECT15):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT16):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT17):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F90
$(OBJECT18):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT19):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT20):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT21):
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
$(OBJECT22):
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
EXEC15 = cgns_to_cfl3dinput
EXEC16 = cgns_readhist
EXEC17 = v6_ronnie_mod
EXEC18 = cfl3d_to_nmf
EXEC19 = cfl3dinp_to_FVBND
EXEC20 = nmf_to_cfl3dinput
EXEC21 = gridswitchijk
EXEC22 = v6inpswitchijk

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

$(EXEC15): $(SOURCE15)  $(OBJECT15) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC15) $(OBJECT15) $(LIBS) $(LLIBS)

$(EXEC16): $(SOURCE16)  $(OBJECT16) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC16) $(OBJECT16) $(LIBS) $(LLIBS)

$(EXEC17): $(SOURCE17)  $(OBJECT17) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC17) $(OBJECT17) $(LIBS) $(LLIBS)

$(EXEC18): $(SOURCE18)  $(OBJECT18) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC18) $(OBJECT18) $(LIBS) $(LLIBS)

$(EXEC19): $(SOURCE19)  $(OBJECT19) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC19) $(OBJECT19) $(LIBS) $(LLIBS)

$(EXEC20): $(SOURCE20)  $(OBJECT20) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC20) $(OBJECT20) $(LIBS) $(LLIBS)

$(EXEC21): $(SOURCE21)  $(OBJECT21) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC21) $(OBJECT21) $(LIBS) $(LLIBS)

$(EXEC22): $(SOURCE22)  $(OBJECT22) $(LIBS)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC22) $(OBJECT22) $(LIBS) $(LLIBS)

$(EXEC): $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) $(EXEC5) $(EXEC6) $(EXEC7) \
	 $(EXEC8) $(EXEC9) $(EXEC10) $(EXEC11) $(EXEC12) $(EXEC13) $(EXEC14) \
         $(EXEC15) $(EXEC16) $(EXEC17) $(EXEC18) $(EXEC19) $(EXEC20) \
         $(EXEC21) $(EXEC22)
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
	@ echo "                      $(DIR)/$(EXEC15)                        "
	@ echo "                      $(DIR)/$(EXEC16)                        "
	@ echo "                      $(DIR)/$(EXEC17)                        "
	@ echo "                      $(DIR)/$(EXEC18)                        "
	@ echo "                      $(DIR)/$(EXEC19)                        "
	@ echo "                      $(DIR)/$(EXEC20)                        "
	@ echo "                      $(DIR)/$(EXEC21)                        "
	@ echo "                      $(DIR)/$(EXEC22)                        "
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
	 $(EXEC8) $(EXEC9) $(EXEC10) $(EXEC11) $(EXEC12) $(EXEC13) $(EXEC14) \
         $(EXEC15) $(EXEC16) $(EXEC17) $(EXEC18) $(EXEC19) $(EXEC20) \
         $(EXEC21) $(EXEC22)

cleanf:
	@touch temp.f
	-rm -f *.f *.F

cleang:
	@touch temp.F
	-rm -f *.F *.F90

cleanh:
	@touch temp.h
	-rm -f *.h

scrub: cleana cleano cleane cleanf cleang cleanh 
