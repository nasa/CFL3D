#     $Id$
#=============================================================================
#
#                   builds the executable for preronnie
#
#=============================================================================

# ***************************** CREATE LINKS *********************************

link: lncode lnhead 

lncode:
	@ echo "        linking source code"
	ln -s  $(PRERONSRC)/main_pre.F .
	ln -s  $(PRERONSRC)/termn8.F .
	ln -s  $(PRERONSRC)/usrint.F .
	ln -s  $(PRERONSRC)/sizer.F .
	ln -s  $(CFLSRC_D)/umalloc.F .
	ln -s  $(CFLSRC_S)/global2.F .
	ln -s  $(CFLSRC_S)/outbuf.F .
	ln -s  $(CFLSRC_S)/ccomplex.F .
	ln -s  $(CFLSRC_S)/my_flush.F .

lnhead:

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# *************************** CREATE EXECUTABLE ******************************

SOURCE = global2.F    termn8.F    usrint.F    umalloc.F \
	 sizer.F      outbuf.F    ccomplex.F  my_flush.F \
	 main_pre.F 

FSRC_SPEC =

OBJECT = $(SOURCE:.F=.o)

FOBJ_SPEC = $(FSRC_SPEC:.F=.o)

PRERON_HEAD  = 

$(OBJECT): $(PRERON_HEAD)
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(FOBJ_SPEC): $(PRERON_HEAD)
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

# the @touch is used to (silently) create some temp files so that irksome
# warning messages are not created on SGI's if there are no *.whatever files

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
