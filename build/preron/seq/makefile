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
	ln -s  $(CFLSRC_S)/umalloc.F .
	ln -s  $(CFLSRC_S)/ufree.F .
	ln -s  $(CFLSRC_S)/global2.F .
	ln -s  $(CFLSRC_S)/outbuf.F .
	ln -s  $(CFLSRC_S)/ccomplex.F .

lnhead:

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# *************************** CREATE EXECUTABLE ******************************

SOURCE = global2.F    termn8.F    usrint.F    umalloc.F   ufree.F \
	 sizer.F      outbuf.F    ccomplex.F  main_pre.F 

OBJECT = $(SOURCE:.F=.o)

PRERON_HEAD  = 

$(OBJECT): $(PRERON_HEAD)
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(EXEC): $(SOURCE) $(OBJECT)
	$(FTN) $(CPPOPT) $(LFLAG) -o $(EXEC) $(OBJECT) $(LLIBS)
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
