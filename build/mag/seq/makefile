#     $Id$
#=============================================================================
#
#                    builds the executable for maggie
#
#=============================================================================

# ***************************** CREATE LINKS **********************************

link: lncode lnhead 

lncode:
	@ echo "         linking source code"
	ln -s $(MAGSRC)/*.F .

lnhead:
	@ echo "         linking header files"
	ln -s $(HEADER)/mag1.h .

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# **************************** CREATE EXECUTABLE *****************************

SOURCE = maggie.F cputim.F wkstn.F

OBJECT = $(SOURCE:.F=.o)

MAG_HEAD  = mag1.h

$(OBJECT): $(MAG_HEAD)
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
	@ echo "                                                              "

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
