#     $Id$
#=============================================================================
#
#                 builds the executable for the block splitter
#
#=============================================================================

# ***************************** CREATE LINKS *********************************

link: lncode lnhead

lncode:
	@ echo "        linking source code"
	ln -s $(SPLITSRC)/*.F .
	ln -s $(CFLSRC_S)/umalloc.F .
	ln -s $(CFLSRC_S)/ufree.F .
	ln -s $(CFLSRC_S)/parser.F .
	ln -s $(CFLSRC_S)/readkey.F .
	ln -s $(CFLSRC_S)/outbuf.F .

lnhead:

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# **************************** CREATE EXECUTABLE *****************************

SOURCE = cfl3d_to_tlns3d.F  mapblk.F      movseg.F     spltsg.F \
	 getdhdr.F          mapin.F       outgrd.F     termn8.F \
	 ipadd.F            mapout.F      renmbr.F     tlns3d_to_cfl3d.F \
	 ipfix.F            mapsplt.F     shortinp.F   splitter.F \
	 iprvsg.F           mkintr.F      sizer.F      umalloc.F \
	 ufree.F            mkseg.F       splt2.F      ingrd.F \
	 readkey.F          parser.F      main.F       outbuf.F 


OBJECT = $(SOURCE:.F=.o)

SPLIT_HEAD  = 

$(OBJECT): $(SPLIT_HEAD)
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(EXEC): $(SOURCE) $(OBJECT)
	$(FTN) $(CPPOPT) $(CPPOPT) $(LFLAG) -o $(EXEC) $(OBJECT) $(LLIBS)
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
