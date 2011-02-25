#     $Id$
#=============================================================================
#
#             builds the executable for the complex block splitter
#
#=============================================================================

# ***************************** CREATE LINKS *********************************

link: lncode lnhead

lncode:
	@ echo "        linking source code"
	ln -s $(SPLITSRC)/*.F .
	ln -s $(CFLSRC_D)/umalloc.F .
	ln -s $(CFLSRC_D)/umalloc_c.F .
	ln -s $(CFLSRC_S)/parser.F .
	ln -s $(CFLSRC_S)/readkey.F .
	ln -s $(CFLSRC_S)/outbuf.F .
	ln -s $(CFLSRC_S)/my_flush.F .

lnhead:

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# **************************** CREATE EXECUTABLE *****************************

SOURCE = cfl3d_to_tlns3d.F  mapblk.F      movseg.F        spltsg.F \
	 getdhdr.F          mapin.F       outgrd_cmplx.F  termn8.F \
	 ipadd.F            mapout.F      renmbr.F        tlns3d_to_cfl3d.F \
	 ipfix.F            mapsplt.F     shortinp.F      splitter_cmplx.F \
	 iprvsg.F           mkintr.F      sizer.F         umalloc.F \
	 mkseg.F            splt2.F       ingrd_cmplx.F \
	 umalloc_c.F        readkey.F     parser.F        outbuf.F \
	 my_flush.F         main_cmplx.F

FSRC_SPEC =

OBJECT = $(SOURCE:.F=.o)

FOBJ_SPEC = $(FSRC_SPEC:.F=.o)

SPLIT_HEAD  = 

$(OBJECT): $(SPLIT_HEAD)
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(FOBJ_SPEC): $(SPLIT_HEAD)
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
