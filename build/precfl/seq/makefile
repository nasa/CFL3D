#     $Id$
#=============================================================================
#
#                     builds the executable for precfl3d
#
#=============================================================================

# ***************************** CREATE LINKS *********************************

link: lncode lnhead

lncode:
	@ echo "        linking source code"
	ln -s $(PRESRC)/*.F .
	ln -s $(CFLSRC_D)/sizer.F .
	ln -s $(CFLSRC_D)/pointers.F .
	ln -s $(CFLSRC_D)/compg2n.F .
	ln -s $(CFLSRC_D)/umalloc.F .
	ln -s $(CFLSRC_D)/termn8.F .
	ln -s $(CFLSRC_S)/outbuf.F .
	ln -s $(CFLSRC_S)/pre_patch.F .
	ln -s $(CFLSRC_S)/pre_blockbc.F .
	ln -s $(CFLSRC_S)/pre_period.F .
	ln -s $(CFLSRC_S)/getdhdr.F .
	ln -s $(CFLSRC_S)/pre_embed.F .
	ln -s $(CFLSRC_S)/global.F .
	ln -s $(CFLSRC_S)/global2.F .
	ln -s $(CFLSRC_S)/rpatch0.F .
	ln -s $(CFLSRC_S)/getibk0.F .
	ln -s $(CFLSRC_S)/cntsurf.F .
	ln -s $(CFLSRC_S)/lead.F .
	ln -s $(CFLSRC_S)/global0.F .
	ln -s $(CFLSRC_S)/readkey.F .
	ln -s $(CFLSRC_S)/parser.F .
	ln -s $(CFLSRC_S)/ccomplex.F .
	ln -s $(CFLSRC_S)/cgnstools.F .
	ln -s $(CFLSRC_S)/setseg.F .
	ln -s $(CFLSRC_S)/my_flush.F .

lnhead:

# ****************************** SUFFIX RULES ********************************

.F.o:
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

# **************************** CREATE EXECUTABLE *****************************

SOURCE = sizer.F       get_cmd_args.F   rpatch0.F        outbuf.F \
	 getibk0.F     global.F         global2.F        pointers.F \
	 compg2n.F     umalloc.F        cntsurf.F \
	 lead.F        pre_patch.F      pre_blockbc.F    pre_embed.F \
	 pre_period.F  usrint.F         termn8.F         getdhdr.F \
	 global0.F     readkey.F        parser.F         ccomplex.F \
	 cgnstools.F   setseg.F         my_flush.F       main.F

FSRC_SPEC =

OBJECT = $(SOURCE:.F=.o)

FOBJ_SPEC = $(FSRC_SPEC:.F=.o)

PRE_HEAD  = 

$(OBJECT): $(PRE_HEAD)
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(FOBJ_SPEC): $(PRE_HEAD)
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
