#
# 			Makefile for the program "ppp"
#			    
#
###############################################################################
#				directories				       
###############################################################################
#
progname=ppp.x
BINDIR=$(PWD)/bin
DIR=$(PWD)
LAPACK=/opt/intel/mkl/lib/intel64
KLROOT=/opt/intel/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=$(MKLROOT)/include/intel64/lp64
mklinc1=$(MKLROOT)/include

FCCFLAG= -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm

#
# ############    Definitons of objectfiles    ###############
#
object= \
$(DIR)/density.o $(DIR)/density_uhf.o $(DIR)/fockmat.o \
$(DIR)/fockmat_uhf.o \
$(DIR)/ijpk.o $(DIR)/matel.o $(DIR)/plblk.o \
$(DIR)/prblk.o $(DIR)/ppp.o $(DIR)/r_atom.o  $(DIR)/r_site.o \
$(DIR)/scf_rhf.o $(DIR)/scf_uhf.o $(DIR)/tij_read.o $(DIR)/vij_cal.o \
$(DIR)/sci.o $(DIR)/twoind.o \
$(DIR)/fourind.o $(DIR)/ham_sci.o $(DIR)/optics.o \
$(DIR)/delatm.o $(DIR)/printr.o $(DIR)/dipmat.o \
$(DIR)/spctrm_a.o $(DIR)/orbden.o  $(DIR)/benzen.o \
$(DIR)/benperp.o $(DIR)/stline.o $(DIR)/symsite.o \
$(DIR)/bond.o $(DIR)/bands.o $(DIR)/foutra_tb.o $(DIR)/tij_readb.o \
$(DIR)/writorb.o $(DIR)/ci_drv.o $(DIR)/ecp.o \
$(DIR)/write_1.o $(DIR)/write_2.o \
$(DIR)/dipint.o $(DIR)/dipout.o $(DIR)/read_ei.o \
$(DIR)/nlo.o $(DIR)/dipcal.o $(DIR)/delo_nlo.o \
$(DIR)/spctrm_1ex.o $(DIR)/tij_gen.o $(DIR)/c60_gen.o
#
############################
# Include files
############################
#
# Variables
#
FCOMP=ifort
#FOPTIONS=-O  
OPTS= -debug -g
#
################################################################################
#                       rules                                                  #
################################################################################
.SUFFIXES: .f90 .f
#.f.o:
#	$(FCOMP) $(FOPTIONS) -c  $< -o $@
#%.o: %.f90
#	$(FCOMP) $(FOPTIONS) -c  $< -o $@
#
#.f90.o:
#%.o: %.f90
#	$(FCOMP) $(OPTS) -c  $< -o $@ 
%.o : %.f90
	$(FCOMP) -c $(OPTS) -I/$(mklinc) -I/$(mklinc1) $<
#
###############################################################################
#			targets						     
###############################################################################
#
ppp.x : $(object)
	$(FCOMP) -o $@ $^ -L$(MKLLIB) $(FCCFLAG)
#.DEFAULT: $(BINDIR)/$(progname)
#$(BINDIR)/$(progname): $(OBJECTS)    
#	@echo  'linking ' $(progname)
#	$(FCOMP) $(OPTS) $(OBJECTS)  -L$(MKLLIB) $(FCCFLAG) -o tmp.tmp
#	mv -f tmp.tmp $(BINDIR)/$(progname)
#	strip $(BINDIR)/$(progname)
#	@echo 'done'
clean:
	@ rm -f $(DIR)/*.o ppp.x
veryclean:
	@ rm -f $(OBJECT) *~
#tags:
#	etags $(DIR)/*.{f,F} $(DIR1)/*.{f,F}  $(LIB)/*.{f,F} $(DIR4)/*.{f,F}
#  
################################################################################
