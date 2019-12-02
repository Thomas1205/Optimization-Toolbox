include ../common/Makefile.common

##############################################################

COMMONDIR = ../common/
INCLUDE = -I ../common/ -I ../ -I . 

all:  stereo.debug stereo.opt example.debug

stereo.debug:  $(DEBUGDIR) stereo.cc ../optlib.debug ../common/lib/commonlib.debug  $(DEBUGDIR)/second_order_stereo.o
	$(LINKER) $(DEBUGFLAGS) $(INCLUDE) stereo.cc ../optlib.debug ../common/lib/commonlib.debug ../$(DEBUGDIR)/factorTRWS.o ../$(DEBUGDIR)/factorMPBP.o ../$(DEBUGDIR)/factorDualOpt.o ../$(DEBUGDIR)/factorChainDualDecomp.o ../$(DEBUGDIR)/sepTRWS.o ../$(DEBUGDIR)/separatorDualOpt.o ../$(DEBUGDIR)/separatorChainDualDecomp.o ../common/$(OPTDIR)/vector.o ../common/$(DEBUGDIR)/matrix.o ../common/$(DEBUGDIR)/storage1D.o ../common/$(DEBUGDIR)/storage2D.o $(DEBUGDIR)/second_order_stereo.o -o $@

stereo.opt: $(OPTDIR) stereo.cc ../optlib.opt ../common/lib/commonlib.opt  $(OPTDIR)/second_order_stereo.o
	$(LINKER) $(OPTFLAGS) $(INCLUDE) stereo.cc ../optlib.opt ../common/lib/commonlib.opt ../$(OPTDIR)/factorTRWS.o ../$(OPTDIR)/factorMPBP.o ../$(OPTDIR)/factorDualOpt.o ../$(OPTDIR)/factorChainDualDecomp.o ../$(OPTDIR)/sepTRWS.o ../$(OPTDIR)/separatorDualOpt.o ../$(OPTDIR)/separatorChainDualDecomp.o ../common/$(OPTDIR)/vector.o ../common/$(OPTDIR)/matrix.o ../common/$(OPTDIR)/storage1D.o ../common/$(OPTDIR)/storage2D.o $(OPTDIR)/second_order_stereo.o -o $@

example.debug: $(DEBUGDIR) example.cc ../optlib.debug ../common/lib/commonlib.debug
	$(LINKER) $(DEBUGFLAGS) $(INCLUDE) example.cc ../optlib.debug ../common/lib/commonlib.debug -o $@

include ../common/Makefile.finish
