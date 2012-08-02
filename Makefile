include common/Makefile.common

##############################################################

COMMONDIR = common
INCLUDE = -I common/ -I thirdparty/ -I . 

# if you have the any of the solvers Gurobi, Cplex or Xpress, please add -DHAS_GUROBI etc. to the INCLUDE options

all:  optlib.debug optlib.opt

optlib.debug: $(DEBUGDIR)/check_submodularity.o $(DEBUGDIR)/nbest.o $(DEBUGDIR)/factorMPBP.o $(DEBUGDIR)/factorDualOpt.o $(DEBUGDIR)/factorChainDualDecomp.o $(DEBUGDIR)/factorTRWS.o $(DEBUGDIR)/function.o $(DEBUGDIR)/submodularity_check.o
	ar rs $@ $(DEBUGDIR)/check_submodularity.o $(DEBUGDIR)/nbest.o $(DEBUGDIR)/factorMPBP.o $(DEBUGDIR)/factorDualOpt.o $(DEBUGDIR)/factorChainDualDecomp.o $(DEBUGDIR)/factorTRWS.o $(DEBUGDIR)/function.o $(DEBUGDIR)/submodularity_check.o

optlib.opt: $(OPTDIR)/check_submodularity.o $(OPTDIR)/nbest.o $(OPTDIR)/factorMPBP.o $(OPTDIR)/factorDualOpt.o $(OPTDIR)/factorChainDualDecomp.o $(OPTDIR)/factorTRWS.o $(OPTDIR)/function.o $(OPTDIR)/submodularity_check.o
	ar rs $@ $(OPTDIR)/check_submodularity.o $(OPTDIR)/nbest.o $(OPTDIR)/factorMPBP.o $(OPTDIR)/factorDualOpt.o $(OPTDIR)/factorChainDualDecomp.o $(OPTDIR)/factorTRWS.o $(OPTDIR)/function.o $(OPTDIR)/submodularity_check.o


include common/Makefile.finish
