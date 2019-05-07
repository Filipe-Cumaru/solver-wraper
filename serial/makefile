# Include the MOAB configuration information so that
# all required flags and libs are populated correctly
include makefile.config

default: all

ALLEXAMPLES = TPFASolver
TRILINOS_LIB_LINKS = -DMYAPP_EPETRA -L/home/facsa/Trilinos-Serial/lib -lamesos -laztecoo -ltrilinosss -ltriutils -lepetra -lepetraext -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkoscore -lkokkoscore

all: $(ALLEXAMPLES)

TPFASolver: TPFASolver.o
	@echo "[CXXLD]  $@"
	${VERBOSE}$(MOAB_CXX) -o $@ main.cpp $< $(MOAB_LIBS_LINK) -std=c++11 -I$PWD -I/usr/include -I/home/facsa/Trilinos-Serial/include $(TRILINOS_LIB_LINKS)

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}
