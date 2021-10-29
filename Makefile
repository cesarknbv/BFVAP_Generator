#SYSTEM     = x86-64_sles10_4.1
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio127/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio127/concert
#CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio126/cplex
#CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio126/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------
CCC  = g++

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
#CCOPT  = -m32 -O3 -fPIC -fexceptions -DIL_STD
CCOPT  = -m64 -O3 -Wall -Wshadow -fPIC -fexceptions -DIL_STD
#CCOPT  = -m64 -O3 -Wall -fPIC -fexceptions -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -g -lpthread

all:
	make all_cpp

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = .
EXSRC         = .
EXINC         = $(EXDIR)/include
EXDATA        = .

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -std=c++11

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class Comparator
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------

CPP_EX = Comparator

all_cpp: $(CPP_EX)

Comparator: Comparator.o
	$(CCC) $(CCFLAGS) $(CPP_EX).o -o $(CPP_EX) $(CCLNFLAGS)
Comparator.o: $(EXSRC)/$(CPP_EX).cpp
	$(CCC) -c $(CCFLAGS) $(EXSRC)/$(CPP_EX).cpp -o $(CPP_EX).o

