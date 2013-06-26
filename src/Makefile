# Makefile for pdb parser

# Necessary for using GNU g++ compiler
CC = g++
#CFLAGS = -g -pg -Wall 
CFLAGS = -g -Wall 
#CFLAGS = -O3 
#CFLAGS = -O3 -funroll-loops -D_USE_GPLUSPLUS

ECHO = echo

HOME = .
SRC = $(HOME)
INC = $(SRC)/include
EXECDIR = ..

CORE = $(SRC)/utils.cpp $(SRC)/stringutils.cpp $(SRC)/vector3.cpp $(SRC)/matrix3.cpp $(SRC)/matrix.cpp $(SRC)/random.cpp $(SRC)/config.cpp $(SRC)/atom.cpp $(SRC)/residue.cpp $(SRC)/chain.cpp $(SRC)/system.cpp $(SRC)/pairint.cpp $(SRC)/ljint.cpp $(SRC)/hardsphereint.cpp $(SRC)/gentypeint.cpp $(SRC)/kbtypeint.cpp $(SRC)/subset.cpp $(SRC)/forcefield.cpp $(SRC)/movebase.cpp $(SRC)/rotate.cpp $(SRC)/translate.cpp $(SRC)/move.cpp $(SRC)/landscape.cpp

###$(SRC)/.cpp 

OBJECTS = $(SRC)/utils.o $(SRC)/stringutils.o $(SRC)/vector3.o $(SRC)/matrix3.o $(SRC)/matrix.o $(SRC)/random.o $(SRC)/config.o $(SRC)/atom.o $(SRC)/residue.o $(SRC)/chain.o $(SRC)/system.o $(SRC)/pairint.o $(SRC)/ljint.o $(SRC)/hardsphereint.o $(SRC)/subset.o $(SRC)/gentypeint.o $(SRC)/kbtypeint.o $(SRC)/forcefield.o $(SRC)/movebase.o $(SRC)/rotate.o $(SRC)/translate.o $(SRC)/move.o $(SRC)/landscape.o

#implicit rules for handling src files (DEC/COMPAQ/HP- and IBM-style)
.SUFFIXES: 
.SUFFIXES: .o .cpp .h
.cpp.o:
	$(CC) $(CFLAGS) -c $<

#implicit rules for handling src files (GNU-style)
%.o : %.cpp %.h
	$(CC) $(CFLAGS) -c $>

DEP = $(INC)/vector3.h $(INC)/atom.h $(INC)/residue.h $(INC)/chain.h $(INC)/system.h $(INC)/config.h $(INC)/subset.h $(INC)/ffbase.h $(INC)/kbtypeint.h $(INC)/ljint.h $(INC)/hardsphereint.h $(INC)/mtxintbase.h $(INC)/pairint.h $(INC)/forcefield.h $(INC)/movebase.h $(INC)/rotate.h $(INC)/translate.h $(INC)/move.h $(INC)/landscape.h 

LOC = $(INC)

orient2: $(OBJECTS) $(SRC)/orient.cpp $(DEP) $(CFG)
	@$(ECHO) "\n Making orient executable..\n";
	$(CC) $(CFLAGS) -o $(EXECDIR)/orient $(OBJECTS) \
	  $(SRC)/orient.cpp -lm -I$(LOC)
	@$(ECHO) "\n Compilation successful.\n";

orient: $(CORE) $(SRC)/orient.cpp $(DEP) $(CFG)
	  @$(ECHO) "\n Making orient executable..\n";
	  $(CC) $(CFLAGS) -o $(EXECDIR)/orient $(CORE) \
          $(SRC)/orient.cpp -lm -I$(LOC)
	  @$(ECHO) "\n Compilation successful.\n";

perturb2: $(OBJECTS) $(SRC)/perturb.cpp $(DEP) $(CFG)
	@$(ECHO) "\n Making perturb executable..\n";
	$(CC) $(CFLAGS) -o $(EXECDIR)/perturb $(OBJECTS) \
	  $(SRC)/perturb.cpp -lm -I$(LOC)
	@$(ECHO) "\n Compilation successful.\n";

rmsd: $(OBJECTS) $(SRC)/rmsd.cpp $(DEP) $(CFG)
	@$(ECHO) "\n Making rmsd executable..\n";
	$(CC) $(CFLAGS) -o $(EXECDIR)/rmsd $(OBJECTS) \
	  $(SRC)/rmsd.cpp -lm -I$(LOC)
	@$(ECHO) "\n Compilation successful.\n";

interfacepdb2: $(OBJECTS) $(SRC)/interfacepdb.cpp $(DEP) $(CFG)
		@$(ECHO) "\n Making interfacepdb executable..\n";
		$(CC) $(CFLAGS) -o $(EXECDIR)/interfacepdb $(OBJECTS) \
	  	$(SRC)/interfacepdb.cpp -lm -I$(LOC)
		@$(ECHO) "\n Compilation successful.\n";

ppinterface2: $(OBJECTS) $(SRC)/ppinterface.cpp $(DEP) $(CFG)
	@$(ECHO) "\n Making ppinterface executable..\n";
	$(CC) $(CFLAGS) -o $(EXECDIR)/ppinterface $(OBJECTS) \
	  $(SRC)/ppinterface.cpp -lm -I$(LOC)
	@$(ECHO) "\n Compilation successful.\n";

ppinterface: $(CORE) $(SRC)/ppinterface.cpp $(DEP) $(CFG)
	  @$(ECHO) "\n Making ppinterface executable..\n";
	  $(CC) $(CFLAGS) -o $(EXECDIR)/ppinterface $(CORE) \
          $(SRC)/ppinterface.cpp -lm -I$(LOC)
	  @$(ECHO) "\n Compilation successful.\n";

parsepdb2: $(OBJECTS) $(SRC)/parsepdb.cpp $(DEP) $(CFG)
	@$(ECHO) "\n Making parsepdb executable..\n";
	$(CC) $(CFLAGS) -o $(EXECDIR)/parsepdb $(OBJECTS) \
	  $(SRC)/parsepdb.cpp -lm -I$(LOC)
	@$(ECHO) "\n Compilation successful.\n";

parsepdb: $(CORE) $(SRC)/parsepdb.cpp $(DEP) $(CFG)
	  @$(ECHO) "\n Making parsepdb executable..\n";
	  $(CC) $(CFLAGS) -o $(EXECDIR)/parsepdb $(CORE) \
          $(SRC)/parsepdb.cpp -lm -I$(LOC)
	  @$(ECHO) "\n Compilation successful.\n";

evalnrg2: $(OBJECTS) $(SRC)/evalnrg.cpp $(DEP) $(CFG)
	@$(ECHO) "\n Making evalnrg executable..\n";
	$(CC) $(CFLAGS) -o $(EXECDIR)/evalnrg $(OBJECTS) \
	  $(SRC)/evalnrg.cpp -lm -I$(LOC)
	@$(ECHO) "\n Compilation successful.\n";

clean: 
	/bin/rm -f *.o


#dependencies for the object files

config.o       	: include/config.h 
vector3.o      	: include/vector3.h include/config.h
matrix3.o      	: include/matrix3.h include/vector3.h include/config.h
matrix.o       	: include/matrix3.h include/config.h
atom.o         	: include/atom.h vector3.o stringutils.o
residue.o      	: include/residue.h include/config.h atom.o vector3.o
chain.o        	: include/atom.h include/residue.h include/chain.h include/config.h vector3.o residue.o subset.o
system.o       	: include/system.h include/config.h vector3.o chain.o subset.o 
pairint.o      	: include/pairint.h include/ffbase.h
ljint.o        	: include/ljint.h include/ffbase.h
hardsphereint.o : include/hardsphereint.h include/ffbase.h
subset.o       	: include/subset.h
forcefield.o   	: include/forcefield.h include/config.h include/subset.h gentypeint.o
gentypeint.o   	: include/gentypeint.h include/mtxintbase.h pairint.o ljint.o ffbase.o
kbtypeint.o    	: include/kbtypeint.h include/mtxintbase.h pairint.o ljint.o ffbase.o matrix.o
rotate.o       	: include/rotate.h include/config.h movebase.o
translate.o    	: include/translate.h include/config.h movebase.o
move.o         	: include/move.h include/config.h rotate.o translate.o movebase.o
landscape.o    	: include/landscape.h system.o
