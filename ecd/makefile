# -*- Mode: makefile -*-
# Copyright (c) 1999 Matthew Wall, all rights reserved
# -----------------------------------------------------------------------------
#   To make all of the examples, do 'make'.  You can compile any one of
# the examples by typing 'make exN' where N is the number of the example you
# want to compile.  See the README for a description of what each example does.
# -----------------------------------------------------------------------------

include makevars

# Set these paths to the location of the GA library and headers.
#GA_INC_DIR= /usr/local/include
#GA_LIB_DIR= /usr/local/lib
OB_SOURCE_DIR= /home/liu/soft/babel/sourceCode/openbabel-2.3.2

#OB_INC_DIR=/home/liu/soft/babel/include/openbabel-2.0
#OB_LIB_DIR= /home/liu/soft/babel/lib
#OB_LIB_DIR2= /home/liu/soft/babel/lib/openbabel/2.3.2

OB_INC_DIR= /home/liu/soft/babel-git/apps/include/openbabel-2.0
OB_LIB_DIR= /home/liu/soft/babel-git/apps/lib
#
GA_SOURCE_DIR= /home/liu/projects/conformation/GA/galib247
#GA_INC_DIR= $(OB_SOURCE_DIR)/include
GA_INC_DIR= $(GA_SOURCE_DIR)
GA_LIB_DIR= $(GA_SOURCE_DIR)/ga

INC_DIRS= -I$(OB_INC_DIR) -I$(GA_INC_DIR)
LIB_DIRS= -lopenbabel -L$(OB_LIB_DIR) -lga -L$(GA_LIB_DIR)


EXS=ecdconf 

.SUFFIXES: .cpp
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INC_DIRS) -c $<

all: $(EXS)

# Use this for non-gnu make
#$(EXS): $$@.o
#	$(CXX) $@.o -o $@ $(LIB_DIRS) -lga -lm $(CXX_LIBS)

# Use this for gnu make
#$(EXS): %: %.o
#	$(CXX) $@.o -o $@ $(LIB_DIRS) -lga -lopenbabel -linchi -lm $(CXX_LIBS)
$(EXS): %: *.o
	$(CXX) $@.o -o $@ $(LIB_DIRS) -lm $(CXX_LIBS)

clean:
	$(RM) $(EXS)
	$(RM) *.o *~ *.bak *.pixie core
	$(RM) test_results.txt test_stats.txt
	$(RM) *.dat
	$(RM) *.out *.exe vc* *.pdb
