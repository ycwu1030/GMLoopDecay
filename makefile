
LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools/x86_64-Darwin/lib

SRCDIR := src
INCDIR := include
OBJDIR := obj
CXX = $(LT)/../bin/f++
FFLAG = -I$(LT)/../include -I$(INCDIR) `gsl-config --cflags`
FLIBS = -L$(LT) -looptools `gsl-config --libs`

INC = $(wildcard $(INCDIR)/*.h)
SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))

.PRECIOUS: $(OBJ)
.PHONY: clean all MKOBJDIR

all: MKOBJDIR GMFormFactor.x

%.x:%.cpp $(OBJ) 
	$(CXX) $(FFLAG) -o $@ $< $(OBJ) $(FLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INC)
	$(CXX) $(FFLAG) -c $< -o $@

MKOBJDIR:
	@mkdir -p $(OBJDIR)

clean:
	rm -f *.x
	rm -f $(OBJDIR)/*.o

