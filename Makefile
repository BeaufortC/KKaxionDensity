ROOTCFLAGS	:= -std=c++1y -m64 -I$(ROOTSYS)/include

ROOTLIBS	:= -L$(ROOTSYS)/lib -lCore -lHist -lGraf -lGpad -lMathCore

CXX = g++
CXXFLAGS     += $(ROOTCFLAGS) 
CINT         = rootcint

LIBS	     = $(ROOTLIBS) $(SYSLIBS) 
LIBS	     += -lgsl -lgslcblas -lm -pthread 

FILESRC = $(wildcard *.cpp)
OBJDIR	= obj
FILEOBJ = $(addprefix $(OBJDIR)/, $(patsubst %.cpp,%.o,$(FILESRC) ))
FILEINC = $(wildcard *.h)

EXESRC =  $(wildcard *.cxx)
EXEOBJ =  $(patsubst %.cxx,%.bin,$(EXESRC) )

default: clean all

all: $(FILEOBJ) $(EXEOBJ)

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/%.o: %.cpp %.h $(OBJSRC)
	@echo -e "" "\n-> Building object: " $@
	$(CXX) -fPIC $(CXXFLAGS) -o $@ -c $<

%.bin: %.cxx $(FILEOBJ)
	@echo -e "" "\n-> Building executable: " $@
	$(CXX) $(CXXFLAGS) -o $@ $? $(LIBS)

all :	 
	@echo -e "" "\n\n------Compilation OK------\n\n" 

clean :
	@echo -e "" "\n-> Deleting compiled files"
	@rm -f $(OBJDIR)/*.o 
	@rm -f *.o
	@rm -f *.bin
