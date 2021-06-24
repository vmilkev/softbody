.RECIPEPREFIX = >

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)

#CXX = icc
CXX = g++

OS = linux
BINEXMPL = ./example_lin/bin

# CXXFLAGS = -g -O0 -qopenmp -traceback -check-pointers=rw -debug all -debug parallel -mcmodel=large -I${MKLROOT}/include -std=c++17 -traceback -Wall -qoverride-limits
#CXXFLAGS = -O3 -qopenmp -mcmodel=large -I${MKLROOT}/include -std=c++17 -D=$(OS) -MMD
#CPPFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
CXXFLAGS = -O3 -fopenmp -I../graphics/include -std=c++17 -D $(OS) -MMD
LNKFLAGS = -L/usr/local/lib -lglfw -lpthread -lGLEW -lGLU -lGL -lrt -lXrandr -lXxf86vm -lXi -lXinerama -lX11 -ldl -L../graphics/bin -loglu

OBJDIR = ./lib
SRCDIR = ./src
BINDIR = ./bin
INCLUDES = ./incl
EXMPL = ./example_lin

SRCS = $(SRCDIR)/dem.cpp $(SRCDIR)/Functions.cpp
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
MAIN = $(BINDIR)/dem.exe

.PHONY: depend 

all: $(MAIN)
>@echo Program compiled

$(MAIN): $(OBJS)
>$(CXX) $(CXXFLAGS) -o $(MAIN) $(OBJS) $(LNKFLAGS)
>cp $(MAIN) $(BINEXMPL)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
>@echo "Compiling: " $@
>$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
>$(RM) $(OBJDIR)/*.o $(OBJDIR)/*.d $(MAIN)
>rm -rf $(EXMPL)/bin/*.exe $(EXMPL)/sample $(EXMPL)/simulation $(EXMPL)/*.txt

depend: $(SRCS)
>makedepend $(INCLUDES) $^

else

CXX = g++
OS = _Win64
BINEXMPL = example_win\bin

CXXFLAGS = -O3 -fopenmp -std=c++17 -D $(OS) -MMD
CPPFLAGS =

OBJDIR = lib
SRCDIR = src
BINDIR = bin

SRCS = $(SRCDIR)\dem.cpp $(SRCDIR)\Functions.cpp
#OBJS = $(patsubst $(SRCDIR)\%.cpp,$(OBJDIR)\%.o,$(SRCS))
OBJS = lib\dem.o lib\Functions.o
MAIN = $(BINDIR)\dem.exe

all: $(MAIN)
>@echo Program compiled

$(MAIN): $(OBJS)
>$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(MAIN) $(OBJS)
>copy $(MAIN) $(BINEXMPL)

lib\dem.o: src\dem.cpp
>@echo "Compiling: " $@
>$(CXX) $(CXXFLAGS) -c -o $@ $<

lib\Functions.o: src\Functions.cpp
>@echo "Compiling: " $@
>$(CXX) $(CXXFLAGS) -c -o $@ $<

#$(OBJDIR)\%.obj: $(SRCDIR)\%.cpp
#lib\%.o: src\%.cpp
#>@echo "Compiling: " $@
#>$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
>del $(OBJDIR)\*.o $(OBJDIR)\*.d $(MAIN)

endif
