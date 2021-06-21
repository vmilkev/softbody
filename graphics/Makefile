.RECIPEPREFIX = >

CXX = g++

CPPFLAGS = -O3 -I./include -MMD

LNKFLAGS = -L/usr/local/lib -lglfw -lpthread -lGLEW -lGLU -lGL -lrt -lXrandr -lXxf86vm -lXi -lXinerama -lX11 -ldl

#define the directive for object files
OBJDIR = ./lib
SRCDIR = ./src
BINDIR = ./bin

# define the C source files
SRCS = $(SRCDIR)/glad.cpp $(SRCDIR)/texture.cpp $(SRCDIR)/shader.cpp $(SRCDIR)/camera.cpp \
	$(SRCDIR)/glwindow.cpp $(SRCDIR)/gprimitive.cpp $(SRCDIR)/cube.cpp $(SRCDIR)/box.cpp \
	$(SRCDIR)/plane.cpp $(SRCDIR)/line.cpp $(SRCDIR)/sphere.cpp $(SRCDIR)/visualizer.cpp

# define the C object files 
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

# define the executable file 
MAIN = $(BINDIR)/liboglu.a

.PHONY: depend 

all: $(MAIN)
>@echo Program compiled

$(MAIN): $(OBJS)
>ar rcs $(MAIN) $(OBJS)
#>$(CXX) $(CPPFLAGS) $(LNKFLAGS) -o $(MAIN) $(OBJS)
#    

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
>@echo "Compiling: " $@
>$(CXX) $(CPPFLAGS) -c -o $@ $<

clean:
>$(RM) $(OBJDIR)/*.o $(MAIN)

depend: $(SRCS)
>makedepend $(INCLUDES) $^
