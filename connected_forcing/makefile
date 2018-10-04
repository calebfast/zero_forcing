# define variables
HDRDIR  = ./
# set options for this machine
# specify which compilers to use for c and linking
CC	= g++
CPP     = g++
LD	= g++
CLIB     = -L/usr/site/gurobi752/lib -lgurobi75  #-L/opt/gurobi604/linux64/lib -lgurobi60 #
CPPLIB   = $(CLIB) -lgurobi_c++ #-L/opt/gurobi562/linux64/lib -lgurobi_c++ $(CLIB)
# compiler flags to be used (set to compile with debugging on)
CFLAGS = -std=c++11 -I$(HDRDIR) -g

# link flags to be used 
LDFLAGS	= -I$(HDRDIR) -g 

# libraries to be linked in
LIBS	=  -lm $(CLIB) $(CPPLIB)

# types of files we are going to construct rules for
.SUFFIXES: .cpp

# rule for .c files
.cpp.o:
	$(CPP) $(CFLAGS) -o $*.o -c $*.cpp

# list of objects to be compiled
OBJS    = \
	find_connected_forcing_set_MTZ.o\


main:$(OBJS) 

	$(LD)  $(LDFLAGS) -o find_connected_forcing_set_MTZ_setup $(OBJS) $(LIBS)
	rm -r $(OBJS)

# what to do if user types "make clean"
clean :
	rm -r $(OBJS)
