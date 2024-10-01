

CXX           = g++
DEFINES       =
CFLAGS        = -pipe -O2 -Wall -Wextra -D_REENTRANT -fPIC $(DEFINES)
CXXFLAGS      = -pipe -O2 -std=c++17 -Wall -Wextra -D_REENTRANT -fPIC $(DEFINES)
INCPATH       = -I.
DISTNAME      = crysmake1.0.0
LINK          = g++
LFLAGS        = -Wl,-O1
LIBS          = $(SUBLIBS) -L/usr/lib64 -lgomp  -lpthread
objdir=objrel

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = main.cpp \
		clcmdlist.cpp \
		databuffer.cpp \
		iofile.cpp \
		stndlstgrain.cpp \
		elements.cpp \
		hcpgrain.cpp 
OBJECTS       = main.o \
		clcmdlist.o \
		databuffer.o \
		iofile.o \
		stndlstgrain.o \
		elements.o \
		hcpgrain.o
TARGET        = crysmake


first: all
####### Build rules

crysmake:  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

all: dirs crysmake


dirs: $(objdir) $(bindir)
	@echo
	@echo -n -e "\e[30;44m COMPILATION ....       \033[0m"
	@echo


$(objdir):
	@echo "create objrel directory"
	mkdir -p $(objdir)


clean: compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Compile

main.o: main.cpp clcmdlist.h \
		iofile.h \
		databuffer.h \
		stndlstgrain.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

clcmdlist.o: clcmdlist.cpp clcmdlist.h \
		iofile.h \
		databuffer.h \
		stndlstgrain.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o clcmdlist.o clcmdlist.cpp

databuffer.o: databuffer.cpp databuffer.h \
		stndlstgrain.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o databuffer.o databuffer.cpp

iofile.o: iofile.cpp iofile.h \
		databuffer.h \
		stndlstgrain.h \
		elements.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o iofile.o iofile.cpp

stndlstgrain.o: stndlstgrain.cpp stndlstgrain.h \
		elements.h \
		clcmdlist.h \
		iofile.h \
		databuffer.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o stndlstgrain.o stndlstgrain.cpp

elements.o: elements.cpp elements.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o elements.o elements.cpp

hcpgrain.o: hcpgrain.cpp hcpgrain.h \
		databuffer.h \
		stndlstgrain.h \
		clcmdlist.h \
		iofile.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o hcpgrain.o hcpgrain.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:

