BIN  = nn

CXX=g++
SHELL=/bin/sh

CPPFLAGS=-O3 -Wall
CFLAGS=$(CPPFLAGS)
LIBS = -lm

SRCS=\
nn.cpp\
mat.cpp\
randmt.cpp

HDRS=\
mat.h\
rand.h

OBJS=\
mat.o\
rand.o

$(BIN): $(OBJS) $(BIN).o
	$(CXX)  $(CFLAGS) $(OBJS) $(BIN).o $(LIBS) -o $(BIN)

$(BIN)oneof: $(OBJS) $(BIN)oneof.o
	$(CXX)  $(CFLAGS) $(OBJS) $(BIN)oneof.o $(LIBS) -o $(BIN)oneof

clean:
	/bin/rm -f *.o a.out
