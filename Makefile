#!/bin/bash

CC=g++
CFLAGS= -ansi -O5 -Wall
LDFLAGS= -ansi -lm -Wall
EXEC=bin/louvain_community bin/louvain_convert bin/louvain_hierarchy
OBJ1= louvain/graph_binary.o louvain/community.o
OBJ2= louvain/graph.o


all: $(EXEC)

bin/louvain_community : $(OBJ1) louvain/main_community.o
	$(CC) -o $@ $^ $(LDFLAGS)

bin/louvain_convert : $(OBJ2) louvain/main_convert.o
	$(CC) -o $@ $^ $(LDFLAGS)

bin/louvain_hierarchy : louvain/main_hierarchy.o
	$(CC) -o $@ $^ $(LDFLAGS)


##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CC) -o $@ -c $< $(CFLAGS)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -f *.o *~ $(EXEC)
