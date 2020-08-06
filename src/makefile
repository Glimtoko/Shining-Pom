CC = mpic++
#CFLAGS = -g -Wall -Wextra -ansi -pedantic -fsanitize=address -std=c++17
#CFLAGS = -g -std=c++17 -fsanitize=address
CFLAGS = -O3 -std=c++17

OBJECTS = hydro_flux.o hydro_timestep.o hydro_MH2D_Construct.o hydro_update.o initialise.o AmrCorePom.o main.o

# Expected library locations
AMREX_DIR := /home/nick/MyLibs/AMReX2D

INC = -I$(AMREX_DIR)/include -pthread
LIB = -L$(AMREX_DIR)/lib -lamrex -pthread

VPATH=.:./hydro:./mesh

BL_NO_FORT=TRUE

.PHONY: clean

shiningpom: $(OBJECTS)
	$(CC) $(CFLAGS) *.o -o shiningpom -lm $(LIB)

clean:
	-rm *.o
	-rm shiningpom


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

