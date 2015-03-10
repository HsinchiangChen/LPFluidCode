CC     = g++
DEBUG  = -g
OMP    = -fopenmp
#OMP =
INCS   = -I ~/Software/lapack-3.4.2/lapacke/include
LIBS   = -L ~/Software/lapack-3.4.2
CFLAGS = -Wall -c -std=c++11 $(DEBUG) $(OMP)
LFLAGS = -Wall $(DEBUG) $(INCS) $(LIBS) $(OMP)
OBJS   = eos.o geometry.o geometry_collision.o hexagonal_packing.o\
		 initializer.o lp_main.o lp_solver.o ls_solver.o\
         neighbour_searcher.o octree.o particle_data.o\
		 particle_viewer.o registrar.o state.o state_collision.o\
	     time_controller.o

all: lp

lp: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o lp -lgomp -llapacke -llapack -lgfortran -lrefblas

eos.o: eos.h eos.cpp
	$(CC) $(CFLAGS) eos.cpp

geometry.o: geometry.h geometry.cpp
	$(CC) $(CFLAGS) geometry.cpp

geometry_collision.o: geometry.h geometry_collision.h geometry_collision.cpp
	$(CC) $(CFLAGS) geometry_collision.cpp

hexagonal_packing.o: hexagonal_packing.h hexagonal_packing.cpp
	$(CC) $(CFLAGS) hexagonal_packing.cpp

initializer.o: initializer.h initializer.cpp geometry.h state.h eos.h hexagonal_packing.h
	$(CC) $(CFLAGS) initializer.cpp

lp_main.o: lp_main.cpp initializer.h neighbour_searcher.h particle_data.h particle_viewer.h\
           lp_solver.h time_controller.h
	$(CC) $(CFLAGS) lp_main.cpp

lp_solver.o: lp_solver.h lp_solver.cpp neighbour_searcher.h eos.h particle_data.h\
           initializer.h ls_solver.h hexagonal_packing.h
	$(CC) $(CFLAGS) lp_solver.cpp

ls_solver.o: ls_solver.h ls_solver.cpp
	$(CC) $(CFLAGS) ls_solver.cpp

neighbour_searcher.o: neighbour_searcher.h neighbour_searcher.cpp octree.h initializer.h
	$(CC) $(CFLAGS) neighbour_searcher.cpp

octree.o: octree.h octree.cpp
	$(CC) $(CFLAGS) octree.cpp

particle_data.o: particle_data.h particle_data.cpp initializer.h
	$(CC) $(CFLAGS) particle_data.cpp

particle_viewer.o: particle_viewer.h particle_viewer.cpp particle_data.h
	$(CC) $(CFLAGS) particle_viewer.cpp

registrar.o: registrar.h registrar.cpp geometry.h state.h geometry_collision.h state_collision.h
	$(CC) $(CFLAGS) registrar.cpp

state.o: state.h state.cpp
	$(CC) $(CFLAGS) state.cpp

state_collision.o: state.h state_collision.h state_collision.cpp
	$(CC) $(CFLAGS) state_collision.cpp

time_controller.o: time_controller.h time_controller.cpp lp_solver.h particle_viewer.h initializer.h
	$(CC) $(CFLAGS) time_controller.cpp

clean:
	rm *.o *~ lp

clean_out:
	rm -r out

tar:
	tar cvzf lp.tar.gz *.h *.cpp makefile

