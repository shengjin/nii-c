#####

CC=mpicc


CFLAGS = -Wall #-Werror 
CFLAGS += -Wextra
CFLAGS += -O2 	#-O3

#The error you are seeing: error: ld returned 1 exit status is from the linker ld (part of gcc that combines the object files) because it is unable to find where the function pow is defined.  Including math.h brings in the declaration of the various functions and not their definition. The def is present in the math library libm.a. You need to link your program with this library so that the calls to functions like pow() are resolved.
CFLAGS += -lm 
 
# tune_sigma.c: 
mtmc: main.c mpi_init.c mpi_batch.c mpi_stack.c mpi_flow.c mpi_tune.c data_loader.c alloc.c rand_func.c readin.c user_logll.c user_prior.c 
	$(CC) main.c mpi_init.c mpi_batch.c mpi_stack.c mpi_flow.c mpi_tune.c data_loader.c alloc.c rand_func.c readin.c user_logll.c user_prior.c $(CFLAGS) -o a.out

clean:
	rm -rf *.o
	rm a.out

