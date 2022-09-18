cgl: cgl.o
        mpicc -std=gnu99 -o cgl -L/software/FFTW/3.3.3-RH7/lib cgl.o -lfftw3_mpi -lfftw3 -lm
cgl.o: cgl.c
        mpicc -std=gnu99 -c -I/software/FFTW/3.3.3-RH7/include cgl.c
clean: 
        rm cgl *.o
