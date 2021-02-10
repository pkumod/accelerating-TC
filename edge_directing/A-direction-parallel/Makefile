
NVCC = nvcc -arch=sm_35 -lcudadevrt -rdc=true -O3 #-g -G --ptxas-options=-v
CFLAGS = -c -Wall -O6 -g
EXEFLAG = -O6
OPENMP = -Xcompiler -fopenmp

library = -lcudadevrt

Peel: peel.cu buildCSR.cu 
	$(NVCC) -o Peel peel.cu $(OPENMP)
.PHONY : clean
clean:
	rm Peel
