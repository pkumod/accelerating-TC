
NVCC = nvcc -arch=sm_35 -lcudadevrt -rdc=true #-G -g # --ptxas-options=-v
CFLAGS = -c -Wall -O6 -g
EXEFLAG = -O6
OPENMP = -Xcompiler -fopenmp

library = -lcudadevrt

TC: tc.cu preTC.cu 
	$(NVCC) -o TC tc.cu $(OPENMP)
.PHONY : clean
clean:
	rm TC
