
# NVCC compiler and flags
NVCC = nvcc
NVCCFLAGS   = -O3 --compiler-options '-fPIC' --compiler-bindir=/usr/bin/gcc -Xcompiler -Wall -arch=sm_61  -lcudart 

# linker options
LFLAGS_CUFFT = -lcufft
LFLAGS_PGPLOT = -L/usr/lib64/pgplot -lpgplot -lcpgplot -lX11


# bin directory
BINDIR = ./bin

CC          = g++
HPDEMO_LIB_CCFLAGS     = -g -O3 -fPIC -shared -lstdc++ -mavx -msse4 \
                     -I. -I$(CUDA_DIR)/include -I/usr/local/include \
                     -L. -L/usr/local/lib \
                     -lhashpipe -lrt -lm

FILTERBANK_OBJECT   = filterbank.o
FILTERBANK_SOURCES  = filterbank.cpp
FILTERBANK_LIB_INCLUDES = filterbank.h

HPDEMO_LIB_TARGET   = crab_hashpipe.so
HPDEMO_LIB_SOURCES  = crab_net_thread.c \
                      crab_gpu_thread.c \
		      crab_output_thread.c\
                      crab_databuf.c
HPDEMO_LIB_INCLUDES = crab_databuf.h \

#GPU_LIB_TARGET = 
#GPU_LIB_SOURCES =
#GPU_LIB_INCLUDES =  
all: $(HPDEMO_LIB_TARGET)

$(FILTERBANK_OBJECT): $(FILTERBANK_SOURCES) $(FILTERBANK_LIB_INCLUDES)
	$(CC) -c $(FILTERBANK_OBJECT) $(FILTERBANK_SOURCES) $(HPDEMO_LIB_CCFLAGS) 
$(HPDEMO_LIB_TARGET): $(HPDEMO_LIB_SOURCES) $(HPDEMO_LIB_INCLUDES)
	$(CC) -o $(HPDEMO_LIB_TARGET) $(HPDEMO_LIB_SOURCES) $(FILTERBANK_OBJECT) $(HPDEMO_LIB_CCFLAGS)
#$(NVCC) $(NVFLAGS) tut5_*.o $<                                            \
#            $(LFLAGS_MATH) $(LFLAGS_CUFFT) $(LFLAGS_PGPLOT) -o $(BINDIR)/$@
#$(FILTERBANK_OBJECT): $(FILTERBANK_SOURCES) $(FILTERBANK_LIB_INCLUDES)
#	$(CC) -c  $(FILTERBANK_SOURCES) $(HPDEMO_LIB_CCFLAGS) 
tags:
	ctags -R .
clean:
	rm -f $(HPDEMO_LIB_TARGET) tags

prefix=/usr/local
LIBDIR=$(prefix)/lib
BINDIR=$(prefix)/bin
install-lib: $(HPDEMO_LIB_TARGET)
	mkdir -p "$(DESTDIR)$(LIBDIR)"
	install -p $^ "$(DESTDIR)$(LIBDIR)"
install: install-lib

.PHONY: all tags clean install install-lib
# vi: set ts=8 noet :
