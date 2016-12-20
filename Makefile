CC = CC
SRCS = $(wildcard *.cpp)
HDRS = $(wildcard *.h)
OBJS = $(SRCS:.cpp=.o)
DIRS = $(subst /, ,$(CURDIR))
PROJ = dbscan_reader 

# HDF5 related
H5Proot=/global/homes/s/sbyna/software/h5part
HDF5root = /global/homes/s/sbyna/software/hdf5
H5CFLAGS = -m64 -DUSE_V4_SSE -DOMPI_SKIP_MPICXX
H5PFLAGS = -I${HDF5root}/include -I${H5Proot}/include
H5LIB = -L. -lm -ldl
H5PLIB = -L${HDF5root}/lib -L${H5Proot}/lib -lH5Part -lhdf5 -lz
H5ADD_FLAGS = -DPARALLEL_IO

APP = $(PROJ)
CFLAGS=-c -w 
LDFLAGS= -O3 -openmp 
LIBS=

all: $(APP)

$(APP): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(APP) $(LIBS) -I. $(H5LIB) $(H5PLIB) $(H5ADD_FLAGS) -ldl

%.o: %.cpp $(HDRS) $(MF)
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ $(H5CFLAGS) $(H5PFLAGS) $(H5ADD_FLAGS)

clean:
	rm -f *.o $(APP)

