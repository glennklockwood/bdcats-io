# HDF5_ROOT = ... #(provided by module load cray-hdf5-parallel)
H5PART_ROOT=$(HOME)/apps.cori-knl/H5Part-1.6.6-intel-18API

CXX       = CC
CXXFLAGS  = -DPARALLEL_IO -I$(H5PART_ROOT)/include -I$(HDF5_ROOT)/include
LDFLAGS   = -L$(H5PART_ROOT)/lib  -L$(HDF5_ROOT)/lib 
LDLIBS    = -lhdf5 -lH5Part -lstdc++

.PHONY: all clean

BINARIES = dbscan_read dbscan_read_dyn
OBJECTS = dbscan_read.o bdats_h5reader.o

all: $(BINARIES)

dbscan_read: $(OBJECTS)

dbscan_read_dyn: LDFLAGS += -dynamic
dbscan_read_dyn: $(OBJECTS)
	$(CC) $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	rm -f $(OBJECTS) $(BINARIES)
