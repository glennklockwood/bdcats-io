/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */
/*                                                                              */
/*   Files:         Have to list the files                                      */
/*                                                                              */
/*   Description:  HDF5 Reader for BDATS: BigData Analytics at Trillion Scale   */
/*                                                                              */
/*   Author:  Suren Byna		                                        */
/*            Research Scientist, Lawrence Berkeley National Lab    		*/
/*            email: SByna@lbl.gov                              		*/
/*                                                                              */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */

#ifndef _HDF5_READER_
#define _HDF5_READER_

#include "stdlib.h"
#include "hdf5.h"
#include "H5Part.h"
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "bdats_data_struct.h"

#define H5READER_DATA H5T_NATIVE_FLOAT
#define H5READER_CLS H5T_NATIVE_INT
#define H5READER_CLS_INT H5T_NATIVE_INT
#define H5READER_CLS_LONG H5T_NATIVE_LONG
#define H5READER_CLS_FLOAT H5T_NATIVE_FLOAT
#define H5READER_FLOAT H5T_NATIVE_FLOAT
#define H5READER_INT   H5T_NATIVE_INT



#define     MAX_DATSET_COUNT 25         // Increase this based m_num_dims;
#define     OBJ_NAME_MAX 244            // Number of characters in different names
// #define          DEBUG 1             // Uncomment to debug; prints the first 100 data values

// Structure to keep information of a dataset
typedef struct dataset_info{
  hid_t    open_file_id;    //file id returned by open
  char     name[244];       //dataset name, which also contains the group path
  hid_t    open_dataset_id; //dataset id returned by open
  hid_t    open_space_id;   //space id
  hsize_t global_size;     //size of the dataset
  hsize_t local_size;      //local size for the dataset
  hsize_t local_offset;    //offset
  int      space_ndims;
  float    *buf;	   
  int      *cls;
  int      datatype;	   // To find the datatype
  int      mpi_rank;
  int      mpi_size;
}dataset_info_t;


//dataset_info_t  data_info_array[MAX_DATSET_COUNT];


void  open_datasets (dataset_info_t *data_info_array, Points *pts);
void read_segment(dataset_info_t * cur_dataset);
void get_metadata_all_datasets (char* filename, dataset_info_t *data_info_array, Points *pts, int mpi_rank, int mpi_size);
void get_metadata_select_datasets (char* filename, dataset_info_t *data_info_array, Points *pts, int mpi_rank, int mpi_size, int dataset_count);
void allocate_memory (Points *pts, bool cls);
void read_driver (char *filename, Points *pts, MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, dataset_info_t *data_info_array, float percent, bool cls);
void read_driver (char *filename, Points *pts, MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, dataset_info_t *data_info_array, float percent, bool cls, bool random);


int read_hdf5_file(char* filename, Points *pts, dataset_info_t*  data_info_array, int dataset_count, bool cls);
int read_hdf5_file(char* filename, Points *pts, dataset_info_t*  data_info_array, int dataset_count, float percent, bool cls);
int read_hdf5_file(char* filename, Points *pts, dataset_info_t*  data_info_array, int dataset_count, float percent, bool cls, bool random);

// functions to write file
H5PartFile* init_h5outfile (char* filename);
H5PartFile* init_h5outfile_append (char* filename);
int write_dataset_BDATS_lIndex (char* dataset_name, BDATS_lIndex* ids, int64_t num_pts, H5PartFile *H5_out_fileid);
int write_dataset_BDATS_gIndex (char* dataset_name, BDATS_gIndex* ids, int64_t num_pts, H5PartFile *H5_out_fileid);
int write_dataset_BDATS_decimal (char* dataset_name, BDATS_decimal* data, int64_t num_pts, H5PartFile *H5_out_fileid);
int close_h5outfile (H5PartFile *H5_out_fileid);

#endif
