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

#include "bdats_h5reader.h"

void  print_help(){
  char *msg="Usage: %s [OPTION] \n\
      	  -h help (--help)\n\
          -f name of the file (only HDF5 file in current version) \n\
          -o output file name (optional: used for binary file) \n\
          -d name of the dataset to read \n  \
          example: dbscan_read -f testf.h5p  \n  \
          example: dbscan_read -f testf.h5p  -d /Step#0/Energy -d /Step#0/x   -d /Step#0/y  -d /Step#0/z -d /Step#0/ux -d /Step#0/uy -d /Step#0/uz \n ";
   fprintf(stdout, msg, "dbscan_read");
}

int main(int argc, char *argv[])
{
	char   filename[OBJ_NAME_MAX], group[OBJ_NAME_MAX], out_filename[OBJ_NAME_MAX];
	int    i, dataset_count = 0;
	int    mpi_size, mpi_rank;

	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Info info = MPI_INFO_NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &mpi_size);
	MPI_Comm_rank(comm, &mpi_rank);


	dataset_info_t  data_info_array[MAX_DATSET_COUNT];
	bool cls = false;
	Points pts;
	Points pts_test;

	static const char *options="f:c:d:o:h";
	extern char *optarg;
	int c;
	while ((c = getopt(argc, argv, options)) != -1) 
	{
		switch (c) 
		{
			case 'f':
			case 'F': 
				strcpy(filename, optarg); 
			case 'o':
			case 'O': 
				strcpy(out_filename, optarg); 
				break;
			case 'c':
			case 'C': 
				strcpy(data_info_array[dataset_count].name, optarg); 
				data_info_array[dataset_count].mpi_rank = mpi_rank;
				data_info_array[dataset_count].mpi_size = mpi_size;
				data_info_array[dataset_count].datatype = H5READER_CLS_INT;
				// data_info_array[dataset_count].datatype = H5READER_CLS_LONG;
				// data_info_array[dataset_count].datatype = H5READER_CLS_FLOAT;
				dataset_count = dataset_count + 1;
				cls = true;
				break;
			case 'd':
			case 'D': 
				strcpy(data_info_array[dataset_count].name, optarg); 
				data_info_array[dataset_count].mpi_rank = mpi_rank;
				data_info_array[dataset_count].mpi_size = mpi_size;
				data_info_array[dataset_count].datatype = H5READER_DATA;
				dataset_count = dataset_count + 1;
			break;
			case 'h':
			case 'H':
				print_help();
			return 1;
			default: 
			break;
		} // switch
	} // while

	float percent_to_read = 10.0;
	bool random = true;
	// Read data: if dataset_count = 0, the the reader will read all columns. 
	read_hdf5_file(filename, &pts, data_info_array, dataset_count, percent_to_read, cls, random); // this will read the class as well
	// read_hdf5_file(filename, &pts, data_info_array, dataset_count, percent_to_read, cls); // this will read the class as well
	// read_hdf5_file(filename, &pts_test, data_info_array, dataset_count);


	// write to a binary file
	/*
	#if 1
	ofstream filew (out_filename, ios::out|ios::binary);	
	int64_t num_points = pts.m_num_local_points;
	int64_t dims = pts.m_num_dims;

	filew.write((char*)&num_points, sizeof(int64_t));
        filew.write((char*)&dims, sizeof(int64_t));
	
	for(int64_t i = 0; i < dims; i++)
		filew.write((char*)&pts.m_data[i][0], num_points * sizeof(BDATS_decimal));

	filew.close();

	#endif
	*/
	
	MPI_Finalize();

	return 0;
}
