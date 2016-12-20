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

void  open_datasets (dataset_info_t *data_info_array, Points *pts)
{
        hid_t   dset_id;
        hsize_t proc_offset;
        int     space_ndims;
        hid_t  dataspace;
        dset_id  = H5Dopen(data_info_array->open_file_id, data_info_array->name, H5P_DEFAULT);
	if (dset_id < 0)
	{
		printf ("Error opening dataset: %s, in open_datasets [H5Dopen]", data_info_array->name);
		exit(1);
	}
	
        dataspace = H5Dget_space(dset_id);
	if (dataspace < 0)
	{
		printf ("Error retrieving dataspace for dataset: %s, in open_datasets [H5Dget_space]", data_info_array->name);
		exit(1);
	}
	
        space_ndims = H5Sget_simple_extent_ndims(dataspace);
	if (space_ndims < 0)
	{
		printf ("Error retrieving ndims for dataset: %s, in open_datasets [H5Sget_simple_extent_ndims]", data_info_array->name);
		exit(1);
	}
	

        proc_offset = data_info_array->mpi_rank * (pts->m_num_global_points/data_info_array->mpi_size);

        data_info_array->open_dataset_id = dset_id;
        data_info_array->local_offset    = proc_offset;
        data_info_array->global_size     = pts->m_num_global_points;
        data_info_array->local_size      = pts->m_num_local_points;
        data_info_array->open_space_id   = dataspace;
        data_info_array->space_ndims     = space_ndims;
}

void read_segment(dataset_info_t * cur_dataset)
{
        hid_t plist2_id, memspace;
	int retval;

        memspace =  H5Screate_simple(cur_dataset->space_ndims, &(cur_dataset->local_size), NULL);
	if (memspace < 0)
	{
		printf ("Error creating dataspace for dataset: %s, in read_segment [H5Screate_simple]", cur_dataset->name);
		exit(1);
	}

        plist2_id = H5Pcreate(H5P_DATASET_XFER);
	if (plist2_id < 0)
	{
		printf ("Error in creating properties for dataset: %s, in read_segment [H5Pcreate]", cur_dataset->name);
		exit(1);
	}
	
	// Set MPI-IO to be collective I/O (two-phase I/O) mode
        retval = H5Pset_dxpl_mpio(plist2_id, H5FD_MPIO_COLLECTIVE);
	if (retval < 0)
	{
		printf ("Error in setting properties for dataset: %s, in read_segment [H5Pset_dxpl_mpio]", cur_dataset->name);
		exit(1);
	}
	

        hsize_t my_offset, my_size;
        my_offset = cur_dataset->local_offset;
        my_size   = cur_dataset->local_size;
        retval = H5Sselect_hyperslab(cur_dataset->open_space_id, H5S_SELECT_SET, &my_offset, NULL, &my_size , NULL);
	if (retval < 0)
	{
		printf ("Error in selecting hyperslab for dataset: %s, in read_segment [H5Sselect_hyperslab]", cur_dataset->name);
		exit(1);
	}
	
	// If the data has to be copied to pts->m_data
	if (cur_dataset->datatype == H5READER_DATA)
	{
		retval = H5Dread(cur_dataset->open_dataset_id, H5T_NATIVE_FLOAT, memspace, cur_dataset->open_space_id, plist2_id, cur_dataset->buf);
	}
	// If the data has to be copied to pts->m_class and pts->m_class is a H5T_NATIVE_LONG type
	else if (cur_dataset->datatype == H5READER_CLS_LONG)
	{
		retval = H5Dread(cur_dataset->open_dataset_id, H5T_NATIVE_LONG, memspace, cur_dataset->open_space_id, plist2_id, cur_dataset->cls);
	}
	// If the data has to be copied to pts->m_class and pts->m_class is a H5T_NATIVE_INT type
	else if (cur_dataset->datatype == H5READER_CLS_INT)
	{
		retval = H5Dread(cur_dataset->open_dataset_id, H5T_NATIVE_INT, memspace, cur_dataset->open_space_id, plist2_id, cur_dataset->cls);
	}
	else if (cur_dataset->datatype == H5READER_CLS_FLOAT)
	{
		retval = H5Dread(cur_dataset->open_dataset_id, H5T_NATIVE_FLOAT, memspace, cur_dataset->open_space_id, plist2_id, cur_dataset->cls);
	}
	if (retval < 0)
	{
		printf ("Error in reading dataset: %s, in read_segment [H5Dread]", cur_dataset->name);
		exit(1);
	}
	
        retval = H5Pclose(plist2_id);
	if (retval < 0)
	{
		printf ("Error in closing properties for dataset: %s, in read_segment [H5Pclose]", cur_dataset->name);
		exit(1);
	}
	
}

void get_metadata_all_datasets (char* filename, dataset_info_t *data_info_array, Points *pts, int mpi_rank, int mpi_size)
{
        int dataset_count;
        H5PartFile* H5P_fileid;
        H5P_fileid = H5PartOpenFileParallel(filename, H5PART_READ, MPI_COMM_WORLD);
        dataset_count = H5PartGetNumDatasets (H5P_fileid);

        int num_timesteps = H5PartGetNumSteps (H5P_fileid);

        // char step_name[OBJ_NAME_MAX];
        // if (H5PartReadStepAttrib(H5P_fileid, "filename", &step_name[0]) == 1)
        // {
                // printf("Read step from file: %s\n", step_name);
        // }

        // h5part_int64_t num_attribs = H5PartGetNumStepAttribs(H5P_fileid);
        // fprintf(stdout, "Number of step attributes in step #%lld: %lld\n", num_timesteps, (long long)num_attribs);

        int index;
        char   dataset_name[OBJ_NAME_MAX];

        pts->m_num_global_points = H5PartGetNumParticles (H5P_fileid);
        pts->m_num_dims = dataset_count;
        for(index=0; index < dataset_count; index++)
        {
                data_info_array[index].mpi_rank = mpi_rank;
                data_info_array[index].mpi_size = mpi_size;
                H5PartGetDatasetName(H5P_fileid, index, dataset_name, OBJ_NAME_MAX);
                // Replace the hard-coded time step name
                sprintf (data_info_array[index].name, "/Step#0/%s", dataset_name);
                // printf("\tDataset[%u] name= %s\n", index, data_info_array[index].name);
        }

	if (mpi_rank == 0) {
	    cout << "pts->m_num_global_points " << pts->m_num_global_points << " pts->m_num_dims " << pts->m_num_dims << endl;
	}

        H5PartCloseFile(H5P_fileid);
}

void get_metadata_select_datasets (char* filename, dataset_info_t *data_info_array, Points *pts, int mpi_rank, int mpi_size, int dataset_count)
{
        H5PartFile* H5P_fileid;
        H5P_fileid = H5PartOpenFileParallel(filename, H5PART_READ, MPI_COMM_WORLD);

        int num_timesteps = H5PartGetNumSteps (H5P_fileid);

        int index;
        char   dataset_name[OBJ_NAME_MAX];

        pts->m_num_global_points = H5PartGetNumParticles (H5P_fileid);
        pts->m_num_dims = dataset_count;
        for(index=0; index < dataset_count; index++)
        {
                data_info_array[index].mpi_rank = mpi_rank;
                data_info_array[index].mpi_size = mpi_size;
                // printf("\tDataset[%u] name= %s\n", index, data_info_array[index].name);
        }

	if (mpi_rank == 0) {
	    cout << "pts->m_num_global_points " << pts->m_num_global_points << " pts->m_num_dims " << pts->m_num_dims << endl;
	}

        H5PartCloseFile(H5P_fileid);
}

void get_metadata_select_datasets_h5 (char* filename, dataset_info_t *data_info_array, Points *pts, int mpi_rank, int mpi_size, int dataset_count)
{
	hid_t H5_fileid, H5_gid, H5_datasetid, H5_typeid, plist_id;

	plist_id = H5Pcreate (H5P_FILE_ACCESS);
	H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	
	H5_fileid = H5Fopen (filename, H5F_ACC_RDONLY, plist_id);

        int index;

	/*
        pts->m_num_global_points = H5PartGetNumParticles (H5P_fileid);
        pts->m_num_dims = dataset_count;
        for(index=0; index < dataset_count; index++)
        {
                data_info_array[index].mpi_rank = mpi_rank;
                data_info_array[index].mpi_size = mpi_size;
                // printf("\tDataset[%u] name= %s\n", index, data_info_array[index].name);
        }

	if (mpi_rank == 0) {
	    cout << "pts->m_num_global_points " << pts->m_num_global_points << " pts->m_num_dims " << pts->m_num_dims << endl;
	}

        H5PartCloseFile(H5P_fileid);
	*/
}

void allocate_memory (Points *pts, bool cls)
{
        BDATS_decimal** data_array;
        int dim;
	int num_dims;

	// If classification dataset is available, pts->m_data will have one dataset less
	if (cls)
		num_dims = pts->m_num_dims - 1;
	else
		num_dims = pts->m_num_dims;
		
	// allocate memory for pts->m_data 
        data_array = (BDATS_decimal**) BDATS_malloc(num_dims * sizeof (BDATS_decimal*));
        for (dim = 0; dim < num_dims; dim++) {
                data_array[dim] = (BDATS_decimal*) BDATS_malloc (pts->m_num_local_points * sizeof (BDATS_decimal));
                if(data_array[dim] == NULL)
                {
                        printf("Memory allocation for data_array[%d] failed... \n", dim);
                        exit(-1);
                }
        }
	
        pts->m_data = data_array;

	// allocate memory for pts->m_class, if cls is true
	if (cls)
	{
		pts->m_class = (int *) BDATS_malloc (pts->m_num_local_points * sizeof (int));
	}	
}

// Compare function for qsort
int qcompare (const void * a, const void * b)
{
	return ( *(int64_t*)a - *(int64_t*)b );
}

// Assign a random percent of the points in the data read by each process to BDATS Points
void read_driver (char *filename, Points *pts, MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, dataset_info_t *data_info_array, float percent, bool cls, bool random)
{
        dataset_info_t *current_dataset;
        hid_t  file_id, plist_id;
        hsize_t         local_size;

        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, info);
        file_id = H5Fopen(filename, H5F_ACC_RDONLY,plist_id);

        double t1, t2;
        MPI_Barrier(MPI_COMM_WORLD);
        t1 =  MPI_Wtime();

        int i, index;
        int64_t last_process_extra_pts = 0;

        // Assuming that each dataset has the same size
        last_process_extra_pts = pts->m_num_global_points % mpi_size;
        if (mpi_rank ==  (mpi_size - 1))
        {
		// Read all the points first; after reading, select the top percent 
                pts->m_num_local_points = (pts->m_num_global_points/mpi_size + last_process_extra_pts);
		#ifdef _BDATS_DEBUG_
		    printf ("MPI_Rank: %d --> m_num_local_points %lld \n", mpi_rank, pts->m_num_local_points); 
		#endif
        }
        else
        {
		// Read all the points first; after reading, select the top percent 
                pts->m_num_local_points = (pts->m_num_global_points/mpi_size);
		#ifdef _BDATS_DEBUG_
		    printf ("MPI_Rank: %d --> m_num_local_points %lld \n", mpi_rank, pts->m_num_local_points); 
		#endif
        }

        // Allocate memory for data (i.e., pts->m_data and pts->m_class [if classification variable is true]).
        allocate_memory (pts, cls);

        // Go through each HDF5 dataset and read data into certain locations.
        for(i = 0; i < pts->m_num_dims; i++)
        {
                current_dataset = &(data_info_array[i]);
                current_dataset->open_file_id = file_id;

                //Open datasets and set hyperslab boundaries
                open_datasets (current_dataset, pts);

		// Suren: 10/02/2015: Changing float* buf to void *buf; typecasting to float *
		if (current_dataset->datatype == H5READER_DATA)
		{
			current_dataset->buf = *&pts->m_data[i];
			if(current_dataset->buf == NULL)
			{
				printf("Memory allocation for reading data fails ! \n");
				exit(-1);
			}
		}
		else 
		{
                	current_dataset->cls = pts->m_class;
			if(current_dataset->cls == NULL)
			{
				printf("Memory allocation for reading class data fails ! \n");
				exit(-1);
			}
		}

                //Parallel read
                read_segment(current_dataset);

		#ifdef _BDATS_DEBUG_
                if (mpi_rank == 0)
                        printf ("Finished reading dataset[%d] \n", i);
		#endif
	
		#if 0	
		#ifdef _BDATS_DEBUG_
		// Debug the contents of pts->m_class 
		
		if (current_dataset->datatype == H5READER_CLS_INT)
		{
                        for(int ki = 0; ki < 100; ki++)
                        {
                                printf ("%i \t", pts->m_class[ki]);
                        }
			printf ("\n");
		}
		#endif
		#endif
        }

	if (random)
	{
		// The number of random points is a percent of total points on an MPI process
		int64_t num_random_points = pts->m_num_local_points * (percent/100.0);
		#ifdef _BDATS_DEBUG_
			printf ("MPI Rank: %d --> number of random points: %lld \n", mpi_rank, num_random_points);
		#endif

		// Allocate memory for location of random points
		int64_t *random_locations = (int64_t *) malloc (num_random_points * sizeof(int64_t));
		if(random_locations == NULL)
		{
			printf("Memory allocation for random_locations variable failed... \n");
			exit(-1);
		}
		int64_t count_rands = 0;
		int64_t temp_rand;
		bool *random_loc_seen = (bool *) malloc (num_random_points * sizeof(int64_t));

		pts->m_vec_global_IDs =  (BDATS_gIndex*) BDATS_malloc (num_random_points * sizeof(BDATS_gIndex));
		if(pts->m_vec_global_IDs == NULL)
		{
			printf("Memory allocation for m_vec_global_IDs failed... \n");
			exit(-1);
		}

		int64_t global_offset = (pts->m_num_global_points/mpi_size) * mpi_rank; 
		#ifdef _BDATS_DEBUG_
			printf ("MPI Rank: %d --> global_offset: %lld \n", mpi_rank, global_offset);
		#endif
	
		// Generate unique random numbers between 0 and pts->m_num_local_points range.
		for (int64_t r_index = 0; r_index < num_random_points; r_index++)
		{
			random_loc_seen[r_index] = false;
		}
	
		srand((unsigned)time(NULL));
		while (count_rands < num_random_points)
		{
			temp_rand = rand () % num_random_points;
			if (!random_loc_seen[temp_rand])
			{
				random_loc_seen[count_rands] = true;
				random_locations[count_rands] = temp_rand;
				count_rands++;
			}
		}

		// Sort the random locations; otherwise, m_data and m_class 
		// data may have overlaps and duplicates
		// TODO: Any better sorting algorithms here needed?
		qsort (random_locations, num_random_points, sizeof(int64_t), qcompare);

		int64_t temp_rand_loc;
		for (int64_t r_index = 0; r_index < num_random_points; r_index++)
		{
			temp_rand_loc = random_locations[r_index]; 
			pts->m_vec_global_IDs[r_index] = global_offset + temp_rand_loc;
			if (cls)
			{
				for(i = 0; i < pts->m_num_dims-1; i++)
				{
					pts->m_data[i][r_index] = pts->m_data[i][temp_rand_loc];
				}
				pts->m_class[r_index] = pts->m_class[temp_rand_loc];
			}
			else
			{
				for(i = 0; i < pts->m_num_dims; i++)
				{
					pts->m_data[i][r_index] = pts->m_data[i][temp_rand_loc];
				}
			}
		}
		// Set the number of pts->m_num_local_points to be equal to the number of random points
		pts->m_num_local_points = num_random_points;
		#ifdef _BDATS_DEBUG_
			printf ("MPI Rank: %d --> pts->m_num_local_points: %lld \n", mpi_rank, pts->m_num_local_points);
		#endif

		// Free the random_locations and random_loc_seen memory
		free (random_locations);
		free (random_loc_seen);

		H5Pclose(plist_id);
		MPI_Barrier(MPI_COMM_WORLD);
		t2 =  MPI_Wtime();

		#ifdef _BDATS_DEBUG_
		if(mpi_rank == 0 ){
			printf("Data read time %f \n", t2 - t1);
		}
		#endif

                // Print the first 100 elements in each dataset
		#if 0
		#ifdef DEBUG
			int64_t j;
			if (!cls)
			{
				for(i = 0 ; i < pts->m_num_dims; i++)
				{
					printf ("\nDataset [%d]: \n", i);
					for(j = 0; j < 100; j++)
					{
						printf ("%f \t", pts->m_data[i][j]);
					}
				}
			}
			printf ("\n");
			// Print the first 100 global indexes
			printf ("\nGlobal indexes of rank [%d]: ", mpi_rank);
			for (index = 0; index < 100; index++)
			{
				printf ("%lld \t", pts->m_vec_global_IDs[index]);
			}
			printf ("\n");
		#endif
		#endif
	
		MPI_Barrier(MPI_COMM_WORLD);
	}
		
	#if 1
        for(i = 0; i < pts->m_num_dims; i++)
        {
                current_dataset = &(data_info_array[i]);
                H5Sclose(current_dataset->open_space_id);
                H5Dclose(current_dataset->open_dataset_id);
        }
	#endif

        H5Fclose(file_id);
        //free (pts->m_data);
	//free (pts->m_class);
        //free (pts->m_vec_global_IDs);
}

// Each MPI process reads the top K percent of the points from the data block the process reads
void read_driver (char *filename, Points *pts, MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, dataset_info_t *data_info_array, float percent, bool cls)
{
        dataset_info_t *current_dataset;
        hid_t  file_id, plist_id;
        hsize_t         local_size;

        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, info);
        file_id = H5Fopen(filename, H5F_ACC_RDONLY,plist_id);

        double t1, t2;
        MPI_Barrier(MPI_COMM_WORLD);
        t1 =  MPI_Wtime();

        int i, index;
        int64_t last_process_extra_pts = 0;

        // Assuming that each dataset has the same size
        last_process_extra_pts = pts->m_num_global_points % mpi_size;
        if (mpi_rank ==  (mpi_size - 1))
        {
		// Read only the top given percentage of points
                pts->m_num_local_points = (pts->m_num_global_points/mpi_size + last_process_extra_pts) * (percent/100.0);
		#ifdef _BDATS_DEBUG_
		    printf ("MPI_Rank: %d --> m_num_local_points %lld \n", mpi_rank, pts->m_num_local_points); 
		#endif
        }
        else
        {
		// Read only the top given percentage of points
                pts->m_num_local_points = (pts->m_num_global_points/mpi_size) * (percent/100.0);
		#ifdef _BDATS_DEBUG_
		    printf ("MPI_Rank: %d --> m_num_local_points %lld \n", mpi_rank, pts->m_num_local_points); 
		#endif
        }

        // Allocate memory for Indexes
        pts->m_vec_global_IDs =  (BDATS_gIndex*) BDATS_malloc (pts->m_num_local_points * sizeof(BDATS_gIndex));
        if(pts->m_vec_global_IDs == NULL)
        {
                printf("Memory allocation for m_vec_global_IDs failed... \n");
                exit(-1);
        }

        int64_t global_offset = (pts->m_num_global_points/mpi_size) * mpi_rank; // pts->m_num_local_points*mpi_rank; 
	
        for (index = 0; index < pts->m_num_local_points; index++)
        {
                pts->m_vec_global_IDs[index] = global_offset + index;
        }

        // Allocate memory for data (i.e., pts->m_data and pts->m_class [if classification variable is true]).
        allocate_memory (pts, cls);

        // Go through each HDF5 dataset and read data into certain locations.
        for(i = 0; i < pts->m_num_dims; i++)
        {
                current_dataset = &(data_info_array[i]);
                current_dataset->open_file_id = file_id;

                //Open datasets and set hyperslab boundaries
                open_datasets (current_dataset, pts);

		// Suren: 10/02/2015: Changing float* buf to void *buf; typecasting to float *
		if (current_dataset->datatype == H5READER_DATA)
		{
			current_dataset->buf = *&pts->m_data[i];
			if(current_dataset->buf == NULL)
			{
				printf("Memory allocation for reading data fails ! \n");
				exit(-1);
			}
		}
		else 
		{
                	current_dataset->cls = pts->m_class;
			if(current_dataset->cls == NULL)
			{
				printf("Memory allocation for reading class data fails ! \n");
				exit(-1);
			}
		}

                //Parallel read
                read_segment(current_dataset);

		#ifdef _BDATS_DEBUG_
                if (mpi_rank == 0)
                        printf ("Finished reading dataset[%d] \n", i);
		#endif
	
		#if 0	
		#ifdef _BDATS_DEBUG_
		// Debug the contents of pts->m_class 
		
		if (current_dataset->datatype == H5READER_CLS_INT)
		{
                        for(int ki = 0; ki < 100; ki++)
                        {
                                printf ("%i \t", pts->m_class[ki]);
                        }
			printf ("\n");
		}
		#endif
		#endif
        }

        H5Pclose(plist_id);
        MPI_Barrier(MPI_COMM_WORLD);
        t2 =  MPI_Wtime();

	#ifdef _BDATS_DEBUG_
        if(mpi_rank == 0 ){
                printf("Data read time %f \n", t2 - t1);
        }
	#endif

                // Print the first 100 elements in each dataset
	#if 0
        #ifdef DEBUG
                int64_t j;
                for(i = 0 ; i < pts->m_num_dims; i++)
                {
                        printf ("\nDataset [%d]: \n", i);
                        for(j = 0; j < 100; j++)
                        {
                                printf ("%f \t", pts->m_data[i][j]);
                        }
                }
                printf ("\n");
                // Print the first 100 global indexes
                printf ("\nGlobal indexes of rank [%d]: ", mpi_rank);
                for (index = 0; index < 100; index++)
                {
                        printf ("%lld \t", pts->m_vec_global_IDs[index]);
                }
                printf ("\n");
        #endif
	#endif
	
        MPI_Barrier(MPI_COMM_WORLD);
	
	#if 1
	
        for(i = 0; i < pts->m_num_dims; i++)
        {
                current_dataset = &(data_info_array[i]);
                H5Sclose(current_dataset->open_space_id);
                H5Dclose(current_dataset->open_dataset_id);
        }
	#endif

        H5Fclose(file_id);
        //free (pts->m_data);
	//free (pts->m_class);
        //free (pts->m_vec_global_IDs);
}

// Assign a random percent of the points in the data read by each process to BDATS Points
int read_hdf5_file(char* filename, Points *pts, dataset_info_t*  data_info_array, int dataset_count, float percent, bool cls, bool random)
{
	int    mpi_size, mpi_rank;

     	MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Info info = MPI_INFO_NULL;

        MPI_Comm_size(comm, &mpi_size);
        MPI_Comm_rank(comm, &mpi_rank);

	if(dataset_count == 0) // read all columns
	        get_metadata_all_datasets (filename, data_info_array, pts, mpi_rank, mpi_size);
	else
	        get_metadata_select_datasets (filename, data_info_array, pts, mpi_rank, mpi_size, dataset_count);

        // Read data
        read_driver (filename, pts, comm, info, mpi_rank, mpi_size, data_info_array, percent, cls, random);

	return 0;
}

// MPI processes read the top k percent of the block of data allocated to each process
int read_hdf5_file(char* filename, Points *pts, dataset_info_t*  data_info_array, int dataset_count, float percent, bool cls)
{
	int    mpi_size, mpi_rank;

     	MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Info info = MPI_INFO_NULL;

        MPI_Comm_size(comm, &mpi_size);
        MPI_Comm_rank(comm, &mpi_rank);

	if(dataset_count == 0) // read all columns
	        get_metadata_all_datasets (filename, data_info_array, pts, mpi_rank, mpi_size);
	else
	        get_metadata_select_datasets (filename, data_info_array, pts, mpi_rank, mpi_size, dataset_count);

        // Read data
        read_driver (filename, pts, comm, info, mpi_rank, mpi_size, data_info_array, percent, cls);

	return 0;
}

// MPI processes read all the records of the hdf5 file
int read_hdf5_file(char* filename, Points *pts, dataset_info_t*  data_info_array, int dataset_count, bool cls)
{
	int    mpi_size, mpi_rank;

     	MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Info info = MPI_INFO_NULL;

        MPI_Comm_size(comm, &mpi_size);
        MPI_Comm_rank(comm, &mpi_rank);

	if(dataset_count == 0) // read all columns
	        get_metadata_all_datasets (filename, data_info_array, pts, mpi_rank, mpi_size);
	else
	        get_metadata_select_datasets (filename, data_info_array, pts, mpi_rank, mpi_size, dataset_count);

        // Read all data
	float percent = 100.0;
        read_driver (filename, pts, comm, info, mpi_rank, mpi_size, data_info_array, percent, cls);

	return 0;
}	 

H5PartFile* init_h5outfile (char* filename)
{
	H5PartFile* H5_out_fileid;
	H5_out_fileid = H5PartOpenFileParallel (filename, H5PART_WRITE | H5PART_FS_LUSTRE, MPI_COMM_WORLD);

	if (H5_out_fileid == NULL)
	{
		printf ("Error in opening file for init_h5outfile dataset \n");
		exit(-1);
	}
	return (H5_out_fileid);
}

H5PartFile* init_h5outfile_append (char* filename)
{
	H5PartFile* H5_out_fileid;
	//H5_out_fileid = H5PartOpenFileParallel (filename, H5PART_WRITE | H5PART_FS_LUSTRE, MPI_COMM_WORLD);
	H5_out_fileid = H5PartOpenFileParallel (filename, H5PART_APPEND | H5PART_FS_LUSTRE, MPI_COMM_WORLD);

	if (H5_out_fileid == NULL)
	{
		printf ("Error in opening file for init_h5outfile dataset \n");
		exit(-1);
	}
	return (H5_out_fileid);
}

int write_dataset_BDATS_lIndex (char* dataset_name, BDATS_lIndex* ids, int64_t num_pts, H5PartFile *H5_out_fileid)
{
	int    mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	int retval;

	double t1, t2;
        MPI_Barrier(MPI_COMM_WORLD);
        t1 =  MPI_Wtime();

	retval = H5PartSetStep(H5_out_fileid, 0);  // Assuming only one time step, with name /Step#0/
	if (retval != H5PART_SUCCESS)
	{
		printf ("Error in H5PartSetStep in init_h5outfile \n");
		exit(-1);
	}

	retval = H5PartSetNumParticles(H5_out_fileid, num_pts);
	if (retval != H5PART_SUCCESS)
	{
		printf ("[Rank: %d] Error in H5PartSetNumParticles in writing dataset: %s \n", mpi_rank, dataset_name);
	}

	retval = H5PartWriteDataInt32(H5_out_fileid, dataset_name, ids);
	if (retval != H5PART_SUCCESS)
	{
		printf ("[Rank: %d] Error in H5PartWriteDataInt32 in writing dataset: %s \n", mpi_rank, dataset_name);
	}

        MPI_Barrier(MPI_COMM_WORLD);
	t2 =  MPI_Wtime();

	#ifdef _BDATS_DEBUG_	
        if(mpi_rank == 0 ){
                printf("Time to write dataset %s: %f \n", dataset_name, t2 - t1);
        }
	#endif	
}

int write_dataset_BDATS_gIndex (char* dataset_name, BDATS_gIndex* ids, int64_t num_pts, H5PartFile *H5_out_fileid)
{
	int    mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	int retval;

	double t1, t2;
        MPI_Barrier(MPI_COMM_WORLD);
        t1 =  MPI_Wtime();

	retval = H5PartSetStep(H5_out_fileid, 0);  // Assuming only one time step, with name /Step#0/
	if (retval != H5PART_SUCCESS)
	{
		printf ("Error in H5PartSetStep in init_h5outfile \n");
		exit(-1);
	}

	retval = H5PartSetNumParticles(H5_out_fileid, num_pts);
	if (retval != H5PART_SUCCESS)
	{
		printf ("[Rank: %d] Error in H5PartSetNumParticles in writing dataset: %s \n", mpi_rank, dataset_name);
	}

	retval = H5PartWriteDataInt64(H5_out_fileid, dataset_name, ids);
	if (retval != H5PART_SUCCESS)
	{
		printf ("[Rank: %d] Error in H5PartWriteDataInt64 in writing dataset: %s \n", mpi_rank, dataset_name);
	}

        MPI_Barrier(MPI_COMM_WORLD);
	t2 =  MPI_Wtime();

	#ifdef _BDATS_DEBUG_	
        if(mpi_rank == 0 ){
                printf("Time to write dataset %s: %f \n", dataset_name, t2 - t1);
        }
	#endif	
}

int write_dataset_BDATS_decimal (char* dataset_name, BDATS_decimal* data, int64_t num_pts, H5PartFile *H5_out_fileid)
{
	int    mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	int retval;

	retval = H5PartSetStep(H5_out_fileid, 0);  // Assuming only one time step, with name /Step#0/
	if (retval != H5PART_SUCCESS)
	{
		printf ("[Rank: %d] Error in H5PartSetStep in writing dataset: %s \n", mpi_rank, dataset_name);
	}

	retval = H5PartSetNumParticles(H5_out_fileid, num_pts);
	if (retval != H5PART_SUCCESS)
	{
		printf ("[Rank: %d] Error in H5PartSetNumParticles in writing dataset: %s \n", mpi_rank, dataset_name);
	}

	double t1, t2;
        MPI_Barrier(MPI_COMM_WORLD);
        t1 =  MPI_Wtime();

	retval = H5PartWriteDataFloat32(H5_out_fileid, dataset_name, data);
	if (retval != H5PART_SUCCESS)
	{
		printf ("[Rank: %d] Error in H5PartWriteDataFloat32 in writing dataset: %s \n", mpi_rank, dataset_name);
	}

        MPI_Barrier(MPI_COMM_WORLD);
	t2 =  MPI_Wtime();
	#ifdef _BDATS_DEBUG_
        if(mpi_rank == 0 ){
                printf("Time to write dataset %s: %f \n", dataset_name, t2 - t1);
        }
	#endif
}

int close_h5outfile (H5PartFile *H5_out_fileid)
{
	int retval = H5PartCloseFile (H5_out_fileid);
	if (retval != H5PART_SUCCESS)
	{
		printf ("Error in H5PartCloseFile in close_h5outfile \n");
		exit (-1);
	}
}
