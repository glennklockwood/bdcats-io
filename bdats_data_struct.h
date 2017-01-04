/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */
/*										*/
/*   Files:         Have to list the files  					*/
/*										*/	
/*   Description:  BDATS: BigData Analytics at Trillion Scale    		*/ 
/*                                                                           	*/
/*   Author:  Md. Mostofa Ali Patwary                                        	*/
/*            Research Scientist, Parallel Computing Lab, Intel Corporation    	*/
/*            email: mostofa.ali.patwary@intel.com                          	*/
/*										*/
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  */

#ifndef _BDATS_DATA_STRUCT_
#define _BDATS_DATA_STRUCT_

#include "bdats_headers.h"

//#define _BDATS_DEBUG_
#define STAT_COMM_TIME

//#define _ADD_RANDOMIZATION_IN_COORDINATES_
#define MEMORY_CRITICAL_DATASET

// general
#define MASTER_NODE 0
#define NOT_ALIGNED_MALLOC
#define PARALLEL_FOR_CHUNK_SIZE 64
#define MAX_MPI_CHUNK_SIZE 134212727 //268425455
#define SPECIFIC_MAX_MPI_CHUNK 134212727 //262144
#define MAX_EXPECTED_CHUNK 20 //6

// partitioning and gathering
#define TOTAL_SAMPLES_PER_BILLION_FOR_MEDIAN 131072 //4096 // change to smaller value later
#define SAMPLES_PER_NODE 256

#define RAND_SEED 1234567
#define PARTITIONING_COMMU_BUF 6 //6 WORK FINE with large dataset //4 //6 //3
#define GATHERING_COMMU_BUF 3 //3 WORK FINE with large dataset//2 //3
#define THREAD_LOCAL_BUFFER 1024
#define MPI_COMM_BEGIN_TAG 100 

// kdtree defines
#define KDT_LEAF_BUCKET_SIZE 31 //255
#define KDT_MAX_NODE_COUNT_MUL_FACTOR 3 //3 WORK FINE with large dataset//2 //4 //16 // if this is close to the above line, BDATS_kdt_cIndex should be int64_t 
#define KDT_FIRST_PART_LEAF_NODES_PER_THREAD 3 // increase for better load balancing 
#define KDT_DEPTH_MUL_FACTOR 4 //16
#define KDT_MAX_SAMPLES 8192 //512 //2048 //1024 //524288 //131072 //1024
#define KDT_MIN_SAMPLES 256
#define KDT_SAMPLES_PERCENT 0.5 // may not be used everywhere
#define KDT_SORT_SIZE 16384

// dbscan defines
#define MAX_EXACT_NEIGHBORS 65536 //8388608 //16384 //8388608 //1048576 //16384 //9216 //8192
#define NO_OWNER 255
#define THREAD_MEM_MUL 4 //4 WORK FINE with large dataset //8 //5 //9 (worked well 16)
#define NODE_MEM_MUL 4 //4 WORK FINE with large dataset //8 //3 //6 // increase sime more to 12
#define CHECK_UNLIMITED -1
#define THREAD_NODE_BUFFER 256
#define DBSCAN_PER_THREAD_QUEUE_MUL 2
#define NO_OWNER 255
#define DBSCAN_LOCAL_COMP_CHUNK 4096
#define LOCK_REDUCTION_FACTOR 40 //10

// parallel sort defines
#define MIN_NUM_LOCAL_PIVOTS 32 // 32 128 // 32
#define MAX_NUM_LOCAL_PIVOTS 256 //256 2048 //32
//#define APP_ASTROPHYSICS
//#define _BDATS_DEBUG_DETAIL
#define PARALLEL_SORT_OBJECT_NOT_DONE
#define DBSCAN_SPECIFIC_SORT

// used only in debug time
#define MAX_NEIGHBOR_EXPECTATION 300000 // 300000 WORK FINE with large dataset

// find cluster IDs
#define NOISE_ID -1
#define MIN_CLUS_SIZE 2 //2 //2 //1000000 //100000 //2 //
#define WRITE_OR_APPEND 1 //20 // in billion, small than this will append

// memory allocation
#ifdef NOT_ALIGNED_MALLOC
#define BDATS_malloc(x) malloc(x)
#define BDATS_free(x) free(x) 
#define BDATS_realloc(x,y) realloc(x,y)
#else
#define BDATS_malloc(x) _mm_malloc(x, 64) 
#define BDATS_free(x) _mm_free(x) 
#endif

#define POW2(x) (1 << (x))
#define LOWER(i) (i<<1)
#define UPPER(i) ((i<<1)+1)


typedef float 	BDATS_decimal;
#define BDATS_decimal_mpi MPI_FLOAT

typedef int64_t BDATS_gIndex; 
#define BDATS_gIndex_mpi MPI_INT64_T

typedef int BDATS_lIndex;
#define BDATS_lIndex_mpi MPI_INT

typedef int	BDATS_kdt_cIndex;
typedef float   BDATS_kdt_decimal;

struct BoundingBox
{
	BDATS_decimal 	m_lower;
	BDATS_decimal	m_upper;
};

struct clus_solution
{
	BDATS_gIndex m_cluster_IDs;
	BDATS_lIndex m_nei_count;	
};

struct Points
{
	int64_t 	m_num_local_points; 	// local points count including gathered
	int64_t 	m_num_dims; 		// dimension of each point
	BDATS_decimal**	m_data;			// coordinates: 2D array, each D for one dimension data together
	int* 		m_class; 
	///BoundingBox* 	m_local_bbox; // local bounding box 
	BoundingBox*    m_global_bbox; // global bounding box 
	
	int*		m_vec_prIDs;		// 32 bit should be good enough 
	BDATS_lIndex*	m_vec_local_IDs;	// could be 32 bit or 64 bit, assuming 64 bit
	BDATS_gIndex*	m_vec_global_IDs;	// 64 bit, no doubt  	
	BDATS_gIndex*	m_vec_cluster_IDs; 

	int64_t		m_num_global_points; 	
	int		m_mpi_comm_current_tag;
	int64_t 	m_num_points_gathered; // points gathered from neighbors
};

#endif
