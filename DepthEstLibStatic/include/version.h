#ifndef __INCLUDE_VERSION_H__
#define __INCLUDE_VERSION_H__

#define VERSION 6.1

#define OUTPUT_COMPUTATIONAL_TIME

#define SUB_PEL_PRECISION
#define SUB_PEL_VERTICAL_PRECISION    //Poznan University of Technology - Owieczka - Increase matching precision in vertical direction
#define ACCURATE_CORRESPONDENCE

#define POZNAN_DYNAMIC_TABLE            //Poznan University of Technology - Owieczka - Dynamic Table insted of Static Array
#define POZNAN_ERRORS_DISTRIBUTED_ALLOC //Poznan University of Technology - Owieczka 

//#define POZNAN_TWO_COSTS

//#define POZNAN_STORE_ERROR            //Poznan University of Technolgy - Owieczka - Store Matching Volument in the file for debug

//#define POZNAN_16BIT_DEPTH

//#define POZNAN_DOUBLE_COST

//#define POZNAN_ZNEAR_ZFAR_SEARCH_RANGE //Poznan University of Technology - Owieczka - Search Range from zNear and zFar

//#define POZNAN_OCC                     //Poznan University of Technology - Owieczka - Use advence occlusion handling
//#define POZNAN_OCC_ALWAYS_DEPTH        //Poznan University of Technology - Owieczka - Use always depth insted of disparity for occlusion
//#define POZNAN_OCC_VERBOSE

#define POZNAN_TWOVIEW_SUPPORT

#define SEOULTECH_CUDA_SUPPORT

#define POZNAN_DEPTHMAP_CHROMA_FORMAT 420 //420

#ifdef POZNAN_DOUBLE_COST
typedef double CostType;
#define COST_MAX 255
#define COST_INF SHRT_MAX
#else
typedef int CostType;
#define COST_MAX 255
#define COST_INF SHRT_MAX
#endif

#ifdef POZNAN_16BIT_DEPTH
typedef unsigned short DepthType;
typedef unsigned char ImageType;
#define MAX_DEPTH (256*256)
#define MAX_LUMA 256
#else
typedef unsigned char DepthType;
typedef unsigned char ImageType;
#define MAX_DEPTH 256
#define MAX_LUMA 256
#endif

#endif


