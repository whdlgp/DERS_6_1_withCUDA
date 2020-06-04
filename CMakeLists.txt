cmake_minimum_required(VERSION 2.8)
project( DERS-6.1 LANGUAGES CXX CUDA)
cmake_policy(SET CMP0015 NEW)

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set( COMMON_LIB_STATIC_SRC )
aux_source_directory( "./CommonLibStatic/src" COMMON_LIB_STATIC_SRC)
set( COMMON_LIB_STATIC_INCLUDE_DIR ./CommonLibStatic/include )

set( DEPTH_EST_LIB_STATIC_SRC )
aux_source_directory( "./DepthEstLibStatic/src" DEPTH_EST_LIB_STATIC_SRC)
set( DEPTH_EST_LIB_STATIC_INCLUDE ./DepthEstLibStatic/include )

set( DEPTH_EST_SRC )
aux_source_directory( "./DepthEst/src" DEPTH_EST_SRC)
set( DEPTH_EST_INCLUDE ./DepthEst/include )

set( CUDA_CUTS_SRC )
aux_source_directory( "./cuda_cuts" CUDA_CUTS_SRC)
set( CUDA_CUTS_INCLUDE ./cuda_cuts )

set( ld_flag opencv-2.4.13.6/install/lib )
set( include_flag opencv-2.4.13.6/install/include )
set( libs 
    opencv_calib3d
    opencv_core
    opencv_features2d
    opencv_flann
    opencv_gpu
    opencv_highgui
    opencv_imgproc
    opencv_legacy
    opencv_ml
    opencv_nonfree
    opencv_objdetect
    opencv_ocl
    opencv_photo
    opencv_stitching
    opencv_superres
    opencv_video
    opencv_videostab
    OpenMP::OpenMP_CXX )

find_package(OpenMP REQUIRED)

include_directories( ${include_flag} ${COMMON_LIB_STATIC_INCLUDE_DIR} ${DEPTH_EST_LIB_STATIC_INCLUDE} ${DEPTH_EST_INCLUDE} ${CUDA_CUTS_INCLUDE} )
link_directories( ${ld_flag} )

add_library( CommonLibStatic STATIC ${COMMON_LIB_STATIC_SRC} )
add_library( DepthEstLibStatic STATIC ${DEPTH_EST_LIB_STATIC_SRC} ${CUDA_CUTS_SRC} )

target_link_libraries(CommonLibStatic ${libs} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries(DepthEstLibStatic ${libs} ${CMAKE_THREAD_LIBS_INIT} )

add_executable( DepthEst ${DEPTH_EST_SRC} )
target_link_libraries(DepthEst ${libs} ${CMAKE_THREAD_LIBS_INIT} DepthEstLibStatic CommonLibStatic )
