/*
* Enhansed Depth Estimation Reference Software
* ****************************************************************************
* Copyright (c) 2011/2013 Poznan University of Technology
*
* Address: 
*   Poznan University of Technology, 
*   Polanka 3 Street, Poznan, Poland
*
* Authors:
*   Krysztof Wegner     <kwegner@multimedia.edu.pl>
*   Olgierd Stankiewicz <ostank@multimedia.edu.pl>
*
* You may use this Software for any non-commercial purpose, subject to the
* restrictions in this license. Some purposes which can be non-commercial are
* teaching, academic research, and personal experimentation. You may also
* distribute this Software with books or other teaching materials, or publish
* the Software on websites, that are intended to teach the use of the 
* Software.
*
* Reference to the following source document:
*
* O. Stankiewicz, K. Wegner, M. Tanimoto, M. Doma?ki, 
* "Enhanced Depth Estimation Reference Software (DERS) for Free-viewpoint Television"
* ISO/IEC JTC1/SC29/WG11 MPEG2013/M31518 October 2013, Geneva, Switzerland
*
* are required in all documents that report any usage of the software.

* You may not use or distribute this Software or any derivative works in any
* form for commercial purposes. Examples of commercial purposes would be
* running business operations, licensing, leasing, or selling the Software, or
* distributing the Software for use with commercial products.
* ****************************************************************************
*/

/*
 * This software DERS(Depth Estimation Reference Software) was originally developed by
 * NAGOYA UNIVERSITY, JAPAN in the course of development of the ISO/IEC JTC1/SC29 WG 11 (MPEG) 3D Video
 * for reference purposes and its performance may not have been optimized.
 *
 * Those intending to use this software module in products are advised that its use may infringe
 * existing patents. ISO/IEC have no liability for use of this software module or modifications thereof.
 *
 * Assurance that the originally developed software module can be used
 *   (1) in the ISO/IEC JTC1/SC29 WG 11 (MPEG) 3D Video once the it is adopted to be used as reference
 *       software; and
 *   (2) to develop the codec for ISO/IEC JTC1/SC29 WG 11 (MPEG) 3D Video.
 *
 * To the extent that NAGOYA UNIVERSITY OR ANY OF ITS AFFILIATES owns patent rights that would be required to
 * make, use, or sell the originally developed software module or portions thereof included in the ISO/IEC
 * JTC1/SC29 WG 11 (MPEG) 3D Video in a conforming product, NAGOYA UNIVERSITY will assure the ISO/IEC that it
 * is willing to negotiate licenses under reasonable and non-discriminatory terms and conditions with
 * applicants throughout the world.
 *
 * NAGOYA UNIVERSITY retains full right to modify and use the code for its own purpose, assign or donate the
 * code to a third party and to inhibit third parties from using the code for products that do not conform
 * to MPEG-related and/or ISO/IEC International Standards.
 *
 * This copyright notice must be included in all copies or derivative works.
 * Copyright (c) ISO/IEC 2008.
 */
/*
 * In addition to the original author, NAGOYA UNIVERSITY Japan, this software was further modified by the
 * following parties (see ders_changes.txt and source for details):
 *
 *   TUT/Nokia Finland, NICT Japan, NTT Japan, GIST Korea, ETRI Korea, Poznan University Poland
 *
 *
 * The related parties retain full right to their code for their own purpose, assign or donate the cooresponding
 * code to another party and to inhibit third parties from using the code for products that do not conform
 * to MPEG-related and/or ISO/IEC International Standards.
 *
 */

#include "version.h"
#include "Estimation.h"
#include "yuv.h"
#include "ParameterDepthEstimation.h"
#include "Estimation.h"

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#ifdef WIN32
#include <crtdbg.h>
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifndef WIN32
#define BYTE unsigned char

#endif

#define CYCLE 2     //Frequency of Graph Cuts

int main(int argc, char *argv[])
{
    unsigned int n;
    CParameterDepthEstimation  cParameter;
    CEstimation cEstimation;
    int iCYCLE = CYCLE;

	//Char curDirectory[100];
	//LPSTR jj;
	//GetCurrentDirectory(100, curDirectory);
	//printf("%s", curDirectory);


#ifdef OUTPUT_COMPUTATIONAL_TIME
    clock_t start, finish, first;
    first = start = clock();
#endif

    printf("Depth Estimation Reference Software (DERS), Version %.1f\n", VERSION);
    printf("MPEG 3DV-FTV\n");
    printf("--------------------------------------------\n\n");

    if ( cParameter.xInit( argc, argv ) != 1 ) return 0;

    if(!cEstimation.xInit(cParameter)) return 10;


    int width = cParameter.getSourceWidth();
    int height = cParameter.getSourceHeight();

    // SHOW CONFIGURATIONS
    printf("Depth Type                   : %d\n", cParameter.getDepthType() );
    printf("Output Depth Map             : %s\n", cParameter.getFileOutputDepth().c_str() );
    printf("Image Resolution             : %dx%d\n", width, height );
    printf("Total Number of Frames       : %d\n", cParameter.getNumberOfFrames() );
    printf("Smoothing Coefficient        : %.2f\n", cParameter.getSmoothCoeff() );
    //Nagoya start
    printf("Smoothing Coefficient 2      : %.2f\n", cParameter.getSmoothCoeff2() );
    //Nagoya end

    CIYuv<DepthType> yuvDepth;
    CIYuv<ImageType> yuvLeft, yuvCenter, yuvRight;
    CIYuv<ImageType> yuvCenterSegment;// Nagoya
    CIYuv<ImageType> yuvCenter_prev;  // GIST
    CIYuv<ImageType> yuvReference;    // ETRI
    CIYuv<DepthType> yuvRefDepth;     // Nagoya

    CIYuv<ImageType> yuvErrors;		//Poznan Univ. - Owieczka

#ifdef SUB_PEL_PRECISION
    CIYuv<ImageType> yuvBuffer;
#ifdef SUB_PEL_VERTICAL_PRECISION
    CIYuv<ImageType> yuvBuffer2;
#endif
#endif

    if(!yuvDepth.resize(height, width, POZNAN_DEPTHMAP_CHROMA_FORMAT) ||
#ifdef SUB_PEL_PRECISION
#ifdef SUB_PEL_VERTICAL_PRECISION
            !yuvLeft.resize(height*cParameter.getVerticalPrecision(), width*cParameter.getPrecision(), 420) ||
            !yuvRight.resize(height*cParameter.getVerticalPrecision(), width*cParameter.getPrecision(), 420) ||
#else
            !yuvLeft.resize(height, width*cParameter.getPrecision(), 420) ||
            !yuvRight.resize(height, width*cParameter.getPrecision(), 420) ||
#endif
#else
            !yuvLeft.resize(height, width, 420) ||
            !yuvRight.resize(height, width, 420) ||
#endif
            !yuvCenter.resize(height, width, 420)
            //Nagoya start
            || !yuvCenterSegment.resize(height, width, 420)
            || !yuvRefDepth.resize(height, width, POZNAN_DEPTHMAP_CHROMA_FORMAT)
            //Nagoya end
            // GIST start
            || !yuvCenter_prev.resize(height, width, 420)
            // GIST end
#ifdef POZNAN_STORE_ERROR
            //Poznan Univ. - Owieczka start
            || !yuvErrors.resize(height,width,420)		//Owieczka
            //Poznan Univ. - Owieczka end
#endif
            )
    {
        fprintf(stderr, "Can't allocate enough memory.\n");
        return 3;
    }
#ifdef SUB_PEL_PRECISION
    switch(cParameter.getPrecision())
    {
    case 1:
        printf("Matching Precision           : INTEGER PIXEL\n");
        break;
    case 2:
        printf("Matching Precision           : HALF PIXEL\n");
        break;
    case 4:
        printf("Matching Precision           : QUATER PIXEL\n");
        break;
    default:
        fprintf(stderr, "Unknown value on Precision\n");
        break;
    }
#ifdef SUB_PEL_VERTICAL_PRECISION
    switch(cParameter.getVerticalPrecision())
    {
    case 1:
        printf("Matching Vertical Precision           : INTEGER PIXEL\n");
        break;
    case 2:
        printf("Matching Vertical Precision           : HALF PIXEL\n");
        break;
    case 4:
        printf("Matching Vertical Precision           : QUATER PIXEL\n");
        break;
    default:
        fprintf(stderr, "Unknown value on Vertical Precision\n");
        break;
    }
#endif

    if(cParameter.getPrecision()!=1
#ifdef SUB_PEL_VERTICAL_PRECISION
      ||cParameter.getVerticalPrecision()!=1
#endif
      )
    {
        if(!yuvBuffer.resize(height, width, 420)
#ifdef SUB_PEL_VERTICAL_PRECISION
          || !yuvBuffer2.resize(height*cParameter.getVerticalPrecision(), width, 420)
#endif
          )
        {
            fprintf(stderr, "Can't allocate enough memory\n");
            return 3;
        }
        if(!yuvLeft.setUpsampleFilter(cParameter.getFilter(), cParameter.getPrecision())) return 4;
        if(!yuvRight.setUpsampleFilter(cParameter.getFilter(), cParameter.getPrecision())) return 4;

        switch(cParameter.getFilter())
        {
        case 0:
            printf("Filter                       : (Bi-)linear Filter\n");
            break;
        case 1:
            printf("Filter                       : (Bi-)Cubic Filter\n");
            break;
        case 2:
            printf("Filter                       : MPEG-4 AVC Interpolation Filter\n");
            break;
        }

#ifdef SUB_PEL_VERTICAL_PRECISION
        //if(!yuvLeft.setVerticalUpsampleFilter(cParameter.getVerticalFilter(), cParameter.getVerticalPrecision())) return 4;
        //if(!yuvRight.setVerticalUpsampleFilter(cParameter.getVerticalFilter(), cParameter.getVerticalPrecision())) return 4;
        if(!yuvBuffer2.setVerticalUpsampleFilter(cParameter.getVerticalFilter(), cParameter.getVerticalPrecision())) return 4;

        switch(cParameter.getVerticalFilter())
        {
        case 0:
            printf("Vertical Filter              : (Bi-)linear Filter\n");
            break;
        case 1:
            printf("Vertical Filter              : (Bi-)Cubic Filter\n");
            break;
        case 2:
            printf("Vertical Filter              : MPEG-4 AVC Interpolation Filter\n");
            break;
        }
#endif

    }
#endif

    //Nagoya start
    if(!yuvCenterSegment.resize(height, width, 444)) return 2;
    //Nagoya end

    FILE *fpi_l, *fpi_c, *fpi_r;
    FILE *fpi_rd;
    FILE *fpo;

	FILE *fpo_bm;
#ifdef POZNAN_STORE_ERROR
    FILE *fpo_e;//Poznan Univ. - Owieczka - Errors File
#endif
    //ETRI start
    int **FirstFrame, **SecondFrame;
#ifdef SUB_PEL_VERTICAL_PRECISION
    int new_height = height * cParameter.getVerticalPrecision();
#else
    int new_height = height;
#endif
#ifdef SUB_PEL_PRECISION
    int new_width  = width * cParameter.getPrecision();
#else
    int new_width = width;
#endif
    if( cParameter.getDEmode() == 2) {
        // Allocate memory
        if (cParameter.getPrecision() != 1
#ifdef SUB_PEL_VERTICAL_PRECISION
          || cParameter.getVerticalPrecision() != 1
#endif          
          ) {
            yuvReference.setUpsampleFilter(cParameter.getFilter(), cParameter.getPrecision());
#ifdef SUB_PEL_VERTICAL_PRECISION
            yuvReference.setVerticalUpsampleFilter(cParameter.getVerticalFilter(), cParameter.getVerticalPrecision());
            yuvReference.resize(height*cParameter.getVerticalPrecision(), width*cParameter.getPrecision(), 420);
#else
            yuvReference.resize(height, width*cParameter.getPrecision(), 420);
#endif
        } else {
            yuvReference.resize(height, width, 420);
        }
        FirstFrame = new int *[new_height];
        SecondFrame = new int *[new_height];
        for (int i=0; i<new_height; i++)
        {
            FirstFrame[i] = new int [new_width];
            SecondFrame[i] = new int [new_width];
        }
        for (int j=0; j<new_height; j++)
            for (int i=0; i<new_width; i++)
            {
                FirstFrame[j][i] = 0;
                SecondFrame[j][i] = 0;
            }
    }
    // ETRI end

	if (cParameter.getFileLeftView().length() > 0) {
		if ((fpi_l = fopen(cParameter.getFileLeftView().c_str(), "rb")) == NULL)
		{
			fprintf(stderr, "Can't open left YUV image file.[%s]\n", cParameter.getFileLeftView().c_str());
			return 5;
		}
	}
	else {
		fpi_l = NULL;
		fprintf(stderr, "No Left image defined\n");
	}

	if (cParameter.getFileCenterView().length() > 0) {
		if ((fpi_c = fopen(cParameter.getFileCenterView().c_str(), "rb")) == NULL)
		{
			fprintf(stderr, "Can't open center YUV image file.[%s]\n", cParameter.getFileCenterView().c_str());
			return 5;
		}
	}
	else {
		fpi_c = NULL;
		fprintf(stderr, "No Center image defined\n");
		return 5;
	}
	if (cParameter.getFileRightView().length() > 0) {
		if ((fpi_r = fopen(cParameter.getFileRightView().c_str(), "rb")) == NULL)
		{
			fprintf(stderr, "Can't open right YUV image file.[%s]\n", cParameter.getFileRightView().c_str());
			return 5;
		}
	}
	else {
		fprintf(stderr, "No Right image defined\n");
		fpi_r = NULL;
		if (fpi_l == NULL) // 중간 말고는 정의가 하나도 안되어 있으면 계산 불가능함. 
			return 5;
	}
    if((fpo=fopen(cParameter.getFileOutputDepth().c_str(), "wb"))==NULL)
    {
        fprintf(stderr, "Can't open output file.[%s]\n", cParameter.getFileOutputDepth().c_str());
        return 6;
    }
    if(cParameter.getDEmode() == 3 && (fpi_rd=fopen(cParameter.getRefDepthFile().c_str(), "rb"))==NULL)
    {
        fprintf(stderr, "Can't open RefDepth YUV file.[%s]\n", cParameter.getRefDepthFile().c_str());
        return 5;
    }

	if ((fpo_bm = fopen("depth_bm.yuv", "wb")) == NULL)  // 중간결과 확인을 위하여 
	{
		fprintf(stderr, "Can't open debug output file.[%s]\n", "depth_bm.yuv");
		return 7;
	}

#ifdef POZNAN_STORE_ERROR
    //Poznan Univ. - Owieczka
    if((fpo_e=fopen(cParameter.getFileOutputErrors().c_str(), "wb"))==NULL)
    {
        fprintf(stderr, "Can't open output errors file.[%s]\n", cParameter.getFileOutputErrors().c_str());
        return 6;
    }
#endif

    /* ---------------------- z values ---------------------- */
    char zv_file_name[] = "z_value.txt"; //z value output file
    FILE *zf;

    if((zf=fopen("z_value.txt", "at"))==NULL)
    {
        printf("Can't open file for z values\n");
        printf("Z values are output only on stdout\n");
    }
    else
    {
        fprintf(zf, "%s(Depth Type %d)\n", cParameter.getFileOutputDepth().c_str(), cParameter.getDepthType());
        if(cParameter.getDepthType())
            fprintf(zf, "Depth map from the origin of 3D space is outputted. \n"); //Indicate z values from the origin of 3D space
        else
            fprintf(zf, "Depth map from camera is outputted. \n"); //Indicate z values from camera
        fprintf(zf, "NearestDepthValue    %.6lf\n", cEstimation.getZnear());
        fprintf(zf, "FarthestDepthValue   %.6lf\n", cEstimation.getZfar());
        fprintf(zf, "Input these values to View Synthesis software.\n\n");
        fclose(zf);
    }
    if(cParameter.getDepthType())
        printf("Depth map from the origin of 3D space is outputted. \n"); //Indicate z values from the origin of 3D space
    else
        printf("Depth map from camera is outputted. \n"); //Indicate z values from camera
    printf("NearestDepthValue    %.6lf\n", cEstimation.getZnear());
    printf("FarthestDepthValue   %.6lf\n", cEstimation.getZfar());
    printf("Input these values to View Synthesis software.\n");
    printf("--------------------------------------------\n\n");

#ifdef OUTPUT_COMPUTATIONAL_TIME
    finish = clock();
    printf( "Initialization: %.4f sec\n", (double)(finish - start) / CLOCKS_PER_SEC);
    start = finish;
#endif
    // GIST start
    bool isfirstframe = true;
    // GIST end

// ==========================================================================================================
    for(n=cParameter.getStartFrame(); n<cParameter.getStartFrame()+cParameter.getNumberOfFrames(); n++) // DT
    {
#ifdef OUTPUT_COMPUTATIONAL_TIME
		//printf("Segmentation: %.4f sec\n", ((double)(clock() - finish)) / ((double)CLOCKS_PER_SEC));
		finish = clock();
#endif

        // ------------- Read YUV 4:2:0 Image -------------
#ifdef SUB_PEL_PRECISION
        if(cParameter.getPrecision()==1
#ifdef SUB_PEL_VERTICAL_PRECISION
          && cParameter.getVerticalPrecision()==1
#endif
          )
        {
            if( (fpi_l != NULL && !yuvLeft.readOneFrame(fpi_l, n)) ||
                (fpi_r != NULL && !yuvRight.readOneFrame(fpi_r, n)) ||
                !yuvCenter.readOneFrame(fpi_c, n)   ) break;
        }
        else
        {
			if (fpi_l != NULL && yuvBuffer.readOneFrame(fpi_l, n)) { // break;
#ifdef SUB_PEL_VERTICAL_PRECISION
				yuvBuffer2.upsamplingVertical(&yuvBuffer);
				yuvLeft.upsampling(&yuvBuffer2);
#else
				yuvLeft.upsampling(&yuvBuffer);
#endif
			} 
			if (fpi_r != NULL && yuvBuffer.readOneFrame(fpi_r, n)) { //break;
#ifdef SUB_PEL_VERTICAL_PRECISION
				yuvBuffer2.upsamplingVertical(&yuvBuffer);
				yuvRight.upsampling(&yuvBuffer2);
#else
				yuvRight.upsampling(&yuvBuffer);
#endif
			}
            if(!yuvCenter.readOneFrame(fpi_c, n)) break;
        }
#else
        if( !yuvLeft.readOneFrame(fpi_l, n) ||
            !yuvRight.readOneFrame(fpi_r, n) ||
            !yuvCenter.readOneFrame(fpi_c,n )) break;
#endif

        //Nagoya start
        if( cParameter.getDEmode() == 3 && !yuvRefDepth.readOneFrame(fpi_rd, n) ) {
            printf("EOF RefDepth\n");
            break;
        }

        if( cParameter.getDEmode() == 1 || cParameter.getDEmode() == 2 ) {
            cEstimation.load_man_images(cParameter.getFileCenterManual().c_str(), n, isfirstframe);
        }
#ifdef OUTPUT_COMPUTATIONAL_TIME
		printf("Prepare: %.4f sec\n", ((double)(clock() - finish)) / ((double)CLOCKS_PER_SEC));
		finish = clock();
#endif
        // ----------- Image Segmentation -----------

		printf("cEstimation.getImageSegmentation()=%d\n", cEstimation.getImageSegmentation());

        if(cEstimation.getImageSegmentation()==1 && !cEstimation.refresh_frame()) {
            yuvCenterSegment.setData444_inIBGR(&yuvCenter);  // YUV420 to BGR444 
            cEstimation.center_image_segmentation(yuvCenterSegment.getData());
        }
#ifdef OUTPUT_COMPUTATIONAL_TIME
		printf("Segmentation: %.4f sec\n", ((double)(clock() - finish)) / ((double)CLOCKS_PER_SEC));
		finish = clock();
#endif
        //Nagoya end

        // ------------- Block Matching -------------
        printf("Block Matching\n");
        // GIST start
        bool temporalonoff = false;
        if(cParameter.getTemporalEnhancement())
            isfirstframe?temporalonoff=false:temporalonoff=true;
        // GIST end
        if(cParameter.getDEmode() == 2) temporalonoff=false;
        if(cParameter.getDEmode() == 1) temporalonoff=true;

#ifdef POZNAN_TWOVIEW_SUPPORT
		printf("Matching Driection: %d\n", cParameter.getDirections());
		if(cParameter.getDirections() == 1)      cEstimation.block_matching(&yuvLeft, NULL, &yuvCenter, yuvCenterSegment.getData(), &yuvCenter_prev, temporalonoff, cParameter.getThreshold() );
        else if(cParameter.getDirections() == 2) cEstimation.block_matching(NULL, &yuvRight, &yuvCenter, yuvCenterSegment.getData(), &yuvCenter_prev, temporalonoff, cParameter.getThreshold() );
        else                                     cEstimation.block_matching(&yuvLeft, &yuvRight, &yuvCenter, yuvCenterSegment.getData(), &yuvCenter_prev, temporalonoff, cParameter.getThreshold() );
#else
        cEstimation.block_matching(&yuvLeft, &yuvRight, &yuvCenter, yuvCenterSegment.getData(), &yuvCenter_prev, temporalonoff, cParameter.getThreshold() );
#endif
      
#ifdef POZNAN_STORE_ERROR
        cEstimation.store_errors(&yuvErrors,fpo_e);//Poznan Univ. - Owieczka
#endif

#ifdef OUTPUT_COMPUTATIONAL_TIME
		printf("BlockMatching: %.4f sec\n", ((double)(clock() - finish)) / ((double)CLOCKS_PER_SEC));
		finish = clock();
#endif
		// depth는 아직 결정 안됨. cEstimation.error에 BM 결과가 저장되고 GC을 수행해야만 결과가 나옴.
		// 따라서 NULL GC를 만들어서 BM의 결과를 확인하게 해야함. 
		//yuvCenter.writeOneFrame(fpo_bm);


        // ------------- Update error cost -------------
        iCYCLE = cEstimation.update_error_cost(&cParameter, &yuvRefDepth, &yuvCenter, isfirstframe, CYCLE, temporalonoff);

        // ------ Depth Estimation by Graph-cut -----
		// OpenCV도 BM이후에 GC을 적용할 수 있을까?

		printf("cEstimation.getImageSegmentation()=%d\n", cEstimation.getImageSegmentation());

        //printf("Graph Cuts\n");
        if(cParameter.getDEmode() != 0) {
            cEstimation.depth_estimation_by_graph_cut_semi(yuvDepth.Y, iCYCLE, &yuvCenter);
        }
#ifdef SEOULTECH_CUDA_SUPPORT 
        else if(cParameter.getCudaCheck() == 1){
            printf("Graph Cut with CUDA: %d\n", cEstimation.getImageSegmentation());
            cEstimation.depth_estimation_by_graph_cut_cuda(yuvDepth.Y, iCYCLE, yuvCenterSegment.getData(), &yuvCenter, cParameter.getCudaDataCoeff(), cParameter.getCudaSmoothCoeff(), cParameter.getCudaStochatic());
        } else if(cParameter.getGraphcutNoAuxCheck() == 1){
            printf("Graph Cut without Auxility node: %d\n", cEstimation.getImageSegmentation());
            cEstimation.depth_estimation_by_graph_cut_no_auxnode(yuvDepth.Y, iCYCLE, yuvCenterSegment.getData(), &yuvCenter);
        }
#endif
        else if(cEstimation.getImageSegmentation()==1) {
#ifdef POZNAN_OCC
          if(cParameter.getOcc()==1) {
            printf("Graph Cuts Seg Occ\n");
            cEstimation.depth_estimation_by_graph_cut_segmentation_occ(yuvDepth.Y, CYCLE, yuvCenterSegment.getData(), &yuvCenter);
          }
          else 
#endif
          {
            printf("Graph Cuts\n");
            cEstimation.depth_estimation_by_graph_cut_segmentation(yuvDepth.Y, CYCLE, yuvCenterSegment.getData(), &yuvCenter);
          }
        } 
#ifdef POZNAN_OCC
        else if(cParameter.getOcc()==1) {
          printf("Graph Cuts Occ\n");
           cEstimation.depth_estimation_by_graph_cut_occ(yuvDepth.Y, iCYCLE, yuvCenterSegment.getData(), &yuvCenter);
        } 
#endif
        else {
            printf("Graph Cuts: %d\n", cEstimation.getImageSegmentation());
		
			// 이걸 생략하면 하나도 값이 안나옴. 왜지? 
            cEstimation.depth_estimation_by_graph_cut(yuvDepth.Y, iCYCLE, yuvCenterSegment.getData(), &yuvCenter);
        }

#ifdef OUTPUT_COMPUTATIONAL_TIME
		printf("Graph Cuts: %.4f sec\n", ((double)(clock() - finish)) / ((double)CLOCKS_PER_SEC));
		finish = clock();
#endif
		yuvDepth.writeOneFrame(fpo_bm);  //  GC 이후의 결과를 보기 위하여 

        //Nagoya start
        // -------------- Plane Fitting -------------
        if(cEstimation.getImageSegmentation()==1 && !cEstimation.refresh_frame())
        {
            printf("Plane Fitting\n");
            cEstimation.plane_fitting(yuvDepth.Y, yuvCenterSegment.getData(), cEstimation.getNumOfSegms()); //Added by SZK
        }
        //Nagoya end
#ifdef OUTPUT_COMPUTATIONAL_TIME
		printf("Plane Fitting : %.4f sec\n", ((double)(clock() - finish)) / ((double)CLOCKS_PER_SEC));
		finish = clock();
#endif
        // -------------- Depth Post processing -------------
	

        if (cParameter.getDEmode() != 2) {
            cEstimation.depth_estimation_post_processing(yuvDepth.Y, &cParameter);
        } else {
            // ETRI start
            // ------------- Read YUV 4:2:0 Image (previous) -------------
            if (!isfirstframe) {
                if (cParameter.getPrecision() == 1
#ifdef SUB_PEL_VERTICAL_PRECISION
                  && cParameter.getVerticalPrecision() == 1
#endif
                  )
                    yuvReference.readOneFrame(fpi_c, n-1);
                else
                {
                    yuvBuffer.readOneFrame(fpi_c, n-1);
#ifdef SUB_PEL_VERTICAL_PRECISION
                    yuvBuffer2.upsamplingVertical(&yuvBuffer);
                    yuvReference.upsampling(&yuvBuffer2);
#else
                    yuvReference.upsampling(&yuvBuffer);
#endif
                }

                for (int j=0; j<new_height; j++)
                    for (int i=0; i<new_width; i++)
                        FirstFrame[j][i] = yuvReference.Y[j][i];

                // ------------- Read YUV 4:2:0 Image (current) -------------
                if (cParameter.getPrecision() == 1
#ifdef SUB_PEL_VERTICAL_PRECISION
                  && cParameter.getVerticalPrecision() == 1
#endif
                  )
                    yuvReference.readOneFrame(fpi_c, n);
                else
                {
                    yuvBuffer.readOneFrame(fpi_c, n);
#ifdef SUB_PEL_VERTICAL_PRECISION
                    yuvBuffer2.upsamplingVertical(&yuvBuffer);
                    yuvReference.upsampling(&yuvBuffer2);
#else
                    yuvReference.upsampling(&yuvBuffer);
#endif
                }

				// 왜 memcpy 안하니?
                for (int j=0; j<new_height; j++)
                    for (int i=0; i<new_width; i++)
                        SecondFrame[j][i] = yuvReference.Y[j][i];
            }
            cEstimation.depth_estimation_post_processing_etri(yuvDepth.Y, &cParameter, FirstFrame, SecondFrame);
            // ETRI end

        }
#ifdef OUTPUT_COMPUTATIONAL_TIME
		printf("Post processing : %.4f sec\n", ((double)(clock() - finish)) / ((double)CLOCKS_PER_SEC));
		finish = clock();
#endif

        // ------------- Write Depth Map ------------
        yuvDepth.writeOneFrame(fpo);


#ifdef OUTPUT_COMPUTATIONAL_TIME
        finish = clock();
        printf( "Frame%03d: %.4f sec\n", n, (double)(finish - start) / CLOCKS_PER_SEC);
        start = finish;
#endif
        // GIST start
        isfirstframe = false;

        for(int y=0; y<height; y++)
        {
            for(int x=0; x<width; x++)
            {
                yuvCenter_prev.Y[y][x] = yuvCenter.Y[y][x];
            }
        }
        // GIST end

        printf("\n");
    } // frame loop
// ==========================================================================================================

    //ETRI start
    if( cParameter.getDEmode() == 2) {
        // Free memory
        for (int i=0; i<new_height; i++)
        {
            delete[] FirstFrame[i];
            delete[] SecondFrame[i];
        }
        delete[] FirstFrame;
        delete[] SecondFrame;
    }
    //ETRI end

	if(fpi_l) fclose(fpi_l);
	if(fpi_c) fclose(fpi_c);
    if(fpi_r) fclose(fpi_r);
    if(fpo)   fclose(fpo);
	if (fpo_bm)   fclose(fpo_bm);
    if(cParameter.getDEmode() == 3) fclose(fpi_rd);
#ifdef POZNAN_STORE_ERROR
    fclose(fpo_e);
#endif
/*
    if(cParameter.getDepthType())
        printf("Depth map from the origin of 3D space is outputted. \n"); //Indicate z values from the origin of 3D space
    else
        printf("Depth map from camera is outputted. \n"); //Indicate z values from camera
    printf("NearestDepthValue    %lf\n", cEstimation.getZnear());
    printf("FarthestDepthValue   %lf\n", cEstimation.getZfar());
    printf("Input these values to View Synthesis software.\n\n");
*/

#ifdef OUTPUT_COMPUTATIONAL_TIME
    finish = clock();
    printf("Total: %.4f sec\n", ((double)(finish-first))/((double)CLOCKS_PER_SEC));
#endif

#if defined(WIN32) && defined(_DEBUG)
    _CrtDumpMemoryLeaks();
#endif

    return 0;
}
