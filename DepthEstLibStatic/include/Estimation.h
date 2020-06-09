#include "graph.h"
#include "version.h"

#ifndef WIN32
#define BYTE unsigned char
#endif

#ifndef _INCLUDE_ESTIMATION_H_
#define _INCLUDE_ESTIMATION_H_

#pragma warning(disable:4996)
#pragma warning(disable:4244)
#pragma warning(disable:4819)

#ifndef _INCLUDE_OPEN_CV_
#define _INCLUDE_OPEN_CV_
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

//#pragma comment(lib, "cv.lib")
//#pragma comment(lib, "cxcore.lib")
//#pragma comment(lib, "highgui.lib")
//#pragma comment(lib, "cvaux.lib")

#endif // _INCLUDE_OPEN_CV_

#include "yuv.h"

class CParameterDepthEstimation;

enum
{
    GC_AUTO=0,        // Automatic depth estimation
    SEMI_D_REFRESH,   // Semi-automatic, refresh frame based on refresh Depth map
    SEMI_D_E_REFRESH, // Semi-automatic, refresh frame based on manual Disparity and Manual Edge map
    SEMI_TEMPORAL,    // Semi-automatic, temporal consistency
};


class CEstimation
{
public:
    CEstimation                     ();
    virtual ~CEstimation            ();

    bool                xInit                   (CParameterDepthEstimation cPrameter);
    void                free_memory             ();

    double              getZnear                    ()      { return m_dZnear; }
    double              getZfar                     ()      { return m_dZfar; }
    int                 getImageSegmentation        ()      { return m_iImageSegmentation; }
    int                 getNumOfSegms               ()      { return m_iNumOfSegms; }

    void                block_matching              (CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter_prev, bool temporalonoff, double threshold);
    void                depth_estimation_by_graph_cut               (DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter);
#ifdef POZNAN_OCC
    void                depth_estimation_by_graph_cut_occ           (DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter);
    void                depth_estimation_by_graph_cut_segmentation_occ  (DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter);
#endif
    void                depth_estimation_by_graph_cut_segmentation  (DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter);
    void                depth_estimation_by_graph_cut_semi          (DepthType **pDepth, int iCycle, CIYuv<ImageType> *yuvCenter);
#ifdef SEOULTECH_CUDA_SUPPORT
    void                depth_estimation_by_graph_cut_cuda          (DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter, double datacoeff, double smoothcoeff, int is_stochatic);
    void                depth_estimation_by_graph_cut_no_auxnode    (DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter);
#endif
    void                center_image_segmentation   (BYTE ***srcSEGM);

    void                plane_fitting               (DepthType **pDepth, BYTE ***srcSEGM, int NumSegm);

    bool                load_man_images             (const char* m_cFileManual, int framenumber, bool isfirstframe);
    bool                refresh_frame               (){return bool(gc_mode == SEMI_D_REFRESH || gc_mode == SEMI_D_E_REFRESH);}
    int                 update_error_cost           (CParameterDepthEstimation *cParameter, CIYuv<DepthType> *yuvRefDepth, CIYuv<ImageType> *yuvCenter, bool isfirstframe, int GC_cycle, bool temporalonoff);
    void                depth_estimation_post_processing      (DepthType **pDepth, CParameterDepthEstimation *cParameter);
    void                depth_estimation_post_processing_etri (DepthType **pDepth, CParameterDepthEstimation *cParameter, int **FirstFrame, int **SecondFrame);
    //Nagoya end

    //Owieczka start
    void                store_errors                (CIYuv<ImageType> *yuvErrors,FILE *fpo_e);
#ifdef POZNAN_OCC
    double              (CEstimation::*func_data_term)   (int label, int pp);
    inline double       data_term                   (int label, int pp);
  

    double              data_term_disparity         (int label, int pp);
    double              data_term_depth             (int label, int pp);

    void                (CEstimation::*func_estimation_occ)   ();
    inline void         estimation_occ              ();

    void                estimation_occ_disparity    ();
    void                estimation_occ_depth        ();
    void                clear_occ                   ();
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
    void                clear_error                 (CostType **error);
#endif
    //Owieczka end

    // Debug tool
    bool                getDistrotion               (DepthType **depth, int d, int m);

protected:

private:
    bool                readCameraParam             (FILE *fp, const char *target_id, CvMat *inMat, CvMat *exMat);
    void                convertCameraParam          (CvMat *exMat_dst, CvMat *exMat_src);
    void                image2world_with_z          (CvMat *mat_Rc2w_invIN_from, CvMat *matEX_c2w_from, CvMat *image, CvMat *world);

    void                caloc_z_range               (double focal_length, double baseline, double disparity_offset, CvMat* matEx_c2w, int iMinDisparity, int iMaxDisparity, unsigned int uiDepthType);
    bool                init_label2depth            (double focal_length, double baseline, double disparity_offset, CvMat* matEx_c2w, int iMinDisparity, int iMaxDisparity, int iMinSearch, int iMaxSearch, unsigned int uiDepthType);
#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
    bool                init_label2depth_z          (int iNumberOfDepthSteps, double dNearestSearchDepthValue, double dFarthestSearchDepthValue, double dNearestDepthValue, double dFarthestDepthValue, unsigned int uiDepthType);
#endif
    bool                init_label2disparity_simple (double ratio[2], int iMinSearch, int iMaxSearch);
    bool                init_label2disparity_from_z (double focal_length, double baseline[2], double disparity_offset[2], CvMat* matEx_c2w, unsigned int uiDepthType);
    bool                init_homography             (CvMat *matIn[3], CvMat *matEx_c2w[3], unsigned int uiDepthType);

    bool                setup_from_search_range     (int iMinSearchRnage, int iMaxSearchRange, int iMinDisparityRange, int iMaxDisparityRange,
                                                     CvMat *matIn[3], CvMat *matEx_c2w[3], unsigned int uiDepthType, unsigned int uiBaselineBasis, unsigned int uiMatchingMethod);
#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
    bool                setup_from_z_search_range   (int iNumberOfDepthSteps, double dNearestSearchDepthValue, double dFarthestSearchDepthValue, double dNearestDepthValue, double dFarthestDepthValue, 
                                                     CvMat *matIn[3], CvMat *matEx_c2w[3], unsigned int uiDepthType, unsigned int uiBaselineBasis, unsigned int uiMatchingMethod);
#endif
    void                print_range                 (double dSearchLevel, int iMinSearchRnage, int iMaxSearchRange, int iMinDisparityRange, int iMaxDisparityRange);
#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
    void                print_z_range               (int iDepthSteps, double dNearestSearchDepthValue, double dFarthestSearchDepthValue, double dNearestDepthValue, double dFarthestDepthValue);
#endif

    //Nagoya start
    void                block_matching_homography_pixel     (CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM);
    void                block_matching_homography_window    (CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM);

    void                block_matching_disparity_pixel      (CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM);
    void                block_matching_disparity_window     (CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM);
    //Nagoya end
    void                block_matching_disparity_soft       (CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM); // Poznan
    int                 update_error_cost_etri              (CParameterDepthEstimation *cParameter, IplImage* ipl_Ycurr); // added by MW

    int                 m_iWidth;
    int                 m_iHeight;
    int                 m_iPicsize;
    int                 m_iPrecision;    // Precion �� SearchLevel�� ���̴� ����?
    double              m_dSearchLevel; //Owieczka - support fast rought search

    int                 m_iVerticalPrecision;

    int                 m_iMatchingBlock; //Added by SZK
    int                 m_iImageSegmentation; //Added by SZK
    int                 m_iSegmentationMethod; //Added by SZK
    int                 m_iMaxCluster; //Added by SZK
    //Poznan start
    double              m_dSoftDistanceCoeff;
    double              m_dSoftColorCoeff;
    int                 m_uiSoftBlockWidth;
    int                 m_uiSoftBlockHeight;
   //Poznan end

    int                 m_iNumOfLabels;  // dispairty range ũ�� = maxdisparity - mindisparity 

    CvMat               **m_matH_V2L;
    CvMat               **m_matH_V2R;
    int                 *m_aiLabel2Disparity[2];
    DepthType           *m_acLabel2Depth;
#ifdef POZNAN_DYNAMIC_TABLE
    int                 *depth2label;
#else
    int                 depth2label[MAX_DEPTH];
#endif

    double              m_dZnear;
    double              m_dZfar;

    double              m_dLambda;
    double              m_dRho; //Added by SZK

    int                 m_iMaxWidth;
    int                 m_iMaxHeight;
    int                 m_iWidth_minus1;
    int                 m_iHeight_minus1;
    int                 m_iNumOfSegms; //Added by SZK

    //Owieczka
#ifdef POZNAN_DYNAMIC_TABLE
    double              *m_aiEdgeCost; //int --> double
    double              *m_aiEdgeCost2; //Added by SZK
    int                 *buf_byte_abs;
#else
    double              m_aiEdgeCost[MAX_DEPTH]; //int --> double
    double              m_aiEdgeCost2[MAX_DEPTH]; //Added by SZK
    int                 buf_byte_abs[MAX_DEPTH*2-1];
#endif
    int                 *byte_abs;

    Graph::node_id      *nodes;
    Graph::node_id      *auxiliary;
#ifdef POZNAN_TWO_COSTS
    CostType                 **errors_left;
    CostType                 **errors_right;
#else
    CostType                 **errors;
#endif
    DepthType           *labels;

    long double         **matrix; //Added by SZK
    long double         **Dmatrix; //Added by SZK
    long double         **inv_matrix; //Added by SZK
    long double         **prod_matrix; //Added by SZK
    long double         **check_matrix; //Added by SZK
    unsigned int        *segmlabel; //Added by SZK
//  int                 *statistics_label;
    int                 *reliability_matrix;

//  double              *variance; //Added by SZK
//  double              *average; //Added by SZK

    DepthType            *labels_prev; // GIST

    int                 **PrevTempB;  //ETRI

    void                (CEstimation::*func_block_matching)(CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM);
    void                getMotionMap_block(IplImage* src1, IplImage* src2, IplImage* dest, int block, double th);

    //Nagoya start
    IplImage*           ipl_manedge;
    IplImage*           ipl_mandisp;
    IplImage*           ipl_manstatic;
#ifdef POZNAN_OCC
    IplImage*           ipl_manoccl; //Owieczka
    IplImage*           ipl_manoccr; //Owieczka
    DepthType           *disparity_map_left;
    DepthType           *disparity_map_right;
#endif
    IplImage*           ipl_Yprev;
    IplImage*           ipl_Drefresh;
    BYTE                gc_mode;
    double              m_dOffset;
    double              m_dRefRatio;
    double              m_dRefBaseline;
    double              m_dRefFocal_length;

    //Nagoya end

};

#endif // _INCLUDE_ESTIMATION_H_
