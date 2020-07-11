#include "Estimation.h"
#include "ParameterDepthEstimation.h"
#include "yuv.h"
#include <time.h>
#include "PostProcessing.h" //ETRI
#include "CudaCuts.h"
#include "push_relabel_graphcut.hpp"
#include <vector>

#ifndef SAFE_RELEASE_IMAGE
#define SAFE_RELEASE_IMAGE(p) { if((p)!=NULL){ cvReleaseImage(&(p)); (p)=NULL; } }
#endif
#ifndef SAFE_RELEASE_MAT
#define SAFE_RELEASE_MAT(p) { if((p)!=NULL){ cvReleaseMat(&(p)); (p)=NULL; } }
#endif
#ifndef SAFE_FREE
#define SAFE_FREE(p) { if((p)!=NULL){ free(p); (p)=NULL; } }
#endif

#define MINSEL(a,b) (((a)<(b))?(a):(b))

CEstimation::CEstimation()
{
#ifdef POZNAN_DYNAMIC_TABLE
    depth2label = NULL;
    m_aiEdgeCost  = NULL;
    m_aiEdgeCost2 = NULL;
    buf_byte_abs = NULL;
#endif
    m_matH_V2L = NULL;
    m_matH_V2R = NULL;
    m_aiLabel2Disparity[0] = NULL;
    m_acLabel2Depth = NULL;
    nodes = NULL;
    auxiliary = NULL;
#ifdef POZNAN_TWO_COSTS
    errors_left = NULL; //Owieczka
    errors_right = NULL; //Owieczka
#else
    errors = NULL;
#endif
    labels = NULL;
    //Nagoya start
    segmlabel = NULL;
    reliability_matrix = NULL;

    ipl_manedge = NULL;
    ipl_mandisp = NULL;
    ipl_manstatic = NULL;
#ifdef POZNAN_OCC
    ipl_manoccl = NULL; //Owieczka
    ipl_manoccr = NULL;
    disparity_map_left = NULL;
    disparity_map_right = NULL;
#endif
    ipl_Yprev = NULL;

    //Nagoya end
    PrevTempB = NULL; // ETRI
    // GIST start
    labels_prev = NULL;
    // GIST end
}

CEstimation::~CEstimation()
{
    free_memory();
}

void CEstimation::free_memory()
{
#ifdef POZNAN_DYNAMIC_TABLE
    if(depth2label!=NULL) //Owieczka
    {
        free(depth2label);
        depth2label = NULL;
    }
    if(m_aiEdgeCost!=NULL)
    {
        free(m_aiEdgeCost);
        m_aiEdgeCost = NULL;
    }
    if(m_aiEdgeCost2!=NULL)
    {
        free(m_aiEdgeCost2);
        m_aiEdgeCost2 = NULL;
    }
    if(buf_byte_abs!=NULL)
    {
        free(buf_byte_abs);
        buf_byte_abs = NULL;
    }
#endif
    if(m_matH_V2L!=NULL)
    {
        for(int i=0; i<m_iNumOfLabels; i++)
        {
            cvReleaseMat(&m_matH_V2L[i]);
            cvReleaseMat(&m_matH_V2R[i]);
        }
        free(m_matH_V2L);
        m_matH_V2L=m_matH_V2R=NULL;
    }
    if(m_aiLabel2Disparity[0]!=NULL)
    {
        free(m_aiLabel2Disparity[0]);
        m_aiLabel2Disparity[0]=NULL;
    }
    if(m_acLabel2Depth!=NULL)
    {
        free(m_acLabel2Depth);
        m_acLabel2Depth=NULL;
    }
    if(nodes!=NULL)
    {
        free(nodes);
        nodes=NULL;
    }
    if(auxiliary!=NULL)
    {
        free(auxiliary);
        auxiliary=NULL;
    }
#ifdef POZNAN_TWO_COSTS
    if((errors_left!=NULL)&&(errors_right!=NULL))
#else
    if(errors!=NULL)
#endif
    {
#ifdef POZNAN_ERRORS_DISTRIBUTED_ALLOC
#ifdef POZNAN_TWO_COSTS
        for(int i=0; i<m_iNumOfLabels; i++)
        {
            //SAFE_FREE(errors_left[i]);
            if(errors_left[i]!=NULL)
            {
                free(errors_left[i]);
                errors_left[i]= NULL;
            }
            if(errors_right[i]!=NULL)
            {
                free(errors_right[i]);
                errors_right[i]= NULL;
            }
        }
#else
        for(int i=0; i<m_iNumOfLabels; i++)  // �� �Ѳ����� �Ҵ��� �ȹ���? 
        {
            if(errors[i]!=NULL) 
            {
                free(errors[i]);
                errors[i]= NULL;
            }
        }
#endif
#else
#ifdef POZNAN_TWO_COSTS
        if(errors_left[0]!=NULL) free(errors_left[0]);
        free(errors_left);
        errors_left=NULL;
        if(errors_right[0]!=NULL) free(errors_right[0]);
        free(errors_right);
        errors_right=NULL;
#else
        if(errors[0]!=NULL) free(errors[0]);
        free(errors);
        errors=NULL;
#endif
#endif
    }
    if(labels!=NULL)
    {
        free(labels);
        labels=NULL;
    }

    //Nagoya start
    if(segmlabel != NULL)
    {
//      if(segmlabel[0] != NULL) free(segmlabel[0]);
        free(segmlabel);
        segmlabel = NULL;
    }
    SAFE_RELEASE_IMAGE(ipl_manedge)
    SAFE_RELEASE_IMAGE(ipl_mandisp)
    SAFE_RELEASE_IMAGE(ipl_manstatic)
#ifdef POZNAN_OCC
    SAFE_RELEASE_IMAGE(ipl_manoccr) //Owieczka
    SAFE_RELEASE_IMAGE(ipl_manoccl) //Owieczka
    SAFE_FREE(disparity_map_left)
    SAFE_FREE(disparity_map_right)
#endif
    SAFE_RELEASE_IMAGE(ipl_Yprev)
    //Nagoya end

    //ETRI start
    if(PrevTempB != NULL) {
        for (int i=0; i<m_iHeight; i++) {
            delete[] PrevTempB[i];
            PrevTempB[i] = NULL;
        }
        delete[] PrevTempB;
        PrevTempB = NULL;
    }
    //ETRI end

    // GIST start
    if(labels_prev!=NULL)
    {
        free(labels_prev);
        labels_prev=NULL;
    }
    // GIST end

}

bool CEstimation::xInit(CParameterDepthEstimation cParameter)
{
    int i, pos;
    FILE *fp;
    CvMat* matIn[4];        // intrinsic parameter of left/right/center camera 3x3 matrix
    CvMat* matEx_c2w[4];    // extrinsic parameter of left/right/center camera 3x4 matrix
    free_memory();

    m_iHeight = m_iMaxHeight = int(cParameter.getSourceHeight());
    m_iWidth = m_iMaxWidth = int(cParameter.getSourceWidth());
    m_iPicsize = m_iHeight*m_iWidth;
    m_iHeight_minus1 = m_iHeight - 1;
    m_iWidth_minus1 = m_iWidth - 1;
    m_dLambda = cParameter.getSmoothCoeff();
    m_dRho    = cParameter.getSmoothCoeff2(); //Nagoya
    gc_mode   = 0; // Nagoya

#ifdef SUB_PEL_PRECISION
    m_iPrecision = int(cParameter.getPrecision());
    m_dSearchLevel = double(cParameter.getSearchLevel());
    m_iMaxWidth *= m_iPrecision;
#ifdef SUB_PEL_VERTICAL_PRECISION
    m_iVerticalPrecision = int(cParameter.getVerticalPrecision());
    m_iMaxHeight *= m_iVerticalPrecision;
#endif
#endif
    //Soft-Segmentation by Owieczka - Krzysztof Wegner - PUT Poznan
    m_dSoftDistanceCoeff = cParameter.getSoftDistanceCoeff();
    m_dSoftColorCoeff    = cParameter.getSoftColorCoeff();
    m_uiSoftBlockWidth   = int(cParameter.getSoftBlockWidth());
    m_uiSoftBlockHeight  = int(cParameter.getSoftBlockHeight());
    //Poznan end

    m_iMatchingBlock = cParameter.getMatchingBlock(); //Nagoya
    m_iImageSegmentation = cParameter.getImageSegmentation(); //Nagoya
    m_iSegmentationMethod = cParameter.getSegmentationMethod(); //Nagoya
    m_iMaxCluster = cParameter.getMaxCluster(); //Nagoya

    // read camera parameters
    int NofCameraParam = 3;
    if (cParameter.getDEmode() == 3)
      NofCameraParam = 4;

    if((fp = fopen(cParameter.getFileCamPara().c_str(), "r"))==NULL)
    {
        fprintf(stderr, "Can't open the file of camera parameters\n");
        return false;
    }
    for(i=0; i<NofCameraParam; i++)
    {
        matIn[i] = cvCreateMat(3, 3, CV_64F); //intrinsic parameter of left/right/virtual camera 3x3 matrix
        matEx_c2w[i] = cvCreateMat(3, 4, CV_64F); //extrinsic parameter of left/right/virtual camera 3x4 matrix
    }
	// �� ������ Left-Right-Centor? �׸��� ��Ī�� �� Virtual?
	// getDEmode == 3�̸� depth ī�޶�?
    if( !readCameraParam(fp, cParameter.getLeftCameraName().c_str(), matIn[0], matEx_c2w[0]) ||
        !readCameraParam(fp, cParameter.getRightCameraName().c_str(), matIn[1], matEx_c2w[1]) ||
        !readCameraParam(fp, cParameter.getCenterCameraName().c_str(), matIn[2], matEx_c2w[2]) ||
        ((cParameter.getDEmode() == 3) && !readCameraParam(fp, cParameter.getRefDepthCameraName().c_str(), matIn[3], matEx_c2w[3]))

    )
    {
        fprintf(stderr, "Can't get camera parameter(s)\n");
        fclose(fp);
        for(i=0; i<NofCameraParam; i++)
        {
            cvReleaseMat(&matIn[i]);
            cvReleaseMat(&matEx_c2w[i]);
        }
        return false;
    }
    fclose(fp);

    if(cParameter.getBaselineBasis() < 4)
    {
#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
        if( !setup_from_z_search_range(cParameter.getNumberOfDepthSteps(), cParameter.getNearestSearchDepthValue(), cParameter.getFarthestSearchDepthValue(), cParameter.getNearestDepthValue(), cParameter.getFarthestDepthValue(), matIn, matEx_c2w, cParameter.getDepthType(), cParameter.getBaselineBasis(), cParameter.getMatchingMethod()) )
#else
        if( !setup_from_search_range(cParameter.getMinValDisSearchRange(), cParameter.getMaxValDisSearchRange(), cParameter.getMinValDisRange(), cParameter.getMaxValDisRange(), matIn, matEx_c2w, cParameter.getDepthType(), cParameter.getBaselineBasis(), cParameter.getMatchingMethod()) )
#endif
        {
            for(i=0; i<NofCameraParam; i++)
            {
                cvReleaseMat(&matIn[i]);
                cvReleaseMat(&matEx_c2w[i]);
            }
            return false;
        }
    }

    if (cParameter.getDEmode() == 3) {
        double tmp0 = cvmGet(matEx_c2w[2], 0, 3) - cvmGet(matEx_c2w[0], 0, 3);
        double tmp1 = cvmGet(matEx_c2w[2], 0, 3) - cvmGet(matEx_c2w[1], 0, 3);

        m_dRefBaseline     = cvmGet(matEx_c2w[2], 0, 3)-cvmGet(matEx_c2w[3], 0, 3);
        if( m_dRefBaseline <= tmp1 )
            m_dRefRatio = fabs(tmp1/m_dRefBaseline);
        else
            m_dRefRatio = fabs(tmp0/m_dRefBaseline);

        m_dOffset          = cvmGet(matIn[2], 0, 2)-cvmGet(matIn[3], 0, 2);
        m_dRefFocal_length = cvmGet(matIn[3], 0, 0);
    }

    for(i=0; i<NofCameraParam; i++)
    {
        cvReleaseMat(&matIn[i]);
        cvReleaseMat(&matEx_c2w[i]);
    }

    // memory allocation
    if( (nodes = (Graph::node_id *)malloc(m_iPicsize*sizeof(Graph::node_id))) == NULL ||
            (auxiliary = (Graph::node_id *)malloc((m_iPicsize<<1)*sizeof(Graph::node_id))) == NULL  ||
#ifdef POZNAN_TWO_COSTS
            (errors_left = (CostType **)malloc(m_iNumOfLabels*sizeof(CostType *))) == NULL ||
            (errors_right = (CostType **)malloc(m_iNumOfLabels*sizeof(CostType *))) == NULL ||
#else
            (errors = (CostType **)malloc(m_iNumOfLabels*sizeof(CostType *))) == NULL ||
#endif
#ifdef POZNAN_OCC
            (disparity_map_left = (DepthType *)malloc(m_iPicsize*sizeof(DepthType)))== NULL ||
            (disparity_map_right = (DepthType *)malloc(m_iPicsize*sizeof(DepthType)))== NULL ||
#endif
            (labels = (DepthType *)malloc(m_iPicsize*sizeof(DepthType))) == NULL ||
            (segmlabel = (unsigned int *)malloc(m_iWidth * m_iHeight*sizeof(unsigned int))) == NULL
            // GIST start
            || (labels_prev = (DepthType *)malloc(m_iPicsize*sizeof(DepthType))) == NULL
            // GIST end
        )
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
#ifdef POZNAN_ERRORS_DISTRIBUTED_ALLOC
#ifdef POZNAN_TWO_COSTS
    for(i=0, pos=m_iPicsize; i<m_iNumOfLabels; i++, pos+=m_iPicsize)
    {
        if(((errors_left[i]      = (CostType *)malloc(m_iPicsize*sizeof(CostType))) == NULL)
        ||( (errors_right[i]     = (CostType *)malloc(m_iPicsize*sizeof(CostType))) == NULL))
        {
            fprintf(stderr, "Can't allocate enough memory error row %d\n",i);
            return false;
        }
    }
#else
    for(i=0, pos=m_iPicsize; i<m_iNumOfLabels; i++, pos+=m_iPicsize) //�� �Ѳ����� �Ҵ��� ���ұ��?
    {
        if( (errors[i]      = (CostType *)malloc(m_iPicsize*sizeof(CostType))) == NULL)
        {
            fprintf(stderr, "Can't allocate enough memory error row %d\n",i);
            return false;
        }
    }
#endif
#else
#ifdef POZNAN_TWO_COSTS
    if( (errors_left[0]      = (CostType *)malloc(m_iNumOfLabels*m_iPicsize*sizeof(CostType))) == NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
    if( (errors_right[0]      = (CostType *)malloc(m_iNumOfLabels*m_iPicsize*sizeof(CostType))) == NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }

    for(i=1, pos=m_iPicsize; i<m_iNumOfLabels; i++, pos+=m_iPicsize)
    {
        errors_left[i]      = &errors_left[0][pos];
        errors_right[i]      = &errors_right[0][pos];
    }    
#else
    if( (errors[0]      = (CostType *)malloc(m_iNumOfLabels*m_iPicsize*sizeof(CostType))) == NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }

    for(i=1, pos=m_iPicsize; i<m_iNumOfLabels; i++, pos+=m_iPicsize)
    {
        errors[i]      = &errors[0][pos];
    }    
#endif
#endif

    // Nagoya start
    ipl_Yprev    = cvCreateImage(cvSize(m_iWidth, m_iHeight), IPL_DEPTH_8U, 1);
    ipl_Drefresh = cvCreateImage(cvSize(m_iWidth, m_iHeight), IPL_DEPTH_8U, 1);
    if( ipl_Yprev == NULL || ipl_Drefresh == NULL) {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
    // Nagoya end

    //ETRI start
    if (cParameter.getDEmode() == 2) {
#ifdef SUB_PEL_PRECISION
        int new_width  = m_iWidth * cParameter.getPrecision();
#else
        int new_width = m_iWidth;
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
        int new_height = m_iHeight * cParameter.getVerticalPrecision();
#else
        int new_height = m_iHeight;
#endif

        PrevTempB = new int *[new_height];
        for (int i=0; i<new_height; i++) {
            PrevTempB[i] = new int [new_width];
        }
        for (int j=0; j<new_height; j++) {
            for (int i=0; i<new_width; i++) {
                PrevTempB[j][i] = 0;
            }
        }
    }
    //ETRI end

#ifdef POZNAN_DYNAMIC_TABLE
    m_aiEdgeCost  = (double*) malloc(MAX_DEPTH*sizeof(double));
    m_aiEdgeCost2 = (double*) malloc(MAX_DEPTH*sizeof(double));
    buf_byte_abs = (int*) malloc((MAX_DEPTH*2-1)*sizeof(int));
#endif
    byte_abs=&buf_byte_abs[MAX_DEPTH-1];
    for(i=0; i<MAX_DEPTH; i++)
    {
        byte_abs[-i]=byte_abs[i]=i;
        m_aiEdgeCost[i]=(int)(m_dLambda*i);
        m_aiEdgeCost2[i]=(int)(m_dRho * m_dLambda * i); //Nagoya
    }

    return true;
}


void CEstimation::print_range(double dSearchLevel, int iMinSearchRange, int iMaxSearchRange, int iMinDisparityRange, int iMaxDisparityRange)
{
#ifdef SUB_PEL_PRECISION
    if (dSearchLevel==1)
    {
        printf("Disparity step               : INTEGER PIXEL\n");
        printf("Disparity Search Range       : [%d:%d] pel\n", iMinSearchRange, iMaxSearchRange);
        printf("Disparity Range on Depth Map : [%d:%d] pel\n", iMinDisparityRange, iMaxDisparityRange);
    }
    else if (dSearchLevel==2)
    {
        printf("Disparity step               : HALF PIXEL\n");
        printf("Disparity Search Range       : [%.1f:%.1f] pel\n", ((double)iMinSearchRange)/((double)dSearchLevel), ((double)iMaxSearchRange)/((double)dSearchLevel));
        printf("Disparity Range on Depth Map : [%.1f:%.1f] pel\n", ((double)iMinDisparityRange)/((double)dSearchLevel), ((double)iMaxDisparityRange)/((double)dSearchLevel));
    }
    else if (dSearchLevel==4)
    {
        printf("Disparity step               : QUARTER PIXEL\n");
        printf("Disparity Search Range       : [%.2f:%.2f] pel\n", ((double)iMinSearchRange)/((double)dSearchLevel), ((double)iMaxSearchRange)/((double)dSearchLevel));
        printf("Disparity Range on Depth Map : [%.2f:%.2f] pel\n", ((double)iMinDisparityRange)/((double)dSearchLevel), ((double)iMaxDisparityRange)/((double)dSearchLevel));
    }
    else if (dSearchLevel==0.5)
    {
        printf("Disparity step               : 2 PIXEL\n");
        printf("Disparity Search Range       : [%.2f:%.2f] pel\n", ((double)iMinSearchRange)/((double)dSearchLevel), ((double)iMaxSearchRange)/((double)dSearchLevel));
        printf("Disparity Range on Depth Map : [%.2f:%.2f] pel\n", ((double)iMinDisparityRange)/((double)dSearchLevel), ((double)iMaxDisparityRange)/((double)dSearchLevel));
    }
    else if (dSearchLevel==0.25)
    {
        printf("Disparity step               : 4 PIXEL\n");
        printf("Disparity Search Range       : [%.2f:%.2f] pel\n", ((double)iMinSearchRange)/((double)dSearchLevel), ((double)iMaxSearchRange)/((double)dSearchLevel));
        printf("Disparity Range on Depth Map : [%.2f:%.2f] pel\n", ((double)iMinDisparityRange)/((double)dSearchLevel), ((double)iMaxDisparityRange)/((double)dSearchLevel));
    }
    else
    {
        fprintf(stderr, "Unknown number on SearchLevel [%lf]\n", dSearchLevel);
        return;
    }
#else
    printf("Disparity step               : INTEGER PIXEL\n");
    printf("Disparity Search Range       : [%d:%d] pel\n", iMinSearchRange, iMaxSearchRange);
    printf("Disparity Range on Depth Map : [%d:%d] pel\n", iMinDisparityRange, iMaxDisparityRange);
#endif
}

#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
void CEstimation::print_z_range(int iDepthSteps, double dNearestSearchDepthValue, double dFarthestSearchDepthValue, double dNearestDepthValue, double dFarthestDepthValue)
{
    printf("Z steps              : %d\n",iDepthSteps);
    printf("Z Search Range       : [%f:%f] \n", dNearestSearchDepthValue, dFarthestSearchDepthValue);
    printf("Z Range on Depth Map : [%f:%f] \n", dNearestDepthValue, dFarthestDepthValue);
}
#endif

bool CEstimation::setup_from_search_range(int iMinSearchRange, int iMaxSearchRange, int iMinDisparityRange, int iMaxDisparityRange, CvMat *matIn[3], CvMat *matEx_c2w[3], unsigned int uiDepthType, unsigned int uiBaselineBasis, unsigned int uiMatchingMethod)
{
    double baseline[2], ratio[2], disparity_offset[2];

#ifdef SUB_PEL_PRECISION
    iMinSearchRange *= m_dSearchLevel;   // ? Precision vs Search level (1: integer, 2: half-pel, 4: quarter-pel
    iMaxSearchRange *= m_dSearchLevel;
    iMinDisparityRange *= m_dSearchLevel;
    iMaxDisparityRange *= m_dSearchLevel;
#endif

	// 255 ���� ū MaxDisparity �� ��� �̷��� ǥ�� �ϴ� ������ ����.
    if(iMaxDisparityRange < iMinDisparityRange)
    {
        printf("*** AUTOMATICAL PARAMETER MODIFICATION ***\n");
#ifdef SUB_PEL_PRECISION
        printf("  MaximumValueOfDisparityRange : %.2f -> %.2f\n", double(iMaxDisparityRange)/double(m_dSearchLevel), double(iMinDisparityRange + 255)/double(m_dSearchLevel));
#else
        printf("  MaximumValueOfDisparityRange : %d -> %d\n", iMaxDisparityRange, iMinDisparityRange + 255);
#endif
        printf("  MaximumValueOfDisparityRange must be larger than or equal to MinimumValueOfDisparityRange\n");
        iMaxDisparityRange = iMinDisparityRange + MAX_DEPTH-1;
    }

	// �ּ� �˻� ���� ���� ����   
    if(iMinSearchRange < iMinDisparityRange)
    {
        printf("*** AUTOMATICAL PARAMETER MODIFICATION ***\n");
#ifdef SUB_PEL_PRECISION
        printf("  MinimumValueOfdisparitySearchRange : %.2f -> %.2f\n", double(iMinSearchRange)/double(m_dSearchLevel), double(iMinDisparityRange)/double(m_dSearchLevel));
#else
        printf("  MinimumValueOfdisparitySearchRange : %d -> %d\n", iMinSearchRange, iMinDisparityRange);
#endif
        printf("  MinimumValueOfDisparitySearchRange must be larger than or equal to MinimumValueOfDisparityRange\n");
        iMinSearchRange = iMinDisparityRange;
    }
	// �ִ� �˻� ���� ���� ���� 
    if(iMaxSearchRange > iMaxDisparityRange)
    {
        printf("*** AUTOMATICAL PARAMETER MODIFICATION ***\n");
#ifdef SUB_PEL_PRECISION
        printf("  MaximumValueOfDisparitySearchRange : %.2f -> %.2f\n", double(iMaxSearchRange)/double(m_dSearchLevel), double(iMaxDisparityRange)/double(m_dSearchLevel));
#else
        printf("  MaximumValueOfDisparitySearchRange : %.2f -> %.2f\n", iMaxSearchRange, iMaxDisparityRange);
#endif
        printf("  MaximumValueOfDisparitySearchRange must be less than or equal to MaximumValueOfDisparityRange\n");
        iMaxSearchRange = iMaxDisparityRange;
    }

	// ���� ��ġ ������ ǥ�� ���� 
    print_range(m_dSearchLevel, iMinSearchRange, iMaxSearchRange, iMinDisparityRange, iMaxDisparityRange);

    m_iNumOfLabels = iMaxSearchRange - iMinSearchRange + 1;  // �˻��� disparity �� (subpixel ����)
    if((iMaxDisparityRange - iMinDisparityRange + 1)>MAX_DEPTH)
    {
        fprintf(stderr, "The number of depth candidates is too large.\n");
        return false;
    }

	// ��������: BaselineBasis
    // set candidates of disparity
    baseline[0] = cvmGet(matEx_c2w[2], 0, 3)-cvmGet(matEx_c2w[0], 0, 3);  // Center.x - Left.x = left BL
    baseline[1] = cvmGet(matEx_c2w[1], 0, 3)-cvmGet(matEx_c2w[2], 0, 3);  // Right.x  - Center.x  = right baseline
    switch(uiBaselineBasis)
    {
    case 0: // minimum baseline
        if(fabs(baseline[0])<fabs(baseline[1]))
            uiBaselineBasis = 0;
        else
            uiBaselineBasis = 1;
        break;
    case 1: // maximum baseline
        if(fabs(baseline[0])>fabs(baseline[1]))
            uiBaselineBasis = 0;
        else
            uiBaselineBasis = 1;
        break;
    case 2: // left baseline
        uiBaselineBasis = 0;
        break;
    case 3: // right baseline
        uiBaselineBasis = 1;
        break;
    default:
        fprintf(stderr, "Unknown number on BaselineBasis [%d]\n", uiBaselineBasis);
        return false;
    }

    ratio[0] = baseline[0]/baseline[uiBaselineBasis];
    ratio[1] = baseline[1]/baseline[uiBaselineBasis];
    disparity_offset[0] = cvmGet(matIn[2], 0, 2) - cvmGet(matIn[0], 0, 2); // intrinsic, Center.offset - Left.offset
    disparity_offset[1] = cvmGet(matIn[1], 0, 2) - cvmGet(matIn[2], 0, 2); //  

	// (m_dZnear, m_dZfar) ��� 
    caloc_z_range(cvmGet(matIn[2], 0, 0), baseline[uiBaselineBasis], disparity_offset[uiBaselineBasis], matEx_c2w[2], iMinDisparityRange, iMaxDisparityRange, uiDepthType);

    if(!init_label2depth(cvmGet(matIn[2], 0, 0), baseline[uiBaselineBasis], disparity_offset[uiBaselineBasis], matEx_c2w[2], iMinDisparityRange, iMaxDisparityRange, iMinSearchRange, iMaxSearchRange, uiDepthType))
        return false;

#ifdef ACCURATE_CORRESPONDENCE
    switch(uiMatchingMethod)
    {
    case 0: // not defined
        if(!init_label2disparity_simple(ratio, iMinSearchRange, iMaxSearchRange))
            return false;
#ifdef POZNAN_OCC
        func_data_term = &CEstimation::data_term_disparity;
        func_estimation_occ = &CEstimation::estimation_occ_disparity;
#endif
        if(m_iMatchingBlock==3)
        func_block_matching = &CEstimation::block_matching_disparity_window;
        else
        func_block_matching = &CEstimation::block_matching_disparity_pixel;
        break;
    case 1: // disparity-based
        if(!init_label2disparity_from_z(cvmGet(matIn[2], 0, 0), baseline, disparity_offset, matEx_c2w[2], uiDepthType))
            return false;
#ifdef POZNAN_OCC
        func_data_term = &CEstimation::data_term_disparity;
        func_estimation_occ = &CEstimation::estimation_occ_disparity;
#endif
        if(m_iMatchingBlock==3)
        func_block_matching = &CEstimation::block_matching_disparity_window;
        else
        func_block_matching = &CEstimation::block_matching_disparity_pixel;
        break;
    case 2: // homography-based
    case 4: // epipolar_line-based
        if(!init_homography(matIn, matEx_c2w, uiDepthType))
            return false;
#ifdef POZNAN_OCC
        func_data_term = &CEstimation::data_term_depth;
        func_estimation_occ = &CEstimation::estimation_occ_depth;
        //if(!init_label2disparity_from_z(cvmGet(matIn[2], 0, 0), baseline, disparity_offset, matEx_c2w[2], uiDepthType))
        //    return false;
#endif
        /*
        if(m_iImageSegmentation==1 && m_iMatchingBlock==1)
//      func_block_matching = &CEstimation::block_matching_homography_segmentation_pixel;
        func_block_matching = &CEstimation::block_matching_homography_pixel;
        else if (m_iImageSegmentation==1 && m_iMatchingBlock==3)
//      func_block_matching = &CEstimation::block_matching_homography_segmentation_window;
        func_block_matching = &CEstimation::block_matching_homography_window;
        else if (m_iImageSegmentation==0 && m_iMatchingBlock==3)
        func_block_matching = &CEstimation::block_matching_homography_window;
        else if (m_iImageSegmentation==0 && m_iMatchingBlock==1)
        func_block_matching = &CEstimation::block_matching_homography_pixel;
        */
        //Nagoya start
        if(m_iMatchingBlock==3)
        func_block_matching = &CEstimation::block_matching_homography_window;
        else
        func_block_matching = &CEstimation::block_matching_homography_pixel;
        //Nagoya end
        break;
    //Poznan start
    case 3: // Soft Segment
        if(!init_label2disparity_simple(ratio, iMinSearchRange, iMaxSearchRange))
            return false;
        func_block_matching = &CEstimation::block_matching_disparity_soft;
#ifdef POZNAN_OCC
        func_data_term = &CEstimation::data_term_disparity;
        func_estimation_occ = &CEstimation::estimation_occ_disparity;
#endif
    //Poznan end
        break;
    default:
        fprintf(stderr, "Unknown value on MatchingMethod\n");
        return false;
        break;
    }
#else
    if(!init_label2disparity_simple(ratio, iMinSearchRange, iMaxSearchRange))
        return false;
    func_block_matching = &CEstimation::block_matching_disparity;
#endif

    return true;
}

#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
bool CEstimation::setup_from_z_search_range(int iNumberOfDepthSteps, double dNearestSearchDepthValue, double dFarthestSearchDepthValue, double dNearestDepthValue, double dFarthestDepthValue, CvMat *matIn[3], CvMat *matEx_c2w[3], unsigned int uiDepthType, unsigned int uiBaselineBasis, unsigned int uiMatchingMethod)
{
    m_dZnear = dNearestDepthValue;
    m_dZfar  = dFarthestDepthValue;

    /* Owieczka Fix - Ensure that ZSearchRange is inside ZRange
    if(iMinSearchRange < iMinDisparityRange)
    {
        printf("*** AUTOMATICAL PARAMETER MODIFICATION ***\n");
#ifdef SUB_PEL_PRECISION
        printf("  MinimumValueOfdisparitySearchRange : %.2f -> %.2f\n", double(iMinSearchRange)/double(m_dSearchLevel), double(iMinDisparityRange)/double(m_dSearchLevel));
#else
        printf("  MinimumValueOfdisparitySearchRange : %d -> %d\n", iMinSearchRange, iMinDisparityRange);
#endif
        printf("  MinimumValueOfDisparitySearchRange must be larger than or equal to MinimumValueOfDisparityRange\n");
        iMinSearchRange = iMinDisparityRange;
    }

    if(iMaxSearchRange > iMaxDisparityRange)
    {
        printf("*** AUTOMATICAL PARAMETER MODIFICATION ***\n");
#ifdef SUB_PEL_PRECISION
        printf("  MaximumValueOfDisparitySearchRange : %.2f -> %.2f\n", double(iMaxSearchRange)/double(m_dSearchLevel), double(iMaxDisparityRange)/double(m_dSearchLevel));
#else
        printf("  MaximumValueOfDisparitySearchRange : %.2f -> %.2f\n", iMaxSearchRange, iMaxDisparityRange);
#endif
        printf("  MaximumValueOfDisparitySearchRange must be less than or equal to MaximumValueOfDisparityRange\n");
        iMaxSearchRange = iMaxDisparityRange;
    }//*/

    //print_range(m_dSearchLevel, iMinSearchRange, iMaxSearchRange, iMinDisparityRange, iMaxDisparityRange);
    print_z_range(iNumberOfDepthSteps, dNearestSearchDepthValue, dFarthestSearchDepthValue, m_dZnear, m_dZfar);

    //m_iNumOfLabels = iMaxSearchRange - iMinSearchRange + 1;
    m_iNumOfLabels = iNumberOfDepthSteps;

    if(m_iNumOfLabels > MAX_DEPTH)
    {
        fprintf(stderr, "The number of depth candidates is too large.\n");
        return false;
    }

    //To sko?zy�em
    //caloc_z_range(cvmGet(matIn[2], 0, 0), baseline[uiBaselineBasis], disparity_offset[uiBaselineBasis], matEx_c2w[2], iMinDisparityRange, iMaxDisparityRange, uiDepthType);

    //if(!init_label2depth(cvmGet(matIn[2], 0, 0), baseline[uiBaselineBasis], disparity_offset[uiBaselineBasis], matEx_c2w[2], iMinDisparityRange, iMaxDisparityRange, iMinSearchRange, iMaxSearchRange, uiDepthType))
    if(!init_label2depth_z(iNumberOfDepthSteps, dNearestSearchDepthValue, dFarthestSearchDepthValue, dNearestDepthValue, dFarthestDepthValue, uiDepthType))
        return false;

#ifdef ACCURATE_CORRESPONDENCE
    switch(uiMatchingMethod)
    {
    case 0: // not defined
        //Fix Owieczka
        /*
        if(!init_label2disparity_simple(ratio, iMinSearchRange, iMaxSearchRange))
            return false;
#ifdef POZNAN_OCC
        func_data_term = &CEstimation::data_term_disparity;
        func_estimation_occ = &CEstimation::estimation_occ_disparity;
#endif
        if(m_iMatchingBlock==3)
        func_block_matching = &CEstimation::block_matching_disparity_window;
        else
        func_block_matching = &CEstimation::block_matching_disparity_pixel;
        break;
    case 1: // disparity-based
        if(!init_label2disparity_from_z(cvmGet(matIn[2], 0, 0), baseline, disparity_offset, matEx_c2w[2], uiDepthType))
            return false;
#ifdef POZNAN_OCC
        func_data_term = &CEstimation::data_term_disparity;
        func_estimation_occ = &CEstimation::estimation_occ_disparity;
#endif
        if(m_iMatchingBlock==3)
        func_block_matching = &CEstimation::block_matching_disparity_window;
        else
        func_block_matching = &CEstimation::block_matching_disparity_pixel;
        break;
        //*/
    case 2: // homography-based
    case 4: // epipolar_line-based
        if(!init_homography(matIn, matEx_c2w, uiDepthType))
            return false;
#ifdef POZNAN_OCC
        func_data_term = &CEstimation::data_term_depth;
        func_estimation_occ = &CEstimation::estimation_occ_depth;
        //if(!init_label2disparity_from_z(cvmGet(matIn[2], 0, 0), baseline, disparity_offset, matEx_c2w[2], uiDepthType))
        //    return false;
#endif
        /*
        if(m_iImageSegmentation==1 && m_iMatchingBlock==1)
//      func_block_matching = &CEstimation::block_matching_homography_segmentation_pixel;
        func_block_matching = &CEstimation::block_matching_homography_pixel;
        else if (m_iImageSegmentation==1 && m_iMatchingBlock==3)
//      func_block_matching = &CEstimation::block_matching_homography_segmentation_window;
        func_block_matching = &CEstimation::block_matching_homography_window;
        else if (m_iImageSegmentation==0 && m_iMatchingBlock==3)
        func_block_matching = &CEstimation::block_matching_homography_window;
        else if (m_iImageSegmentation==0 && m_iMatchingBlock==1)
        func_block_matching = &CEstimation::block_matching_homography_pixel;
        */
        //Nagoya start
        if(m_iMatchingBlock==3)
        func_block_matching = &CEstimation::block_matching_homography_window;
        else
        func_block_matching = &CEstimation::block_matching_homography_pixel;
        //Nagoya end
        break;
    //Poznan start
    case 3: // Soft Segment
        //Fix Owieczka
        /*
        if(!init_label2disparity_simple(ratio, iMinSearchRange, iMaxSearchRange))
            return false;
        func_block_matching = &CEstimation::block_matching_disparity_soft;
#ifdef POZNAN_OCC
        func_data_term = &CEstimation::data_term_disparity;
        func_estimation_occ = &CEstimation::estimation_occ_disparity;
#endif
    //Poznan end
        break;
    default:
        fprintf(stderr, "Unknown value on MatchingMethod\n");
        return false;
        //*/
        break;
    }
#else
    if(!init_label2disparity_simple(ratio, iMinSearchRange, iMaxSearchRange))
        return false;
    func_block_matching = &CEstimation::block_matching_disparity;
#endif

    return true;
}
#endif

//
// disparity Range �� baseline ī�޶� focal length�� ���� Znear/far ��� 
// Z = f * baseline / disparity 
// (baseline*f)/Z => disparity  
void CEstimation::caloc_z_range(double focal_length, double baseline, double disparity_offset, CvMat* matEx_c2w, int iMinDisparity, int iMaxDisparity, unsigned int uiDepthType)
{
#ifndef SUB_PEL_PRECISION
    m_dZnear = focal_length * baseline /
                                ( double(iMaxDisparity) + disparity_offset );
    m_dZfar  = focal_length * baseline /
                                ( double(iMinDisparity) + disparity_offset );
#else
    m_dZnear = focal_length * baseline /
                                ( double(iMaxDisparity)/double(m_dSearchLevel) + disparity_offset );
    m_dZfar  = focal_length * baseline /
                                ( double(iMinDisparity)/double(m_dSearchLevel) + disparity_offset );

    // m_dZfar/baseline =  focal_length/minDisparity   // ��κ� 0 
	// m_dZnear/baseline =  focal_length/maxDisparity  
	// 1 m
	//
	//  * Znear
	//  |   1m 
	//  | 
	//  L-cam      R-cam
	//   ����        ���� 
	//     ���� dispairty  
	//  1 m /0.3m =   focal_length   / 200
	//  200/0.3 = f = 600 ....??

#endif

    if(uiDepthType)
    {
        m_dZnear = cvmGet(matEx_c2w, 2, 2) * m_dZnear + cvmGet(matEx_c2w, 2, 3);  // ���� 0 �ƴѰ���?
        m_dZfar  = cvmGet(matEx_c2w, 2, 2) * m_dZfar  + cvmGet(matEx_c2w, 2, 3);
    }
}

bool CEstimation::init_label2depth(double focal_length, double baseline, double disparity_offset, CvMat* matEx_c2w, int iMinDisparity, int iMaxDisparity, int iMinSearch, int iMaxSearch, unsigned int uiDepthType)
{
    int l, iDisparity;
    if((m_acLabel2Depth = (DepthType *)malloc(m_iNumOfLabels*sizeof(DepthType)))==NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
    memset(m_acLabel2Depth, 0, m_iNumOfLabels*sizeof(DepthType));

    if(uiDepthType)
    {
        double rfl, mp;
        double dMaxDisparity, dMinDisparity, dCurDisparity;

        dMaxDisparity = double(iMaxDisparity);
        dMinDisparity = double(iMinDisparity);

#ifdef SUB_PEL_PRECISION
        dMaxDisparity /= double(m_dSearchLevel);
        dMinDisparity /= double(m_dSearchLevel);
#endif

        rfl = cvmGet(matEx_c2w, 2, 2) * focal_length * baseline;
        mp = (MAX_DEPTH-1) * ( rfl + (dMaxDisparity+disparity_offset)*cvmGet(matEx_c2w, 2, 3) ) / ( dMaxDisparity - dMinDisparity );
        for(l=0, iDisparity=iMinSearch; iDisparity<=iMaxSearch; l++, iDisparity++)
        {
            dCurDisparity = double(iDisparity);
#ifdef SUB_PEL_PRECISION
            dCurDisparity /= double(m_dSearchLevel);
#endif
//          m_acLabel2Depth[l] = BYTE( mp * ( dCurDisparity - dMinDisparity ) / ( rfl + (dCurDisparity+disparity_offset)*cvmGet(matEx_c2w, 2, 3) ) + 0.5 );
            m_acLabel2Depth[l] = (DepthType)( mp * ( dCurDisparity - dMinDisparity ) / ( rfl + (dCurDisparity+disparity_offset)*cvmGet(matEx_c2w, 2, 3) ) + 0.5 );
        }
    }
    else
    {
        for(l=0, iDisparity=iMinSearch; iDisparity<=iMaxSearch; l++, iDisparity++)
        {
            m_acLabel2Depth[l] = (DepthType) ((MAX_DEPTH-1)*((double)(iDisparity-iMinDisparity))/((double)(iMaxDisparity-iMinDisparity)) + 0.5);
        }
    }
    //Invert LUT for semi-automatic mode
#ifdef POZNAN_DYNAMIC_TABLE
    depth2label = (int*) malloc(MAX_DEPTH*sizeof(int));
#endif
    memset(depth2label, 0, MAX_DEPTH*sizeof(int));
    for(l=1; l<m_iNumOfLabels; l++) {
        BYTE a = m_acLabel2Depth[l-1];
        BYTE b = m_acLabel2Depth[l]-1;
        for(int l2 = a; l2<=b; l2++) {
          depth2label[l2] = l-1;
        }
    }
    depth2label[MAX_DEPTH-1] = depth2label[MAX_DEPTH-2]+1;
    return true;
}

#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
bool CEstimation::init_label2depth_z(int iNumberOfDepthSteps, double dNearestSearchDepthValue, double dFarthestSearchDepthValue, double dNearestDepthValue, double dFarthestDepthValue, unsigned int uiDepthType)
{
    int l, iDisparity;
    if((m_acLabel2Depth = (DepthType *)malloc(m_iNumOfLabels*sizeof(DepthType)))==NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
    memset(m_acLabel2Depth, 0, m_iNumOfLabels*sizeof(DepthType));

    if(uiDepthType)
    {
        //Fix Owieczka
        printf("Not supported DepthType\n");
        return false;
    }
    else
    {
        for(l=0; l < iNumberOfDepthSteps; l++)
        {
            m_acLabel2Depth[l] = (DepthType) ((MAX_DEPTH-1)*((double)(l))/((double)(iNumberOfDepthSteps)) + 0.5);
        }
    }
    //Invert LUT for semi-automatic mode
#ifdef POZNAN_DYNAMIC_TABLE
    depth2label = (int*) malloc(MAX_DEPTH*sizeof(int));
#endif
    memset(depth2label, 0, MAX_DEPTH*sizeof(int));
    for(l=1; l<m_iNumOfLabels; l++) {
        BYTE a = m_acLabel2Depth[l-1];
        BYTE b = m_acLabel2Depth[l]-1;
        for(int l2 = a; l2<=b; l2++) {
          depth2label[l2] = l-1;
        }
    }
    depth2label[MAX_DEPTH-1] = depth2label[MAX_DEPTH-2]+1;
    return true;
}
#endif

bool CEstimation::init_label2disparity_simple(double ratio[2], int iMinSearch, int iMaxSearch)
{
    int l, iDisparity;
    double dDisparity;

    if((m_aiLabel2Disparity[0] = (int *)malloc(2*m_iNumOfLabels*sizeof(int)))==NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
    m_aiLabel2Disparity[1] = &m_aiLabel2Disparity[0][m_iNumOfLabels];

    for(l=0, iDisparity=iMinSearch; iDisparity<=iMaxSearch; l++, iDisparity++)
    {
        dDisparity = double(iDisparity)*ratio[0];
#ifdef SUB_PEL_PRECISION
        dDisparity *= double(m_iPrecision)/double(m_dSearchLevel);
#endif
        m_aiLabel2Disparity[0][l] = (int)(dDisparity+0.5);

        dDisparity = double(iDisparity)*ratio[1];
#ifdef SUB_PEL_PRECISION
        dDisparity *= double(m_iPrecision)/double(m_dSearchLevel);
#endif
        m_aiLabel2Disparity[1][l] = (int)(dDisparity+0.5);
    }

    return true;
}

// Lookup ���̺� ���� 
// m_aiLabel2Disparity[left/right][��Ī������ disparity] = ���� disparity
//
bool CEstimation::init_label2disparity_from_z(double focal_length, double baseline[2], double disparity_offset[2], CvMat* matEx_c2w, unsigned int uiDepthType)
{
    int l, i;
    double d, temp1, temp2;
    double dDisparity;

    if((m_aiLabel2Disparity[0] = (int *)malloc(2*m_iNumOfLabels*sizeof(int)))==NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
    m_aiLabel2Disparity[1] = &m_aiLabel2Disparity[0][m_iNumOfLabels];

    temp1 = 1.0/m_dZnear - 1.0/m_dZfar;
    temp2 = 1.0/m_dZfar;
    for(l=0; l<m_iNumOfLabels; l++)
    {
        d = 1.0 / ( double(m_acLabel2Depth[l])*temp1/(MAX_DEPTH-1) + temp2 );
        if(uiDepthType) d = (d - cvmGet(matEx_c2w, 2, 3)) / cvmGet(matEx_c2w, 2, 2);

        for(i=0; i<2; i++)
        {
            dDisparity = focal_length * baseline[i] / d - disparity_offset[i];
#ifdef SUB_PEL_PRECISION
            dDisparity *= double(m_iPrecision);
#endif
            m_aiLabel2Disparity[i][l] = (int)(dDisparity+0.5);
        }
    }

    return true;
}

/*
bool CEstimation::init_homography2(CvMat *matIn[3], CvMat *matEx_c2w[3], unsigned int uiDepthType)
{
    int i, j, k, l;
    double val, z;
    CvMat *matProj_w2i[2];
    CvMat *matEx_w2c;
    CvMat *matInvIn_from, *mat_Rc2w_InvIn_from;
    double temp1, temp2;
    CvMat *src_points, *dst_points[2], *image, *world;

    if((m_matH_V2L = (CvMat **)malloc(m_iNumOfLabels*2*sizeof(CvMat *)))==NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
    m_matH_V2R = &m_matH_V2L[m_iNumOfLabels];

    matEx_w2c = cvCreateMat(3, 4, CV_64F);
    for(i=0; i<2; i++)
    {
        matProj_w2i[i] = cvCreateMat(3, 4, CV_64F);
        convertCameraParam(matEx_w2c, matEx_c2w[i]);
        cvmMul(matIn[i], matEx_w2c, matProj_w2i[i]);  // Proj = inMat_c2i * exMat_w2c
    }
    cvReleaseMat(&matEx_w2c);

    matInvIn_from               = cvCreateMat(3, 3, CV_64F);
    mat_Rc2w_InvIn_from = cvCreateMat(3, 3, CV_64F);
    cvInvert(matIn[2], matInvIn_from, CV_LU);
    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            val=0.0;
            for(k=0; k<3; k++)
            {
                val += cvmGet(matEx_c2w[2], i, k)*cvmGet(matInvIn_from, k, j);
            }
            cvmSet(mat_Rc2w_InvIn_from, i, j, val);
        }
    }

    src_points = cvCreateMat(4, 2, CV_64F);
    dst_points[0] = cvCreateMat(4, 2, CV_64F);
    dst_points[1] = cvCreateMat(4, 2, CV_64F);
    image = cvCreateMat(3, 1, CV_64F);
    world = cvCreateMat(4, 1, CV_64F);

    temp1 = 1.0/m_dZnear - 1.0/m_dZfar;
    temp2 = 1.0/m_dZfar;

    for(l=0; l<m_iNumOfLabels; l++)
    {
        m_matH_V2L[l] = cvCreateMat(3, 3, CV_64F);
        m_matH_V2R[l] = cvCreateMat(3, 3, CV_64F);

        z = 1.0 / ( double(m_acLabel2Depth[l])*temp1/255.0 + temp2 );
        if(uiDepthType==0)
            z = cvmGet(matEx_c2w[2], 2, 2) * z + cvmGet(matEx_c2w[2], 2, 3);

        cvmSet(world, 2, 0, z);
        cvmSet(world, 3, 0, 1.0);

        cvmSet(image, 0, 0, 0.0); cvmSet(src_points, 0, 0, 0.0);
        cvmSet(image, 1, 0, 0.0); cvmSet(src_points, 0, 1, 0.0);
        cvmSet(image, 2, 0, 1.0);
        image2world_with_z(mat_Rc2w_InvIn_from, matEx_c2w[2], image, world);
        for(i=0; i<2; i++)
        {
            cvmMul(matProj_w2i[i], world, image);
            cvmSet(dst_points[i], 0, 0, cvmGet(image, 0, 0) / cvmGet(image, 2, 0));
            cvmSet(dst_points[i], 0, 1, cvmGet(image, 1, 0) / cvmGet(image, 2, 0));
        }

        cvmSet(image, 0, 0, 0.0);                   cvmSet(src_points, 1, 0, 0.0);
        cvmSet(image, 1, 0, double(m_iHeight-1));   cvmSet(src_points, 1, 1, double(m_iHeight-1));
        cvmSet(image, 2, 0, 1.0);
        image2world_with_z(mat_Rc2w_InvIn_from, matEx_c2w[2], image, world);
        for(i=0; i<2; i++)
        {
            cvmMul(matProj_w2i[i], world, image);
            cvmSet(dst_points[i], 1, 0, cvmGet(image, 0, 0) / cvmGet(image, 2, 0));
            cvmSet(dst_points[i], 1, 1, cvmGet(image, 1, 0) / cvmGet(image, 2, 0));
        }

        cvmSet(image, 0, 0, double(m_iWidth-1));    cvmSet(src_points, 2, 0, double(m_iWidth-1));
        cvmSet(image, 1, 0, double(m_iHeight-1));   cvmSet(src_points, 2, 1, double(m_iHeight-1));
        cvmSet(image, 2, 0, 1.0);
        image2world_with_z(mat_Rc2w_InvIn_from, matEx_c2w[2], image, world);
        for(i=0; i<2; i++)
        {
            cvmMul(matProj_w2i[i], world, image);
            cvmSet(dst_points[i], 2, 0, cvmGet(image, 0, 0) / cvmGet(image, 2, 0));
            cvmSet(dst_points[i], 2, 1, cvmGet(image, 1, 0) / cvmGet(image, 2, 0));
        }

        cvmSet(image, 0, 0, double(m_iWidth-1));    cvmSet(src_points, 3, 0, double(m_iWidth-1));
        cvmSet(image, 1, 0, 0.0);                   cvmSet(src_points, 3, 1, 0.0);
        cvmSet(image, 2, 0, 1.0);
        image2world_with_z(mat_Rc2w_InvIn_from, matEx_c2w[2], image, world);
        for(i=0; i<2; i++)
        {
            cvmMul(matProj_w2i[i], world, image);
            cvmSet(dst_points[i], 3, 0, cvmGet(image, 0, 0) / cvmGet(image, 2, 0));
            cvmSet(dst_points[i], 3, 1, cvmGet(image, 1, 0) / cvmGet(image, 2, 0));
        }

        cvFindHomography(src_points, dst_points[0], m_matH_V2L[l]);
        cvFindHomography(src_points, dst_points[1], m_matH_V2R[l]);
    }

    cvReleaseMat(&src_points);
    cvReleaseMat(&dst_points[0]);
    cvReleaseMat(&dst_points[1]);
    cvReleaseMat(&image);
    cvReleaseMat(&world);
    cvReleaseMat(&matInvIn_from);
    cvReleaseMat(&mat_Rc2w_InvIn_from);
    cvReleaseMat(&matProj_w2i[0]);
    cvReleaseMat(&matProj_w2i[1]);

    return true;
}
//*/


//Poznan - Krzysztof Wegner - Owieczka - Fix the homography for Epipolar line search
bool CEstimation::init_homography(CvMat *matIn[3], CvMat *matEx_c2w[3], unsigned int uiDepthType)
{
    int i, l;//j, k
    double z; //val
    CvMat *matProj_w2i[3];
    CvMat *matProj_w2i_x[3];
    CvMat *matH[2];
    CvMat *matEx_w2c;
    //CvMat *matInvIn_from, *mat_Rc2w_InvIn_from;
    CvMat *matInvProj;
    double temp1, temp2;
    CvMat *src_points, *dst_points[2], *image, *world;

    if((m_matH_V2L = (CvMat **)malloc(m_iNumOfLabels*2*sizeof(CvMat *)))==NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
        return false;
    }
    m_matH_V2R = &m_matH_V2L[m_iNumOfLabels];

    //z2*m2 = K2*(R2|T3)*i(R3|T3)*iK3*z3*m3

    matEx_w2c = cvCreateMat(3, 4, CV_64F);
    for(i=0; i<3; i++)
    {
        matProj_w2i[i] = cvCreateMat(3, 4, CV_64F);
        convertCameraParam(matEx_w2c, matEx_c2w[i]);
        cvmMul(matIn[i], matEx_w2c, matProj_w2i[i]);  // Proj = inMat_c2i * exMat_w2c

        matProj_w2i_x[i] = cvCreateMat(4, 4, CV_64F);
        for(int y=0;y<3;y++)
        for(int x=0;x<4;x++)
          cvmSet(matProj_w2i_x[i],y,x,cvmGet(matProj_w2i[i],y,x));
        for(int x=0;x<3;x++)
          cvmSet(matProj_w2i_x[i],3,x,0);
        cvmSet(matProj_w2i_x[i],3,3,1.0);
    }
    cvReleaseMat(&matEx_w2c);

    matInvProj = cvCreateMat(4,4,CV_64F);

    cvmInvert(matProj_w2i_x[2],matInvProj);
    for(i=0;i<2;i++)
    {
      matH[i] = cvCreateMat(3, 4, CV_64F);
      cvmMul(matProj_w2i[i],matInvProj,matH[i]);
    }

    cvReleaseMat(&matInvProj);

    src_points = cvCreateMat(4, 2, CV_64F);
    dst_points[0] = cvCreateMat(4, 2, CV_64F);
    dst_points[1] = cvCreateMat(4, 2, CV_64F);
    image = cvCreateMat(3, 1, CV_64F);
    world = cvCreateMat(4, 1, CV_64F);

    temp1 = 1.0/m_dZnear - 1.0/m_dZfar;
    temp2 = 1.0/m_dZfar;

    for(l=0; l<m_iNumOfLabels; l++)
    {
        m_matH_V2L[l] = cvCreateMat(3, 3, CV_64F);
        m_matH_V2R[l] = cvCreateMat(3, 3, CV_64F);
        z = 1.0 / ( double(m_acLabel2Depth[l])*temp1/(MAX_DEPTH-1) + temp2 );
        if(uiDepthType==0)
            z = cvmGet(matEx_c2w[2], 2, 2) * z + cvmGet(matEx_c2w[2], 2, 3);
        
        cvmSet(world, 0, 0, 0.0);
        cvmSet(world, 1, 0, 0.0);
        cvmSet(world, 2, 0, 1.0);
        cvmSet(world, 3, 0, 1.0/z);

        cvmSet(image, 0, 0, 0.0); cvmSet(src_points, 0, 0, 0.0);
        cvmSet(image, 1, 0, 0.0); cvmSet(src_points, 0, 1, 0.0);
        cvmSet(image, 2, 0, 1.0);
        //image2world_with_z(mat_Rc2w_InvIn_from, matEx_c2w[2], image, world);
        
        for(i=0; i<2; i++)
        {
            cvmMul(matH[i],world,image);
            //cvmMul(matProj_w2i[i], world, image);
            cvmSet(dst_points[i], 0, 0, cvmGet(image, 0, 0) / cvmGet(image, 2, 0));
            cvmSet(dst_points[i], 0, 1, cvmGet(image, 1, 0) / cvmGet(image, 2, 0));
        }

        cvmSet(world, 0, 0, 0.0);
        cvmSet(world, 1, 0, double(m_iHeight-1));
        cvmSet(world, 2, 0, 1.0);
        cvmSet(world, 3, 0, 1.0/z);
        
        cvmSet(image, 0, 0, 0.0);                   cvmSet(src_points, 1, 0, 0.0);
        cvmSet(image, 1, 0, double(m_iHeight-1));   cvmSet(src_points, 1, 1, double(m_iHeight-1));
        cvmSet(image, 2, 0, 1.0);
        //image2world_with_z(mat_Rc2w_InvIn_from, matEx_c2w[2], image, world);
        for(i=0; i<2; i++)
        {
            cvmMul(matH[i],world,image);
            //cvmMul(matProj_w2i[i], world, image);
            cvmSet(dst_points[i], 1, 0, cvmGet(image, 0, 0) / cvmGet(image, 2, 0));
            cvmSet(dst_points[i], 1, 1, cvmGet(image, 1, 0) / cvmGet(image, 2, 0));
        }

        cvmSet(world, 0, 0, double(m_iWidth-1));
        cvmSet(world, 1, 0, double(m_iHeight-1));
        cvmSet(world, 2, 0, 1.0);
        cvmSet(world, 3, 0, 1.0/z);

        cvmSet(image, 0, 0, double(m_iWidth-1));    cvmSet(src_points, 2, 0, double(m_iWidth-1));
        cvmSet(image, 1, 0, double(m_iHeight-1));   cvmSet(src_points, 2, 1, double(m_iHeight-1));
        cvmSet(image, 2, 0, 1.0);
        //image2world_with_z(mat_Rc2w_InvIn_from, matEx_c2w[2], image, world);
        for(i=0; i<2; i++)
        {
            cvmMul(matH[i],world,image);
            //cvmMul(matProj_w2i[i], world, image);
            cvmSet(dst_points[i], 2, 0, cvmGet(image, 0, 0) / cvmGet(image, 2, 0));
            cvmSet(dst_points[i], 2, 1, cvmGet(image, 1, 0) / cvmGet(image, 2, 0));
        }

        cvmSet(world, 0, 0, double(m_iWidth-1));
        cvmSet(world, 1, 0, 0.0);
        cvmSet(world, 2, 0, 1.0);
        cvmSet(world, 3, 0, 1.0/z);

        cvmSet(image, 0, 0, double(m_iWidth-1));    cvmSet(src_points, 3, 0, double(m_iWidth-1));
        cvmSet(image, 1, 0, 0.0);                   cvmSet(src_points, 3, 1, 0.0);
        cvmSet(image, 2, 0, 1.0);
        //image2world_with_z(mat_Rc2w_InvIn_from, matEx_c2w[2], image, world);
        for(i=0; i<2; i++)
        {
            cvmMul(matH[i],world,image);
            //cvmMul(matProj_w2i[i], world, image);
            cvmSet(dst_points[i], 3, 0, cvmGet(image, 0, 0) / cvmGet(image, 2, 0));
            cvmSet(dst_points[i], 3, 1, cvmGet(image, 1, 0) / cvmGet(image, 2, 0));
        }

        cvFindHomography(src_points, dst_points[0], m_matH_V2L[l]);
        cvFindHomography(src_points, dst_points[1], m_matH_V2R[l]);
    }

    cvReleaseMat(&src_points);
    cvReleaseMat(&dst_points[0]);
    cvReleaseMat(&dst_points[1]);
    cvReleaseMat(&image);
    cvReleaseMat(&world);
    //cvReleaseMat(&matInvIn_from);
    //cvReleaseMat(&mat_Rc2w_InvIn_from);
    cvReleaseMat(&matProj_w2i[0]);
    cvReleaseMat(&matProj_w2i[1]);

    return true;
}


//
// Wrapping/Agent �Լ�. ���� �Լ��� function pointer�� �����.
//
void CEstimation::block_matching(CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter_prev, bool temporalonoff, double threshold)
{
	// �� memset�� �Ⱦ���?
    for (int i=0; i<m_iNumOfLabels; i++)
        for (int j=0; j<m_iHeight*m_iWidth; j++){
#ifdef POZNAN_TWO_COSTS
            errors_left[i][j] = 0;
            errors_right[i][j] = 0;
#else
            errors[i][j] = 0;
#endif
        }

    if (gc_mode != SEMI_D_REFRESH) {
        printf("\nTemporal Enhancement: %s\n", temporalonoff?"ON":"OFF");
        (this->*func_block_matching)(yuvLeft, yuvRight, yuvCenter, srcSEGM);
    }
}

bool CEstimation::readCameraParam(FILE *fp, const char *target_id, CvMat *inMat, CvMat *exMat)
{
    char id[64];
    int read=0;
    double gomi[2];

    if(fseek(fp, 0, SEEK_SET)) return false;

    while( fscanf(fp, "%s", id) == 1 )
    {
        if( strcmp(target_id, id) == 0 )
        {
            read += fscanf(fp, "%lf %lf %lf", &inMat->data.db[0], &inMat->data.db[1], &inMat->data.db[2]);
            read += fscanf(fp, "%lf %lf %lf", &inMat->data.db[3], &inMat->data.db[4], &inMat->data.db[5]);
            read += fscanf(fp, "%lf %lf %lf", &inMat->data.db[6], &inMat->data.db[7], &inMat->data.db[8]);
            read += fscanf(fp, "%lf %lf", &gomi[0], &gomi[1]);
            read += fscanf(fp, "%lf %lf %lf %lf", &exMat->data.db[0], &exMat->data.db[1], &exMat->data.db[ 2], &exMat->data.db[ 3]);
            read += fscanf(fp, "%lf %lf %lf %lf", &exMat->data.db[4], &exMat->data.db[5], &exMat->data.db[ 6], &exMat->data.db[ 7]);
            read += fscanf(fp, "%lf %lf %lf %lf", &exMat->data.db[8], &exMat->data.db[9], &exMat->data.db[10], &exMat->data.db[11]);
            if(read!=23) return false;
            return true;
        }
    }

    return false;
}

void CEstimation::convertCameraParam(CvMat *exMat_dst, CvMat *exMat_src)
{
    int i, j;
    double val;

    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            //cvmSet(exMat_dst, i, j, cvmGet(exMat_src, j, i));
          cvmSet(exMat_dst, i, j, cvmGet(exMat_src, i, j));
        }
    }

    for(i=0; i<3; i++)
    {
        val = 0.0;
        for(j=0; j<3; j++)
        {
            val -= cvmGet(exMat_dst, i, j)*cvmGet(exMat_src, j, 3);
        }
        cvmSet(exMat_dst, i, 3, val);
    }
}

void CEstimation::image2world_with_z(CvMat *mat_Rc2w_invIN_from, CvMat *matEX_c2w_from, CvMat *image, CvMat *world)
{
    CvMat *temp;
    double s;

    temp = cvCreateMat(3, 1, CV_64F);

    cvmMul(mat_Rc2w_invIN_from, image, temp);

    s = ( cvmGet(world, 2, 0) - cvmGet(matEX_c2w_from, 2, 3) ) / cvmGet(temp, 2, 0);

    cvmSet(world, 0, 0, s*cvmGet(temp, 0, 0) + cvmGet(matEX_c2w_from, 0, 3));
    cvmSet(world, 1, 0, s*cvmGet(temp, 1, 0) + cvmGet(matEX_c2w_from, 1, 3));

    cvReleaseMat(&temp);
}

//Poznan start
//Soft-Segmentation by Owieczka - Krzysztof Wegner - PUT Poznan
void CEstimation::block_matching_disparity_soft(CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM)
{
    int i, j, d, pp;
    int target_pixel_u;
    CostType error_l, error_r;
    double*** w_l1;
    double*** w_r1;
    double*** w_c1;
    double* exptab;
    int wx = m_uiSoftBlockWidth/2;
    int wy = m_uiSoftBlockHeight/2;
    double gc = m_dSoftColorCoeff;
    double gd = m_dSoftDistanceCoeff;
    int a,b,aa;
    int ii,jj,kk;
    double right_denominator;
    double left_denominator;
    double right_nominator;
    double left_nominator;
    double dd,ddd;

    printf("Soft Seg\n");

    //Initiation soft-segment weights
    exptab = (double*) calloc(256*256*3,sizeof(double));
    for(i=0;i<256*256*3;i++)
        exptab[i] = exp(-sqrt((double)i)/gc);

    w_l1 = (double***) calloc(m_iWidth*m_iPrecision,sizeof(double**));
    for(i=0;i<m_iWidth*m_iPrecision;i++)
    {
        w_l1[i] = (double**) calloc(2*wy+1,sizeof(double*));
        for(j=0;j<wy*2+1;j++)
            w_l1[i][j] = (double*) calloc((2*wx+1),sizeof(double));
    }//*/

    //Initiation soft-segment weights
    w_r1 = (double***) calloc(m_iWidth*m_iPrecision,sizeof(double**));
    for(i=0;i<m_iWidth*m_iPrecision;i++)
    {
        w_r1[i] = (double**) calloc(2*wy+1,sizeof(double*));
        for(j=0;j<wy*2+1;j++)
            w_r1[i][j] = (double*) calloc((2*wx+1),sizeof(double));
    }

    w_c1 = (double***) calloc(m_iWidth,sizeof(double**));
    for(i=0;i<m_iWidth;i++)
    {
        w_c1[i] = (double**) calloc(2*wy+1,sizeof(double*));
        for(j=0;j<wy*2+1;j++)
            w_c1[i][j] = (double*) calloc((2*wx+1),sizeof(double));
    }


    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++)
        {
            //Compute weight of soft segment for every pixel
            for(jj=0;jj<=wy*2;jj++)
            for(ii=0;ii<=wx*2;ii++)
            {
                b = min(max(j+jj-wy,0),m_iHeight-1);
                a = min(max(i+ii-wx,0),m_iWidth-1);

                dd = exp(-sqrt((double)(i-a)*(i-a)+(double)(j-b)*(j-b))/gd);

                for(kk=0;kk<m_iPrecision;kk++)
                {
                 aa = min(max((i+ii-wx)*m_iPrecision+kk,0),m_iWidth*m_iPrecision-1);
#ifdef POZNAN_TWOVIEW_SUPPORT
                 if(yuvRight)
#endif
                 w_r1[i*m_iPrecision+kk][jj/*+wy*/][ii/*+wx*/] = exptab[(yuvRight->Y[j][i*m_iPrecision+kk]-yuvRight->Y[b][aa])*(yuvRight->Y[j][i*m_iPrecision+kk]-yuvRight->Y[b][aa])
                                                       /*+(yuvRight->U[j][i]-yuvRight->U[b][a])*(yuvRight->U[j][i]-yuvRight->U[b][a])
                                                       +(yuvRight->V[j][i]-yuvRight->V[b][a])*(yuvRight->V[j][i]-yuvRight->V[b][a])*/];
#ifdef POZNAN_TWOVIEW_SUPPORT
                 if(yuvLeft)
#endif
                 w_l1[i*m_iPrecision+kk][jj/*+wy*/][ii/*+wx*/] = exptab[(yuvLeft->Y[j][i*m_iPrecision+kk]-yuvLeft->Y[b][aa])*(yuvLeft->Y[j][i*m_iPrecision+kk]-yuvLeft->Y[b][aa])
                                                       /*+(yuvLeft->U[j][i]-yuvRight->U[b][a])*(yuvRight->U[j][i]-yuvRight->U[b][a])
                                                        +(yuvLeft->V[j][i]-yuvRight->V[b][a])*(yuvRight->V[j][i]-yuvRight->V[b][a])*/];
                 ddd = exp(-sqrt((double)(i-aa/(double)m_iPrecision)*(i-aa/(double)m_iPrecision)+(double)(j-b)*(j-b))/gd);
                 w_r1[i*m_iPrecision+kk][jj/*+wy*/][ii/*+wx*/] *= ddd;
                 w_l1[i*m_iPrecision+kk][jj/*+wy*/][ii/*+wx*/] *= ddd;
                }
                w_c1[i][jj/*+wy*/][ii/*+wx*/] = exptab[(yuvCenter->Y[j][i]-yuvCenter->Y[b][a])*(yuvCenter->Y[j][i]-yuvCenter->Y[b][a])
                                                      /*+(yuvCenter->U[j][i]-yuvCenter->U[b][a])*(yuvCenter->U[j][i]-yuvCenter->U[b][a])
                                                      +(yuvCenter->V[j][i]-yuvCenter->V[b][a])*(yuvCenter->V[j][i]-yuvCenter->V[b][a])*/];


                w_c1[i][jj/*+wy*/][ii/*+wx*/] *= dd;
            }
        }
        for(i=0; i<m_iWidth; i++, pp++)
        {
            for(d=0; d<m_iNumOfLabels; d++)
            {
                right_nominator = 0.0;
                right_denominator = 0.0;
                left_nominator = 0.0;
                left_denominator = 0.0;
                //Eveluate of similarity metric - compute soft-segment SAD
                for(jj=0;jj<=wy*2;jj++)
                for(ii=0;ii<=wx*2;ii++)
                {
#ifdef POZNAN_TWOVIEW_SUPPORT
                    if(yuvLeft)
                    {
#endif
#ifdef SUB_PEL_PRECISION
                    target_pixel_u = (i+ii-wx)*m_iPrecision + m_aiLabel2Disparity[0][d];
#else
                    target_pixel_u = (i+ii-wx) + m_aiLabel2Disparity[0][d];
#endif
                    if(target_pixel_u>=0 && target_pixel_u<m_iMaxWidth && j+jj-wy>=0 && j+jj-wy<m_iHeight && i+ii-wx>=0 && i+ii-wx<m_iWidth)
                    {
                        error_l = byte_abs[yuvLeft->Y[j+jj-wy][target_pixel_u] - yuvCenter->Y[j+jj-wy][i+ii-wx]];
#ifdef SUB_PEL_PRECISION
                        target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d]; //Niepotrzebne liczone n razy
#else
                        target_pixel_u = i + m_aiLabel2Disparity[0][d];
#endif
                        if(target_pixel_u>=0 && target_pixel_u<m_iMaxWidth)
                        {
                            left_nominator += w_l1[target_pixel_u][jj][ii]*w_c1[i][jj][ii]*error_l;
                            left_denominator +=w_l1[target_pixel_u][jj][ii]*w_c1[i][jj][ii];
                        }

                    }
                    else
                        error_l = 255;

#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                    if(yuvRight)
                    {
#endif
#ifdef SUB_PEL_PRECISION
                    target_pixel_u = (i+ii-wx)*m_iPrecision - m_aiLabel2Disparity[1][d];
#else
                    target_pixel_u = (i+ii-wx) - m_aiLabel2Disparity[1][d];
#endif
                    if(target_pixel_u>=0 && target_pixel_u<m_iMaxWidth && j+jj-wy>=0 && j+jj-wy<m_iHeight && i+ii-wx>=0 && i+ii-wx<m_iWidth)
                    {
                        error_r = byte_abs[yuvRight->Y[j+jj-wy][target_pixel_u] - yuvCenter->Y[j+jj-wy][i+ii-wx]];
#ifdef SUB_PEL_PRECISION
                        target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
#else
                        target_pixel_u = i - m_aiLabel2Disparity[1][d];
#endif
                        if(target_pixel_u>=0 && target_pixel_u<m_iMaxWidth)
                        {
                            right_nominator += w_r1[target_pixel_u][jj][ii]*w_c1[i][jj][ii]*error_r;
                            right_denominator +=w_r1[target_pixel_u][jj][ii]*w_c1[i][jj][ii];
                        }
                    }
                    else
                        error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
#endif
                }
                //printf("%f/%f %f/%f\n",left_nominator,left_denominator,right_nominator,right_denominator);
                if(left_denominator>0)
                    error_l = left_nominator/left_denominator;
                else
                    error_l = 255;
                if(right_denominator>0)
                    error_r = right_nominator/right_denominator;
                else
                    error_r = 255;
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
                //errors[d][pp] = error_r;
            }//d
        }//i
        //printf(".");
		printf("\r%4d/%4d", j, m_iHeight);
    }//j
    printf("END\n");


    //Clear memory
    for(i=0;i<m_iWidth*m_iPrecision;i++)
    {
        for(j=0;j<wy*2+1;j++)
        {
            free(w_l1[i][j]);
            free(w_r1[i][j]);
        }
        free(w_l1[i]);
        free(w_r1[i]);
    }
    for(i=0;i<m_iWidth;i++)
    {
        for(j=0;j<wy*2+1;j++)
        {
            free(w_c1[i][j]);
        }
        free(w_c1[i]);
    }
    free(w_l1);
    free(w_r1);
    free(w_c1);

    free(exptab);

    return;
}
// Poznan end

void CEstimation::block_matching_homography_pixel(CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM)
{
    int i, j, d, pp;
    CvMat* pix[2];
    int target_pixel_u, target_pixel_v;
    CostType error_l, error_r;


    pix[0] = cvCreateMat(3, 1, CV_64F);
    pix[1] = cvCreateMat(3, 1, CV_64F);
    cvmSet(pix[0], 2, 0, 1.0);

    for(j=pp=0; j<m_iHeight; j++)
    {
        cvmSet(pix[0], 1, 0, j);
        for(i=0; i<m_iWidth; i++, pp++)
        {
            cvmSet(pix[0], 0, 0, i);
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight)
                {
                    error_l = byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]];
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif


                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight)
                {
                    error_r = byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]];
                }
                else
                {
                    error_r = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif


            }//d
        }//i
        //printf(".");
		printf("\r%4d/%4d", j, m_iHeight);
    }//j


    printf("END\n");
    fflush(stdout);

    cvReleaseMat(&pix[0]);
    cvReleaseMat(&pix[1]);

    return;
}

void CEstimation::block_matching_homography_window(CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM)
{
    int i, j, d, pp;
    CvMat* pix[2];
    int target_pixel_u, target_pixel_v;
    CostType error_l, error_r;
    int ci, cj, cu, cv, i1, im1, j1, jm1;
    int target_pixel_u_um1, target_pixel_v_um1;
    int target_pixel_u_u1, target_pixel_v_u1;
    int target_pixel_u_v1, target_pixel_v_v1;
    int target_pixel_u_vm1, target_pixel_v_vm1;
    int target_pixel_u_um1v1, target_pixel_v_um1v1;
    int target_pixel_u_um1vm1, target_pixel_v_um1vm1;
    int target_pixel_u_u1v1, target_pixel_v_u1v1;
    int target_pixel_u_u1vm1, target_pixel_v_u1vm1;

    pix[0] = cvCreateMat(3, 1, CV_64F);
    pix[1] = cvCreateMat(3, 1, CV_64F);
    cvmSet(pix[0], 2, 0, 1.0);

    CvMat* pix_um1[2];
    CvMat* pix_u1[2];
    CvMat* pix_vm1[2];
    CvMat* pix_v1[2];
    CvMat* pix_um1vm1[2];
    CvMat* pix_um1v1[2];
    CvMat* pix_u1vm1[2];
    CvMat* pix_u1v1[2];

    pix_um1[0] = cvCreateMat(3, 1, CV_64F);
    pix_um1[1] = cvCreateMat(3, 1, CV_64F);
    pix_u1[0] = cvCreateMat(3, 1, CV_64F);
    pix_u1[1] = cvCreateMat(3, 1, CV_64F);
    pix_vm1[0] = cvCreateMat(3, 1, CV_64F);
    pix_vm1[1] = cvCreateMat(3, 1, CV_64F);
    pix_v1[0] = cvCreateMat(3, 1, CV_64F);
    pix_v1[1] = cvCreateMat(3, 1, CV_64F);
    pix_um1v1[0] = cvCreateMat(3, 1, CV_64F);
    pix_um1v1[1] = cvCreateMat(3, 1, CV_64F);
    pix_um1vm1[0] = cvCreateMat(3, 1, CV_64F);
    pix_um1vm1[1] = cvCreateMat(3, 1, CV_64F);
    pix_u1vm1[0] = cvCreateMat(3, 1, CV_64F);
    pix_u1vm1[1] = cvCreateMat(3, 1, CV_64F);
    pix_u1v1[0] = cvCreateMat(3, 1, CV_64F);
    pix_u1v1[1] = cvCreateMat(3, 1, CV_64F);

    cvmSet(pix_um1[0],  2, 0, 1.0);
    cvmSet(pix_u1[0],  2, 0, 1.0);
    cvmSet(pix_vm1[0],  2, 0, 1.0);
    cvmSet(pix_v1[0],  2, 0, 1.0);
    cvmSet(pix_um1v1[0], 2, 0, 1.0);
    cvmSet(pix_um1vm1[0], 2, 0, 1.0);
    cvmSet(pix_u1v1[0], 2, 0, 1.0);
    cvmSet(pix_u1vm1[0], 2, 0, 1.0);

////////////////////////////////////////////////////////////////////
    pp=0;
    //j=0
    j=0, j1=1;
        //i=0
        i=0, ci=0, i1=1;
        cvmSet(pix[0],   1, 0, j);
        cvmSet(pix_u1[0], 1, 0, j);
        cvmSet(pix_v1[0],  1, 0, j1);
        cvmSet(pix_u1v1[0], 1, 0, j1);

        cvmSet(pix[0],   0, 0, i);
        cvmSet(pix_v1[0], 0, 0, i);
        cvmSet(pix_u1[0],  0, 0, i1);
        cvmSet(pix_u1v1[0], 0, 0, i1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2L[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2L[d], pix_u1v1[0], pix_u1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#endif
                if(target_pixel_u >= 0 && target_pixel_u_v1 >= 0 && target_pixel_u_u1 < m_iMaxWidth  && target_pixel_u_u1v1 < m_iMaxWidth
                && target_pixel_v >= 0 && target_pixel_v_u1 >= 0 && target_pixel_v_v1 < m_iMaxHeight && target_pixel_v_u1v1 < m_iMaxHeight)
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.2 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.1 * (byte_abs[yuvLeft->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                    + byte_abs[yuvLeft->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                    + byte_abs[yuvLeft->Y[target_pixel_v_u1v1][target_pixel_u_u1v1] - yuvCenter->Y[j1][i1]])
                             + 0.25 * (byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif

                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2R[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2R[d], pix_u1v1[0], pix_u1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u_v1 >= 0 && target_pixel_u_u1 < m_iMaxWidth && target_pixel_u_u1v1 < m_iMaxWidth
                && target_pixel_v >= 0 && target_pixel_v_u1 >= 0 && target_pixel_v_v1 < m_iMaxHeight && target_pixel_v_u1v1 < m_iMaxHeight)
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.2 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.1 * (byte_abs[yuvRight->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                    + byte_abs[yuvRight->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                    + byte_abs[yuvRight->Y[target_pixel_v_u1v1][target_pixel_u_u1v1] - yuvCenter->Y[j1][i1]])
                             + 0.25 * (byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;

#endif

#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif

            }//d
            pp+=1;
        //i=0 end
////////////////////////////////////////////////////////////////////
        //0<i<m_iWidth-1
        for(i=1; i<m_iWidth-1; i++, pp++)
        {
            cvmSet(pix[0],   1, 0, j);
            cvmSet(pix_u1[0], 1, 0, j);
            cvmSet(pix_um1[0], 1, 0, j);
            cvmSet(pix_v1[0],  1, 0, j1);
            cvmSet(pix_u1v1[0], 1, 0, j1);
            cvmSet(pix_um1v1[0], 1, 0, j1);

            ci=i>>1, i1=i+1, im1 =i-1;
            cvmSet(pix[0],   0, 0, i);
            cvmSet(pix_v1[0], 0, 0, i);
            cvmSet(pix_u1[0],  0, 0, i1);
            cvmSet(pix_u1v1[0], 0, 0, i1);
            cvmSet(pix_um1[0],  0, 0, im1);
            cvmSet(pix_um1v1[0], 0, 0, im1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2L[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2L[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2L[d], pix_u1v1[0], pix_u1v1[1]);
                cvmMul(m_matH_V2L[d], pix_um1v1[0], pix_um1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u_um1 >= 0 && target_pixel_u_um1v1>=0 && target_pixel_u_v1 >=0 && target_pixel_u_u1 < m_iMaxWidth && target_pixel_u_u1v1 < m_iMaxWidth
                && target_pixel_v>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_u1>=0 && target_pixel_v_um1>=0 && target_pixel_v_u1v1<m_iMaxHeight && target_pixel_v_um1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.15 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.07 * (byte_abs[yuvLeft->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_um1v1][target_pixel_u_um1v1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1v1][target_pixel_u_u1v1] - yuvCenter->Y[j1][i1]])
                             + 0.25 * (byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2R[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2R[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2R[d], pix_u1v1[0], pix_u1v1[1]);
                cvmMul(m_matH_V2R[d], pix_um1v1[0], pix_um1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u_um1 >= 0 && target_pixel_u_um1v1>=0 && target_pixel_u_v1 >=0 && target_pixel_u_u1 < m_iMaxWidth && target_pixel_u_u1v1 < m_iMaxWidth
                && target_pixel_v>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_u1>=0 && target_pixel_v_um1>=0 && target_pixel_v_u1v1<m_iMaxHeight && target_pixel_v_um1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.15 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.07 * (byte_abs[yuvRight->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_um1v1][target_pixel_u_um1v1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_u1v1][target_pixel_u_u1v1] - yuvCenter->Y[j1][i1]])
                             + 0.25 * (byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;
#endif

#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        }//i
////////////////////////////////////////////////////////////////////
        //i=m_iWidth-1
        i=m_iWidth-1, im1=i-1;
        cvmSet(pix[0],   1, 0, j);
        cvmSet(pix_um1[0], 1, 0, j);
        cvmSet(pix_v1[0],  1, 0, j1);
        cvmSet(pix_um1v1[0], 1, 0, j1);

        cvmSet(pix[0],   0, 0, i);
        cvmSet(pix_v1[0], 0, 0, i);
        cvmSet(pix_um1[0],  0, 0, im1);
        cvmSet(pix_um1v1[0], 0, 0, im1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2L[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2L[d], pix_um1v1[0], pix_um1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u<m_iMaxWidth && target_pixel_u_v1<m_iMaxWidth && target_pixel_u_um1>=0 && target_pixel_u_um1v1>=0
                && target_pixel_v>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_um1>=0 && target_pixel_v_um1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.2 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.1 * (byte_abs[yuvLeft->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                    + byte_abs[yuvLeft->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                    + byte_abs[yuvLeft->Y[target_pixel_v_um1v1][target_pixel_u_um1v1] - yuvCenter->Y[j1][im1]])
                             + 0.25 * (byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2R[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2R[d], pix_um1v1[0], pix_um1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u<m_iMaxWidth && target_pixel_u_v1<m_iMaxWidth && target_pixel_u_um1>=0 && target_pixel_u_um1v1>=0
                && target_pixel_v>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_um1>=0 && target_pixel_v_um1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.2 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.1 * (byte_abs[yuvRight->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                    + byte_abs[yuvRight->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j+1][i]]
                                    + byte_abs[yuvRight->Y[target_pixel_v_um1v1][target_pixel_u_um1v1] - yuvCenter->Y[j1][im1]])
                             + 0.25 * (byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        pp += 1;
        //printf(".");
		printf("\r%4d/%4d", j, m_iHeight);
        fflush(stdout);
        //i=m_iWidth-1 end
    //j=0 end

////////////////////////////////////////////////////////////////////
    //0<j<m_iheight-1
    for(j=1; j<m_iHeight-1; j++)
    {
        j1=j+1, jm1=j-1;
        //i=0
        i=0, i1=i+1;
        cvmSet(pix[0],   1, 0, j);
        cvmSet(pix_u1[0], 1, 0, j);
        cvmSet(pix_vm1[0],  1, 0, jm1);
        cvmSet(pix_u1vm1[0], 1, 0, jm1);
        cvmSet(pix_v1[0],  1, 0, j1);
        cvmSet(pix_u1v1[0], 1, 0, j1);

        cvmSet(pix[0],   0, 0, i);
        cvmSet(pix_vm1[0], 0, 0, i);
        cvmSet(pix_v1[0], 0, 0, i);
        cvmSet(pix_u1[0],  0, 0, i1);
        cvmSet(pix_u1vm1[0], 0, 0, i1);
        cvmSet(pix_u1v1[0], 0, 0, i1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2L[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2L[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2L[d], pix_u1vm1[0], pix_u1vm1[1]);
                cvmMul(m_matH_V2L[d], pix_u1v1[0], pix_u1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#endif


                if(target_pixel_u >= 0 && target_pixel_u_vm1>=0 && target_pixel_u_v1>=0 && target_pixel_u_u1<m_iMaxWidth && target_pixel_u_u1vm1<m_iMaxWidth && target_pixel_u_u1v1<m_iMaxWidth
                && target_pixel_v >= 0 && target_pixel_v_vm1>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_u1>=0 && target_pixel_v_u1vm1>=0 && target_pixel_v_u1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.15 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.07 * (byte_abs[yuvLeft->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1vm1][target_pixel_u_u1vm1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1v1][target_pixel_u_u1v1] - yuvCenter->Y[j1][i1]])
                             + 0.25 * (byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2R[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2R[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2R[d], pix_u1vm1[0], pix_u1vm1[1]);
                cvmMul(m_matH_V2R[d], pix_u1v1[0], pix_u1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u_vm1>=0 && target_pixel_u_v1>=0 && target_pixel_u_u1<m_iMaxWidth && target_pixel_u_u1vm1<m_iMaxWidth && target_pixel_u_u1v1<m_iMaxWidth
                && target_pixel_v >= 0 && target_pixel_v_vm1>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_u1>=0 && target_pixel_v_u1vm1>=0 && target_pixel_v_u1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.15 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.07 * (byte_abs[yuvRight->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_u1vm1][target_pixel_u_u1vm1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_u1v1][target_pixel_u_u1v1] - yuvCenter->Y[j1][i1]])
                             + 0.25 * (byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
            pp+=1;
        //i=0 end
////////////////////////////////////////////////////////////////////
        //0<i<m_iWidth-1
        for(i=1; i<m_iWidth-1; i++, pp++)
        {
            i1=i+1, im1=i-1;
            cvmSet(pix[0],   1, 0, j);
            cvmSet(pix_um1[0], 1, 0, j);
            cvmSet(pix_u1[0], 1, 0, j);
            cvmSet(pix_vm1[0],  1, 0, jm1);
            cvmSet(pix_um1vm1[0], 1, 0, jm1);
            cvmSet(pix_u1vm1[0], 1, 0, jm1);
            cvmSet(pix_v1[0],  1, 0, j1);
            cvmSet(pix_um1v1[0], 1, 0, j1);
            cvmSet(pix_u1v1[0], 1, 0, j1);

            cvmSet(pix[0],   0, 0, i);
            cvmSet(pix_vm1[0], 0, 0, i);
            cvmSet(pix_v1[0], 0, 0, i);
            cvmSet(pix_u1[0],  0, 0, i1);
            cvmSet(pix_u1vm1[0], 0, 0, i1);
            cvmSet(pix_u1v1[0], 0, 0, i1);
            cvmSet(pix_um1[0],  0, 0, im1);
            cvmSet(pix_um1vm1[0], 0, 0, im1);
            cvmSet(pix_um1v1[0], 0, 0, im1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2L[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2L[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2L[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2L[d], pix_u1vm1[0], pix_u1vm1[1]);
                cvmMul(m_matH_V2L[d], pix_um1vm1[0], pix_um1vm1[1]);
                cvmMul(m_matH_V2L[d], pix_u1v1[0], pix_u1v1[1]);
                cvmMul(m_matH_V2L[d], pix_um1v1[0], pix_um1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1>=0 && target_pixel_u_v1>=0 && target_pixel_u_u1<m_iMaxWidth && target_pixel_u_um1>=0
                    && target_pixel_u_u1vm1<m_iMaxWidth && target_pixel_u_um1vm1>=0 && target_pixel_u_u1v1<m_iMaxWidth && target_pixel_u_um1v1>=0
                    && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight
                    && target_pixel_v_vm1>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_u1>=0 && target_pixel_v_um1>=0
                    && target_pixel_v_u1vm1>=0 && target_pixel_v_um1vm1>=0 && target_pixel_v_u1v1<m_iMaxHeight && target_pixel_v_um1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.1 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.05 * (byte_abs[yuvLeft->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_um1vm1][target_pixel_u_um1vm1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_um1v1][target_pixel_u_um1v1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1vm1][target_pixel_u_u1vm1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1v1][target_pixel_u_u1v1] - yuvCenter->Y[j1][i1]])
                             + 0.25 * (byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2R[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2R[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2R[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2R[d], pix_u1vm1[0], pix_u1vm1[1]);
                cvmMul(m_matH_V2R[d], pix_um1vm1[0], pix_um1vm1[1]);
                cvmMul(m_matH_V2R[d], pix_u1v1[0], pix_u1v1[1]);
                cvmMul(m_matH_V2R[d], pix_um1v1[0], pix_um1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1v1 = (int)(cvmGet(pix_u1v1[1], 0, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1v1 = (int)(cvmGet(pix_u1v1[1], 1, 0)/cvmGet(pix_u1v1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1>=0 && target_pixel_u_v1>=0 && target_pixel_u_u1<m_iMaxWidth && target_pixel_u_um1>=0
                    && target_pixel_u_u1vm1<m_iMaxWidth && target_pixel_u_um1vm1>=0 && target_pixel_u_u1v1<m_iMaxWidth && target_pixel_u_um1v1>=0
                    && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight
                    && target_pixel_v_vm1>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_u1>=0 && target_pixel_v_um1>=0
                    && target_pixel_v_u1vm1>=0 && target_pixel_v_um1vm1>=0 && target_pixel_v_u1v1<m_iMaxHeight && target_pixel_v_um1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.1 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.05 * (byte_abs[yuvRight->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_um1vm1][target_pixel_u_um1vm1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_um1v1][target_pixel_u_um1v1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_u1vm1][target_pixel_u_u1vm1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_u1v1][target_pixel_u_u1v1] - yuvCenter->Y[j1][i1]])
                             + 0.25 * (byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;

#endif
#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        }//i
    //0<i<m_iWidth-1 end

////////////////////////////////////////////////////////////////////
  //i=m_iWidth-1

        i=m_iWidth-1, im1=i-1;
        cvmSet(pix[0],   1, 0, j);
        cvmSet(pix_um1[0], 1, 0, j);
        cvmSet(pix_vm1[0],  1, 0, jm1);
        cvmSet(pix_um1vm1[0], 1, 0, jm1);
        cvmSet(pix_v1[0],  1, 0, j1);
        cvmSet(pix_um1v1[0], 1, 0, j1);

        cvmSet(pix[0],   0, 0, i);
        cvmSet(pix_vm1[0], 0, 0, i);
        cvmSet(pix_v1[0], 0, 0, i);
        cvmSet(pix_um1[0],  0, 0, im1);
        cvmSet(pix_um1vm1[0], 0, 0, im1);
        cvmSet(pix_um1v1[0], 0, 0, im1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2L[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2L[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2L[d], pix_um1vm1[0], pix_um1vm1[1]);
                cvmMul(m_matH_V2L[d], pix_um1v1[0], pix_um1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1<m_iMaxWidth && target_pixel_u_v1<m_iMaxWidth && target_pixel_u_um1>=0
                    && target_pixel_u_um1vm1>=0 && target_pixel_u_um1v1>=0
                    && target_pixel_v >= 0
                    && target_pixel_v_vm1>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_um1>=0
                    && target_pixel_v_um1vm1>=0 && target_pixel_v_um1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.15 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.07 * (byte_abs[yuvLeft->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_um1vm1][target_pixel_u_um1vm1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_um1v1][target_pixel_u_um1v1] - yuvCenter->Y[j1][im1]])
                             + 0.25 * (byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2R[d], pix_v1[0], pix_v1[1]);
                cvmMul(m_matH_V2R[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2R[d], pix_um1vm1[0], pix_um1vm1[1]);
                cvmMul(m_matH_V2R[d], pix_um1v1[0], pix_um1v1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)*m_iPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_v1 = (int)(cvmGet(pix_v1[1], 0, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1v1 = (int)(cvmGet(pix_um1v1[1], 0, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_v1 = (int)(cvmGet(pix_v1[1], 1, 0)/cvmGet(pix_v1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1v1 = (int)(cvmGet(pix_um1v1[1], 1, 0)/cvmGet(pix_um1v1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1<m_iMaxWidth && target_pixel_u_v1<m_iMaxWidth && target_pixel_u_um1>=0
                    && target_pixel_u_um1vm1>=0 && target_pixel_u_um1v1>=0
                    && target_pixel_v >= 0
                    && target_pixel_v_vm1>=0 && target_pixel_v_v1<m_iMaxHeight && target_pixel_v_um1>=0
                    && target_pixel_v_um1vm1>=0 && target_pixel_v_um1v1<m_iMaxHeight
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.15 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.07 * (byte_abs[yuvRight->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_v1][target_pixel_u_v1] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_um1vm1][target_pixel_u_um1vm1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvRight->Y[target_pixel_v_um1v1][target_pixel_u_um1v1] - yuvCenter->Y[j1][im1]])
                             + 0.25 * (byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;
#endif

#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        pp+=1;
        //i=m_iWidth-1 end
        //printf(".");
		printf("\r%4d/%4d", j, m_iHeight);
        fflush(stdout);
        } //j
////////////////////////////////////////////////////////////////////
    //j=m_iHeight-1
    j=m_iHeight-1, jm1=j-1;

    //i=0
    i=0, i1=i+1;
        cvmSet(pix[0],   1, 0, j);
        cvmSet(pix_u1[0], 1, 0, j);
        cvmSet(pix_vm1[0],  1, 0, jm1);
        cvmSet(pix_u1vm1[0], 1, 0, jm1);

        cvmSet(pix[0],   0, 0, i);
        cvmSet(pix_vm1[0], 0, 0, i);
        cvmSet(pix_u1[0],  0, 0, i+1);
        cvmSet(pix_u1vm1[0], 0, 0, i+1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2L[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2L[d], pix_u1vm1[0], pix_u1vm1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1>=0 && target_pixel_u_u1<m_iMaxWidth && target_pixel_u_u1vm1<m_iMaxWidth
                    && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight
                    && target_pixel_v_vm1>=0 && target_pixel_v_u1>=0 && target_pixel_v_u1vm1>=0
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.2 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.1 * (byte_abs[yuvLeft->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                    + byte_abs[yuvLeft->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                    + byte_abs[yuvLeft->Y[target_pixel_v_u1vm1][target_pixel_u_u1vm1] - yuvCenter->Y[jm1][i1]])
                             + 0.25*(byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                   + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2R[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2R[d], pix_u1vm1[0], pix_u1vm1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1>=0 && target_pixel_u_u1<m_iMaxWidth && target_pixel_u_u1vm1<m_iMaxWidth
                    && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight
                    && target_pixel_v_vm1>=0 && target_pixel_v_u1>=0 && target_pixel_v_u1vm1>=0
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.2 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.1 * (byte_abs[yuvRight->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                    + byte_abs[yuvRight->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                    + byte_abs[yuvRight->Y[target_pixel_v_u1vm1][target_pixel_u_u1vm1] - yuvCenter->Y[jm1][i1]])
                             + 0.25 * (byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
            pp+=1;
            //i=0 end
////////////////////////////////////////////////////////////////////
        //0<i<m_iWidth-1

        cvmSet(pix[0],   1, 0, j);
        cvmSet(pix_um1[0], 1, 0, j);
        cvmSet(pix_u1[0], 1, 0, j);
        cvmSet(pix_vm1[0],  1, 0, jm1);
        cvmSet(pix_um1vm1[0], 1, 0, jm1);
        cvmSet(pix_u1vm1[0], 1, 0, jm1);

        for(i=1; i<m_iWidth-1; i++, pp++)
        {
            i1=i+1, im1=i-1;
            cvmSet(pix[0],   0, 0, i);
            cvmSet(pix_vm1[0], 0, 0, i);
            cvmSet(pix_u1[0],  0, 0, i1);
            cvmSet(pix_u1vm1[0], 0, 0, i1);
            cvmSet(pix_um1[0],  0, 0, im1);
            cvmSet(pix_um1vm1[0], 0, 0, im1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2L[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2L[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2L[d], pix_u1vm1[0], pix_u1vm1[1]);
                cvmMul(m_matH_V2L[d], pix_um1vm1[0], pix_um1vm1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1>=0 && target_pixel_u_u1<m_iMaxWidth && target_pixel_u_um1>=0
                    && target_pixel_u_u1vm1<m_iMaxWidth && target_pixel_u_um1vm1>=0
                    && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight
                    && target_pixel_v_vm1>=0 && target_pixel_v_u1>=0 && target_pixel_v_um1>=0
                    && target_pixel_v_u1vm1>=0 && target_pixel_v_um1vm1>=0
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.15 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.07 * (byte_abs[yuvLeft->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_um1vm1][target_pixel_u_um1vm1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvLeft->Y[target_pixel_v_u1vm1][target_pixel_u_u1vm1] - yuvCenter->Y[jm1][i1]])
                             + 0.25 * (byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2R[d], pix_u1[0], pix_u1[1]);
                cvmMul(m_matH_V2R[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2R[d], pix_u1vm1[0], pix_u1vm1[1]);
                cvmMul(m_matH_V2R[d], pix_um1vm1[0], pix_um1vm1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_u1 = (int)(cvmGet(pix_u1[1], 0, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 0, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_u1 = (int)(cvmGet(pix_u1[1], 1, 0)/cvmGet(pix_u1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_u1vm1 = (int)(cvmGet(pix_u1vm1[1], 1, 0)/cvmGet(pix_u1vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1>=0 && target_pixel_u_u1<m_iMaxWidth && target_pixel_u_um1>=0
                    && target_pixel_u_u1vm1<m_iMaxWidth && target_pixel_u_um1vm1>=0
                    && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight
                    && target_pixel_v_vm1>=0 && target_pixel_v_u1>=0 && target_pixel_v_um1>=0
                    && target_pixel_v_u1vm1>=0 && target_pixel_v_um1vm1>=0
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.15 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.07 * (byte_abs[yuvRight->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                      + byte_abs[yuvRight->Y[target_pixel_v_u1][target_pixel_u_u1] - yuvCenter->Y[j][i1]]
                                      + byte_abs[yuvRight->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                      + byte_abs[yuvRight->Y[target_pixel_v_um1vm1][target_pixel_u_um1vm1] - yuvCenter->Y[jm1][im1]]
                                      + byte_abs[yuvRight->Y[target_pixel_v_u1vm1][target_pixel_u_u1vm1] - yuvCenter->Y[jm1][i1]])
                             + 0.25*(byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                                      + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        }//i
        //0<i<m_iWidth-1 end
////////////////////////////////////////////////////////////////////
    //i=m_iWidth-1
    i=m_iWidth-1, im1=i-1;

        cvmSet(pix[0],   1, 0, j);
        cvmSet(pix_um1[0], 1, 0, j);
        cvmSet(pix_vm1[0],  1, 0, jm1);
        cvmSet(pix_um1vm1[0], 1, 0, jm1);

        cvmSet(pix[0],   0, 0, i);
        cvmSet(pix_vm1[0], 0, 0, i);
        cvmSet(pix_um1[0],  0, 0, im1);
        cvmSet(pix_um1vm1[0], 0, 0, im1);

            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
                cvmMul(m_matH_V2L[d], pix[0], pix[1]);
#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2L[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2L[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2L[d], pix_um1vm1[0], pix_um1vm1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1>=0 && target_pixel_u_um1>=0 && target_pixel_u_um1vm1 >=0
                    && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight
                    && target_pixel_v_vm1>=0 && target_pixel_v_um1>=0 && target_pixel_v_um1vm1 >=0
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_l = 0.2 * byte_abs[yuvLeft->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.1 * (byte_abs[yuvLeft->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                                    + byte_abs[yuvLeft->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                                    + byte_abs[yuvLeft->Y[target_pixel_v_um1vm1][target_pixel_u_um1vm1] - yuvCenter->Y[jm1][im1]])
                             + 0.25*(byte_abs[yuvLeft->U[cv][cu] - yuvCenter->U[cj][ci]]
                                    + byte_abs[yuvLeft->V[cv][cu] - yuvCenter->V[cj][ci]]);
                }
                else
                {
                    error_l = 255;
                }
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
                cvmMul(m_matH_V2R[d], pix[0], pix[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif
#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix[1], 2, 0) + 0.5);
#else
                target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
#endif

                cvmMul(m_matH_V2R[d], pix_vm1[0], pix_vm1[1]);
                cvmMul(m_matH_V2R[d], pix_um1[0], pix_um1[1]);
                cvmMul(m_matH_V2R[d], pix_um1vm1[0], pix_um1vm1[1]);

#ifdef SUB_PEL_PRECISION
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)*m_iPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_u_vm1 = (int)(cvmGet(pix_vm1[1], 0, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_u_um1 = (int)(cvmGet(pix_um1[1], 0, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_u_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 0, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#endif

#ifdef SUB_PEL_VERTICAL_PRECISION
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)*m_iVerticalPrecision/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#else
                target_pixel_v_vm1 = (int)(cvmGet(pix_vm1[1], 1, 0)/cvmGet(pix_vm1[1], 2, 0) + 0.5);
                target_pixel_v_um1 = (int)(cvmGet(pix_um1[1], 1, 0)/cvmGet(pix_um1[1], 2, 0) + 0.5);
                target_pixel_v_um1vm1 = (int)(cvmGet(pix_um1vm1[1], 1, 0)/cvmGet(pix_um1vm1[1], 2, 0) + 0.5);
#endif

                if(target_pixel_u >= 0 && target_pixel_u < m_iMaxWidth
                    && target_pixel_u_vm1>=0 && target_pixel_u_um1>=0 && target_pixel_u_um1vm1 >=0
                    && target_pixel_v >= 0 && target_pixel_v < m_iMaxHeight
                    && target_pixel_v_vm1>=0 && target_pixel_v_um1>=0 && target_pixel_v_um1vm1 >=0
                    )
                {
                    ci = i >> 1;
                    cj = j >> 1;
                    cu = target_pixel_u >> 1;
                    cv = target_pixel_v >> 1;
                    error_r = 0.2 * byte_abs[yuvRight->Y[target_pixel_v][target_pixel_u] - yuvCenter->Y[j][i]]
                             + 0.1 * (byte_abs[yuvRight->Y[target_pixel_v_um1][target_pixel_u_um1] - yuvCenter->Y[j][im1]]
                             + byte_abs[yuvRight->Y[target_pixel_v_vm1][target_pixel_u_vm1] - yuvCenter->Y[jm1][i]]
                             + byte_abs[yuvRight->Y[target_pixel_v_um1vm1][target_pixel_u_um1vm1] - yuvCenter->Y[jm1][im1]])
                             + 0.25*(byte_abs[yuvRight->U[cv][cu] - yuvCenter->U[cj][ci]]
                             + byte_abs[yuvRight->V[cv][cu] - yuvCenter->V[cj][ci]]);
                    }
                    else
                    {
                        error_r = 255;
                    }
#ifdef POZNAN_TWOVIEW_SUPPORT
                    }
                    else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                    errors_left[d][pp] = error_l;
                    errors_right[d][pp] = error_r;
#else
                    errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
            pp+=1;
        //j=m_iHeight-1 end
////////////////////////////////////////////////////////////////////

    printf("END\n");

    cvReleaseMat(&pix[0]);
    cvReleaseMat(&pix[1]);

    cvReleaseMat (&pix_u1[0]);
    cvReleaseMat (&pix_u1[1]);
    cvReleaseMat (&pix_um1[0]);
    cvReleaseMat (&pix_um1[1]);
    cvReleaseMat (&pix_v1[0]);
    cvReleaseMat (&pix_v1[1]);
    cvReleaseMat (&pix_vm1[0]);
    cvReleaseMat (&pix_vm1[1]);
    cvReleaseMat (&pix_u1v1[0]);
    cvReleaseMat (&pix_u1v1[1]);
    cvReleaseMat (&pix_u1vm1[0]);
    cvReleaseMat (&pix_u1vm1[1]);
    cvReleaseMat (&pix_um1v1[0]);
    cvReleaseMat (&pix_um1v1[1]);
    cvReleaseMat (&pix_um1vm1[0]);
    cvReleaseMat (&pix_um1vm1[1]);

    return;
}

//
// Pixel �ϳ��� �� ��. BM �Լ� �߿� ���� ���� �Լ��̹Ƿ�, ���⼭���� �м� �ϴ°��� �ٶ���  
// IN: left, center, right
// Out: errors
//
void CEstimation::block_matching_disparity_pixel(CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM)
{
    int i, j, d, pp;
    int target_pixel_u;
    CostType error_l, error_r;

    printf("\nMathing Block Size = %d\n", m_iMatchingBlock);

    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++, pp++) // error[depth][position of pixel], GC�� ���Ͽ� �̷��� ����.  
        {
            for(d=0; d<m_iNumOfLabels; d++)  // disparity search �����׿��� �ϴ� ��� ���  
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d]; // ���е��� index to disparity maaping table 
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
#endif
                if(target_pixel_u>=0 && target_pixel_u<m_iMaxWidth)
                    error_l = byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]];  // ABS-DIFF ��� (���̺����)
                else
                    error_l = 255;  // MAX �� (�� ������ ���Ϸ�����)
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
#endif
                if(target_pixel_u>=0 && target_pixel_u<m_iMaxWidth)
                    error_r = byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]];
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif

#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        }//i
        //printf(".");
		printf("\r%4d/%4d", j, m_iHeight);
        fflush(stdout);
    }//j
    printf("END\n");

    return;
}


//
// Block Matching �˰����� ���� (conventional �� disparity based�� ��� ���)
// LEFT-CENTER-RIGHT YUV 420 �Է� 
// Segmentation ������ ������� ���� ...

void CEstimation::block_matching_disparity_window(CIYuv<ImageType> *yuvLeft, CIYuv<ImageType> *yuvRight, CIYuv<ImageType> *yuvCenter, BYTE ***srcSEGM)
{
    int i, j, d, pp, cj, ci, j1, i1, jm1, im1;
    int target_pixel_u, cu, u1, um1;
    CostType error_l, error_r;

//  memset(error[0], 0, picsize*num_labels);

// NICT 3x3 block matching start
  int BlockSize = 1;
  if ( BlockSize !=0 )
  {
    pp = 0;
    // j=0
    j = 0, cj = 0, j1= j + 1;
        // i=0
        i = 0, ci = 0, i1 = 1;
            for(d=0; d<m_iNumOfLabels; d++) // �� Labels?
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + 1;
#endif
                if(target_pixel_u>=0 && u1<m_iMaxWidth)  // ���ʰ� �� 
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.2 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.1 * ( byte_abs[yuvLeft->Y[j][u1] - yuvCenter->Y[j][i1]]
                                    + byte_abs[yuvLeft->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]]
                                    + byte_abs[yuvLeft->Y[j1][u1] - yuvCenter->Y[j1][i1]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x2 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + 1;
#endif
                if(target_pixel_u>=0 && u1<m_iMaxWidth)  // �����ʰ� �� 
                {
                    cu = target_pixel_u >>1;
                    error_r = 0.2 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.1 * ( byte_abs[yuvRight->Y[j][u1] - yuvCenter->Y[j][i1]]
                                    + byte_abs[yuvRight->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]]
                                    + byte_abs[yuvRight->Y[j1][u1] - yuvCenter->Y[j1][i1]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x2 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
    // i=0 end
    pp += 1;

    // 0<i<m_iWidth - 1
     for(i=1; i<(m_iWidth - 1); i++, pp++)
        {
            ci = i >>1, i1 = i + 1, im1 = i - 1;
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + m_iPrecision, um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + 1, um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.15 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.07 * ( byte_abs[yuvLeft->Y[j][um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvLeft->Y[j][u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvLeft->Y[j1][um1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvLeft->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvLeft->Y[j1][u1] - yuvCenter->Y[j1][i1]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x3 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + m_iPrecision, um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + 1, um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_r = 0.15 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.07 * ( byte_abs[yuvRight->Y[j][um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvRight->Y[j][u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvRight->Y[j1][um1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvRight->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvRight->Y[j1][u1] - yuvCenter->Y[j1][i1]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x3 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        }//i

        // i=m_iWidth - 1
        i = m_iWidth - 1, ci = i >>1, im1 = i - 1;
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && target_pixel_u<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.2 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.1 * ( byte_abs[yuvLeft->Y[j][um1] - yuvCenter->Y[j][im1]]
                                    + byte_abs[yuvLeft->Y[j1][um1] - yuvCenter->Y[j1][im1]]
                                    + byte_abs[yuvLeft->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] );  // 2x2 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && target_pixel_u<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1, um1 = target_pixel_u -1;
                    error_r = 0.2 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.1 * ( byte_abs[yuvRight->Y[j][um1] - yuvCenter->Y[j][im1]]
                                    + byte_abs[yuvRight->Y[j1][um1] - yuvCenter->Y[j1][im1]]
                                    + byte_abs[yuvRight->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x2 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        //printf(".");
		printf("\r%4d/%4d", j, m_iHeight);
        fflush(stdout);
        // i=m_iWidth - 1 end
        pp += 1;
    // j=0 end


    // �ٿ���� ���� ���� ����?
    // 0<j<M_iHeight - 1
    for(j=1; j<(m_iHeight - 1); j++)
    {
        // i=0
        cj = j >>1, j1 = j + 1, jm1 = j - 1;
        i = 0, ci = 0, i1 = 1;
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + 1;
#endif
                if(target_pixel_u>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.15 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.07 * ( byte_abs[yuvLeft->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvLeft->Y[jm1][u1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvLeft->Y[j][u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvLeft->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvLeft->Y[j1][u1] - yuvCenter->Y[j1][i1]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 3x2 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + 1;
#endif
                if(target_pixel_u>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_r = 0.15 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.07 * ( byte_abs[yuvRight->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvRight->Y[jm1][u1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvRight->Y[j][u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvRight->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvRight->Y[j1][u1] - yuvCenter->Y[j1][i1]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 3x2 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
    // i=0 end
    pp += 1;
    // 0<i<m_iWidth - 1
        for(i=1; i<(m_iWidth - 1); i++, pp++)
        {
            ci = i >>1, i1 = i + 1, im1 = i - 1;
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + m_iPrecision, um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + 1, um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.1 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.05 * ( byte_abs[yuvLeft->Y[jm1][um1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvLeft->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvLeft->Y[jm1][u1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvLeft->Y[j][um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvLeft->Y[j][u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvLeft->Y[j1][um1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvLeft->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvLeft->Y[j1][u1] - yuvCenter->Y[j1][i1]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 3x3 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + m_iPrecision, um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + 1, um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_r = 0.1 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.05 * ( byte_abs[yuvRight->Y[jm1][um1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvRight->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvRight->Y[jm1][u1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvRight->Y[j][um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvRight->Y[j][u1] - yuvCenter->Y[j][i1]]
                                     + byte_abs[yuvRight->Y[j1][um1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvRight->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]]
                                     + byte_abs[yuvRight->Y[j1][u1] - yuvCenter->Y[j1][i1]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 3x3 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        }//i

        // i=m_iWidth - 1
        i = m_iWidth - 1, ci = i >>1, im1 = i - 1;
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && target_pixel_u<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.15 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.07 * ( byte_abs[yuvLeft->Y[jm1][um1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvLeft->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvLeft->Y[j][um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvLeft->Y[j1][um1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvLeft->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 3x2 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && target_pixel_u<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_r = 0.15 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.07 * ( byte_abs[yuvRight->Y[jm1][um1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvRight->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvRight->Y[j][um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvRight->Y[j1][um1] - yuvCenter->Y[j1][im1]]
                                     + byte_abs[yuvRight->Y[j1][target_pixel_u] - yuvCenter->Y[j1][i]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 3x2 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        //i=m_iWidth - 1 end
        pp += 1;

		//printf(".");
        printf("\r%4d/%4d",j, m_iHeight);
        fflush(stdout);
    }//j

    // j=m_iHeight - 1
    j = m_iHeight - 1, cj = j >>1, jm1 = j - 1;
        // i=0
        i = 0, ci = 0, i1 = 1;
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + 1;
#endif
                if(target_pixel_u>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.2 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.1 * ( byte_abs[yuvLeft->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                    + byte_abs[yuvLeft->Y[jm1][u1] - yuvCenter->Y[jm1][i1]]
                                    + byte_abs[yuvLeft->Y[j][u1] - yuvCenter->Y[j][i1]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x2 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + 1;
#endif
                if(target_pixel_u>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_r = 0.2 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.1 * ( byte_abs[yuvRight->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                    + byte_abs[yuvRight->Y[jm1][u1] - yuvCenter->Y[jm1][i1]]
                                    + byte_abs[yuvRight->Y[j][u1] - yuvCenter->Y[j][i1]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x2 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        // i=0 end
        pp += 1;

        // 0<i<m_iWidth - 1
        for(i=1; i<(m_iWidth - 1); i++, pp++)
        {
            ci = i >>1, i1 = i + 1, im1 = i - 1;
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + m_iPrecision, um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                u1 = target_pixel_u + 1, um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.15 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.07 * ( byte_abs[yuvLeft->Y[jm1][um1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvLeft->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvLeft->Y[jm1][u1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvLeft->Y[j][um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvLeft->Y[j][u1] - yuvCenter->Y[j][i1]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x3 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + m_iPrecision, um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                u1 = target_pixel_u + 1, um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && u1<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_r = 0.15 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.07 * ( byte_abs[yuvRight->Y[jm1][um1] - yuvCenter->Y[jm1][im1]]
                                     + byte_abs[yuvRight->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                     + byte_abs[yuvRight->Y[jm1][u1] - yuvCenter->Y[jm1][i1]]
                                     + byte_abs[yuvRight->Y[j][um1] - yuvCenter->Y[j][im1]]
                                     + byte_abs[yuvRight->Y[j][u1] - yuvCenter->Y[j][i1]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x3 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        }//i

        // i=m_iWidth - 1
        i = m_iWidth - 1, ci = i >>1, im1 = i - 1;
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
                um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
                um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && target_pixel_u<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_l = 0.2 * byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.1 * ( byte_abs[yuvLeft->Y[jm1][um1] - yuvCenter->Y[jm1][im1]]
                                    + byte_abs[yuvLeft->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                    + byte_abs[yuvLeft->Y[j][um1] - yuvCenter->Y[j][im1]] )
                            + 0.25 * ( byte_abs[yuvLeft->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvLeft->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x2 block matching
                }
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
                um1 = target_pixel_u - m_iPrecision;
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
                um1 = target_pixel_u - 1;
#endif
                if(um1>=0 && target_pixel_u<m_iMaxWidth)
                {
                    cu = target_pixel_u >>1;
                    error_r = 0.2 * byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]]
                            + 0.1 * ( byte_abs[yuvRight->Y[jm1][um1] - yuvCenter->Y[jm1][im1]]
                                    + byte_abs[yuvRight->Y[jm1][target_pixel_u] - yuvCenter->Y[jm1][i]]
                                    + byte_abs[yuvRight->Y[j][um1] - yuvCenter->Y[j][im1]] )
                            + 0.25 * ( byte_abs[yuvRight->U[cj][cu] - yuvCenter->U[cj][ci]]
                                     + byte_abs[yuvRight->V[cj][cu] - yuvCenter->V[cj][ci]] ); // 2x2 block matching
                }
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        // i=m_iWidth - 1 end
        //printf(".");
		printf("\r%4d/%4d", j, m_iHeight);
        fflush(stdout);
    // j=m_iHeight - 1 end

  }
  else
  {
// NICT end

    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++, pp++)
        {
            for(d=0; d<m_iNumOfLabels; d++)
            {
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvLeft)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision + m_aiLabel2Disparity[0][d];
#else
                target_pixel_u = i + m_aiLabel2Disparity[0][d];
#endif
                if(target_pixel_u>=0 && target_pixel_u<m_iMaxWidth)
                    error_l = byte_abs[yuvLeft->Y[j][target_pixel_u] - yuvCenter->Y[j][i]];
                else
                    error_l = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_l = 255;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
                if(yuvRight)
                {
#endif
#ifdef SUB_PEL_PRECISION
                target_pixel_u = i*m_iPrecision - m_aiLabel2Disparity[1][d];
#else
                target_pixel_u = i - m_aiLabel2Disparity[1][d];
#endif
                if(target_pixel_u>=0 && target_pixel_u<m_iMaxWidth)
                    error_r = byte_abs[yuvRight->Y[j][target_pixel_u] - yuvCenter->Y[j][i]];
                else
                    error_r = 255;
#ifdef POZNAN_TWOVIEW_SUPPORT
                }
                else error_r = 255;
#endif
#ifdef POZNAN_TWO_COSTS
                errors_left[d][pp] = error_l;
                errors_right[d][pp] = error_r;
#else
                errors[d][pp] = error_l<error_r?error_l:error_r;
#endif
            }//d
        }//i
        //printf(".");
		printf("\r%4d/%4d", j, m_iHeight);
        fflush(stdout);
    }//j

// NICT start
  }
// NICT end

    printf("END\n");

    return;
}

bool CEstimation::getDistrotion(DepthType **pDepth, int d, int scale)
{
    int i, j, pp, temp;

    if(d<0 || d>=m_iNumOfLabels) return false;

    for(i=pp=0; i<m_iHeight; i++)
    {
        for(j=0; j<m_iWidth; j++, pp++)
        {
#ifdef POZNAN_TWO_COSTS
            temp = MINSEL(errors_left[d][pp],errors_right[d][pp])*scale;
#else
            temp = errors[d][pp]*scale;
#endif
            pDepth[i][j] = (DepthType) (temp<(MAX_DEPTH-1)?temp:(MAX_DEPTH-1));
        }
    }
    return true;
}

//Nagoya start
//-----------------------------------------------------------------------------
// OpenCV1.1 automatically converts colour to grayscale for CV_LOAD_IMAGE_GRAYSCALE
// but OpenCV1.0 crashes, so we do it manually.
IplImage* LoadImageGray(const char* fname)
{
    IplImage* ipl_tmp = NULL;
    IplImage* ipl_out = NULL;

    ipl_tmp = cvLoadImage(fname, CV_LOAD_IMAGE_UNCHANGED);

    if (ipl_tmp  != NULL) {
        if(ipl_tmp->nChannels == 1) {
          ipl_out = ipl_tmp;
        } else {
          ipl_out = cvCreateImage(cvSize(ipl_tmp->width, ipl_tmp->height), IPL_DEPTH_8U, 1);
          cvCvtColor(ipl_tmp, ipl_out, CV_BGR2GRAY);
          SAFE_RELEASE_IMAGE(ipl_tmp)
        }
        printf("Found %s\n", fname); // DEBUG
    }
    return (ipl_out);
}

//-----------------------------------------------------------------------------
//If Semi-automatic mode is selected, this function is called at the start
//of every frame. It checks if manually created input files are available.
// ...MDM000 = Manual Disparity Map or Manual Depth Map
// ...MEM000 = Manual Edge map
// ...MSM000 = Manual Static map
// where 000 indicates frame number 0
bool CEstimation::load_man_images(const char* m_cFileManual, int framenumber, bool isfirstframe)
{   char fname[255];
    IplImage* ipl_new;

    // Manual Static Map
    sprintf(fname, "%sMSM%03d.png", m_cFileManual, framenumber);
    ipl_new = LoadImageGray(fname);
    if (ipl_new != NULL) {
        SAFE_RELEASE_IMAGE(ipl_manstatic)
        ipl_manstatic = ipl_new;
    }

    if (!isfirstframe) gc_mode = SEMI_TEMPORAL;

    // Manual Disparity Map / Manual Depth Map
    sprintf(fname, "%sMDM%03d.png", m_cFileManual, framenumber);
    SAFE_RELEASE_IMAGE(ipl_mandisp)
    ipl_mandisp = LoadImageGray(fname);
    if (isfirstframe && (ipl_mandisp == NULL)){
        fprintf(stderr, "\n\n\n");
        fprintf(stderr, "WARNING: Semi-automatic mode is selected, but no refresh\n");
        fprintf(stderr, "         frame was found for the first frame.\n");
        fprintf(stderr, "(file %s was not found)\n\n",fname);
    }

    // Manual Edge Map
    sprintf(fname, "%sMEM%03d.png", m_cFileManual, framenumber);
    SAFE_RELEASE_IMAGE(ipl_manedge)
    ipl_manedge = LoadImageGray(fname);

    if       ((ipl_mandisp != NULL) && (ipl_manedge != NULL)) {
        gc_mode = SEMI_D_E_REFRESH;
    }else if ((ipl_mandisp != NULL) && (ipl_manedge == NULL)) {
        gc_mode = SEMI_D_REFRESH;
    }

    switch(gc_mode) {
        case SEMI_D_REFRESH:
          printf("\nframe %3d; Mode: SEMI_D_REFRESH\n", framenumber); break;
        case SEMI_D_E_REFRESH:
          printf("\nframe %3d; Mode: SEMI_D_E_REFRESH\n", framenumber); break;
        case SEMI_TEMPORAL:
          printf("\nframe %3d; Mode: SEMI_TEMPORAL\n", framenumber); break;
    }

    return(true);

} //CEstimation::load_man_images

//-----------------------------------------------------------------------------
void CEstimation::getMotionMap_block(IplImage* src1, IplImage* src2,IplImage* dest, int block, double th)
//GIST code
{
    int w = src1->width;
    int h = src1->height;
    IplImage *src1f = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);
    IplImage *src2f = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);

    // MERL-arican
    // For Bilateral smoothing to work, two seperate buffers for input and output must be used.
    // In addition,sigma1 and sigma2 values must be entered
    // Otherwise, the function returns a completely black image.
    cvSmooth(src1, src1f, CV_BILATERAL,0,0,70,3);
    cvSmooth(src2, src2f, CV_BILATERAL,0,0,70,3);

    // GIST start 201001
    for(int y=0; y<h; y+=block)
    {
        for(int x=0; x<w; x+=block)
        {
            double mad = 0;
            int count = 0;

            for(int j=0; j<block; j++)
            {
                for(int i=0; i<block; i++)
                {
                    if(x+i>=0 && x+i<w && y+j>=0 && y+j<h)
                    {
                    mad += abs((BYTE)src1f->imageData[(y+j)*w+(x+i)] - (BYTE)src2f->imageData[(y+j)*w+(x+i)]);
                        count++;
                    }
                }
            }

            mad /= count;
            //mad /= (block*block);

            if(mad >= th) dest->imageData[(y/block)*(w/block)+(x/block)] = (char)255;
            else          dest->imageData[(y/block)*(w/block)+(x/block)] = (char)0;
        }
    }
    // GIST end 201001

    IplConvKernel* element = cvCreateStructuringElementEx(3, 3, 1, 1, CV_SHAPE_RECT, 0);
    cvDilate(dest, dest, element);
    //cvSaveImage("motionmap.bmp", dest);

    // Cleanup
    cvReleaseStructuringElement(&element);
    SAFE_RELEASE_IMAGE(src1f)
    SAFE_RELEASE_IMAGE(src2f)
} //getMotionMap_block


//-----------------------------------------------------------------------------
// Reference�� ���� Depth estimation error �� �����. 
// �������� ���·� Reference depth�� �־���
//
int CEstimation::update_error_cost(CParameterDepthEstimation *cParameter, CIYuv<DepthType> *yuvRefDepth, CIYuv<ImageType> *yuvCenter, bool isfirstframe, int GC_cycle, bool temporalonoff)
{
    int iCycle = GC_cycle;
    int i, j, pp;
    int d;
    BYTE t;
    double slope = 1.00;
    int BlkSize  = 16;
    IplImage* motionmap = NULL;
    IplImage* Ycurr     = NULL;
    double threshold = cParameter->getThreshold();

    memset(labels, 0, m_iPicsize*sizeof(DepthType));

    // Convert current frame to grayscale iplImage
    Ycurr = cvCreateImage(cvSize(m_iWidth, m_iHeight), IPL_DEPTH_8U, 1);
    for(int y=0; y<m_iHeight; y++) {
        for(int x=0; x<m_iWidth; x++){
            Ycurr->imageData[y*m_iWidth+x] = yuvCenter->Y[y][x];
        }
    }

    // Reference depth mode
    if(cParameter->getDEmode() == 3) {
        double Zn = m_dZnear;
        double Zf = m_dZfar;
        IplImage* refDepth = cvCreateImage(cvSize(m_iWidth, m_iHeight),8,1);
        cvSet(refDepth,cvScalarAll(0));
        int ws = refDepth->widthStep;
        int str = m_dRefBaseline<0 ? m_iWidth:0;
        int stp = m_dRefBaseline<0 ? 0:m_iWidth;
        int cnt = m_dRefBaseline<0 ? -1:+1;

        printf("Updating Error cost function (reference depth mode)...\n");
        for(int my = 0; my < m_iHeight; my++)
          for(int mx = str; mx != stp; mx += cnt) {
            BYTE   I  = (BYTE)yuvRefDepth->Y[my][mx];
            double Z  = 1/((I/(MAX_DEPTH-1)) * (1/Zn - 1/Zf) + (1/Zf));
            double dd = (m_dRefFocal_length*m_dRefBaseline/Z)-m_dOffset;
            int td = cvRound(dd);

           if ((mx-td) >= 0 && (mx-td < m_iWidth))
                refDepth->imageData[my*ws + mx-td] = (unsigned char)abs(cvRound(m_dRefRatio*m_iPrecision*dd));
          }
        // update error cost
        for(j = pp = 0; j < m_iHeight; j++) {
            for(i = 0; i < m_iWidth; i++, pp++) {
                d = (BYTE)refDepth->imageData[pp];
                if(d != 0) {
                    for (t=0; t < m_iNumOfLabels; t++) {
                        if(d ==  m_aiLabel2Disparity[0][t])
#ifdef POZNAN_TWO_COSTS
                        {
                            errors_left[t][pp] = 0;
                            errors_right[t][pp] = 0;
                        }
#else
                            errors[t][pp] = 0;
#endif
                        else
#ifdef POZNAN_TWO_COSTS
                        {
                            errors_left[t][pp] += abs(d - m_aiLabel2Disparity[0][t])+1;
                            errors_right[t][pp] += abs(d - m_aiLabel2Disparity[1][t])+1; //Not shure whether 1 or 0 it was 0 
                        }
#else
                            errors[t][pp] += abs(d - m_aiLabel2Disparity[0][t])+1;
#endif
                    }
                }
            }
        }
        cvReleaseImage(&refDepth);
        iCycle = 1;

    }else if( cParameter->getDEmode() == 2 && gc_mode == SEMI_TEMPORAL)
    {
        printf("Updating Error cost function (semi-automatic mode 2)...\n");
        update_error_cost_etri(cParameter, Ycurr);
        iCycle = 1;

    } else if(gc_mode == SEMI_TEMPORAL || (gc_mode == GC_AUTO && !isfirstframe && temporalonoff) )
    {
        printf("Motion search...\n");

        // GIST start 201001
        int motionmap_width, motionmap_height;
        motionmap_width  = m_iWidth/BlkSize;
        motionmap_height = m_iHeight/BlkSize;

        if(m_iWidth%BlkSize != 0)   motionmap_width += 1;
        if(m_iHeight%BlkSize != 0)  motionmap_height += 1;

        motionmap = cvCreateImage(cvSize(motionmap_width,motionmap_height), IPL_DEPTH_8U, 1);
        // GIST end 201001

        getMotionMap_block(Ycurr, ipl_Yprev, motionmap, BlkSize, cParameter->getThreshold());

        printf("Updating Error cost function...\n");
        for(j = pp = 0; j < m_iHeight; j++)
        {
          for(i = 0; i < m_iWidth; i++, pp++)
          {

              if( cParameter->getDEmode() == 1 && ipl_manstatic != NULL && ipl_manstatic->imageData[pp] ){
                  t = (BYTE)ipl_Drefresh->imageData[pp];
                  d = depth2label[t];
                  for(t=0; t<m_iNumOfLabels; t++)  {
#ifdef POZNAN_TWO_COSTS
                      errors_left[t][pp] = 2*errors_left[t][pp];
                      errors_right[t][pp] = 2*errors_right[t][pp];
#else
                      errors[t][pp] = 2*errors[t][pp];
#endif
                  }
#ifdef POZNAN_TWO_COSTS
                  errors_left[d][pp]  = 0;
                  errors_right[d][pp]  = 0;
#else
                  errors[d][pp]  = 0;
#endif
              } else if(motionmap->imageData[(j/BlkSize)*(m_iWidth/BlkSize) + (i/BlkSize)] == 0) {
                  for(d=0; d<m_iNumOfLabels; d++)  {
                      if (d == labels_prev[pp])
#ifdef POZNAN_TWO_COSTS
                      {
                          errors_left[d][pp] = 0;
                          errors_right[d][pp] = 0;
                      }
#else
                          errors[d][pp] = 0;
#endif
                      else
#ifdef POZNAN_TWO_COSTS
                      {
                          errors_left[d][pp] += (CostType)slope*abs(d - labels_prev[pp])+1;
                          errors_right[d][pp] += (CostType)slope*abs(d - labels_prev[pp])+1;
                       }   
#else
                          errors[d][pp] += (CostType)slope*abs(d - labels_prev[pp])+1;
#endif
                  }
              }
          }//for i
        } //for j
        iCycle = 1;

    } else if (gc_mode == SEMI_D_E_REFRESH)
    {
        // Use MDM as manual Disparity map
        for(pp = 0; pp< m_iPicsize; pp++) {
            d = m_iPrecision*(DepthType)ipl_mandisp->imageData[pp];
            if(d != 0) {           // If Manual disparity defined
                for (t=0; t < m_iNumOfLabels; t++) {
                    if (d == m_aiLabel2Disparity[0][t]) {
#ifdef POZNAN_TWO_COSTS
                      errors_left[t][pp] = 0;
                      errors_right[t][pp] = 0;
#else
                      errors[t][pp] = 0;
#endif
                    } else {
#ifdef POZNAN_TWO_COSTS
                      errors_left[t][pp] = 2*errors_left[t][pp];
                      errors_right[t][pp] = 2*errors_right[t][pp];
#else
                      errors[t][pp] = 2*errors[t][pp];
#endif
                    }
                }
            }
        }
        iCycle = 1;

    } else if (gc_mode == SEMI_D_REFRESH) {
        // Use MDM as manual DEPTH map
        for(pp = 0; pp< m_iPicsize; pp++) {
            t = (DepthType)ipl_mandisp->imageData[pp];
            d = depth2label[t];
            labels[pp]         = d;
        }
        cvCopy(ipl_mandisp, ipl_Drefresh);
        iCycle = 0;
    }

    //Swap curr/prev image pointers
    IplImage* temp = ipl_Yprev;
    ipl_Yprev = Ycurr;
    Ycurr = temp;
    SAFE_RELEASE_IMAGE(Ycurr)
    SAFE_RELEASE_IMAGE(motionmap)

    return iCycle;
} //CEstimation::update_error_cost

//-----------------------------------------------------------------------------
int CEstimation::update_error_cost_etri(CParameterDepthEstimation *cParameter, IplImage* ipl_Ycurr)
{
    int i, j, x, y, pp;
    int blocksize;

    int width = cParameter->getSourceWidth();
    int height = cParameter->getSourceHeight();

    CIYuv<DepthType> yuvMoveDepth;
    yuvMoveDepth.resize(height, width, POZNAN_DEPTHMAP_CHROMA_FORMAT); // alloc

    switch (cParameter->getMotionSearchBSize())
    {
    case 0:
        blocksize = 4;
        break;
    case 1:
        blocksize = 8;
        break;
    case 2:
        blocksize = 16;
    }

    BYTE **PreviousY, **CurrentY; //Owieczka
    PreviousY = new BYTE *[height];
    CurrentY  = new BYTE *[height];
    for (i=0; i<height; i++)
    {
        PreviousY[i] = new BYTE [width];
        CurrentY[i]  = new BYTE [width];
    }

    // -------------- store color value to find motion vector ----------------
    for (j=0; j<height; j++)
        for (i=0; i<width; i++) PreviousY[j][i] = (BYTE)ipl_Yprev->imageData[j*width+i];

    for (j=0; j<height; j++)
        for (i=0; i<width; i++) CurrentY[j][i] = (BYTE)ipl_Ycurr->imageData[j*width+i];


    /* ---------------------- Depth Frame Estimation by FullSearch ---------------------- */
    int     size = (int)(height*width/blocksize/blocksize);
    int     SearchX, SearchY;
    int     num;

    int     **motion_v1, **motion_v2, **counting;
    motion_v1 = new int *[2];
    motion_v2 = new int *[2];

    for (j=0; j<height; j++)
        for (i=0; i<width; i++)
        {
            yuvMoveDepth.Y[j][i] = 0;
        }

    counting  = new int *[height];
    for (i=0; i<height; i++) counting[i] = new int [width];

    for (i=0; i<2; i++)
    {
        motion_v1[i] = new int [size];
        motion_v2[i] = new int [size];
    }

    for (j=0; j<height; j++)
        for (i=0; i<width; i++)
        {
            counting[j][i] = 0;
        }

    FullSearch(PreviousY, CurrentY, motion_v1[0], motion_v2[0], height, width, blocksize);

    num = 0;
    int min_mv = 10000;
    for (j=0; j<height; j+=blocksize)
        for (i=0; i<width; i+=blocksize)
        {
            SearchX = motion_v1[0][num];
            SearchY = motion_v2[0][num];

            for (y=j; y<j+blocksize; y++)
                for (x=i; x<i+blocksize; x++)
                {
                    if((SearchX == 0) && (SearchY == 0)) counting[y+SearchY][x+SearchX] = 1;
                }
            num++;
        }

    for (j=pp=0; j<height; j++)
        for (i=0; i<width; i++, pp++)
        {
            if(counting[j][i] != 1) {
                yuvMoveDepth.Y[j][i] = MAX_DEPTH-1;                  //Calculate depth
            } else {
              //yuvMoveDepth.Y[j][i] = yuvInitDepth.Y[j][i]; //Previous frame depth
                yuvMoveDepth.Y[j][i] = m_acLabel2Depth[labels_prev[pp]];
            }
        }

    for (i=0; i<height; i++) delete[] counting[i];
    delete[] counting;

    for(j = pp = 0; j < m_iHeight; j++) {
        for(i = 0; i < m_iWidth; i++, pp++) {
            if(yuvMoveDepth.Y[j][i] != (MAX_DEPTH-1)) {
                for(int d=0; d<m_iNumOfLabels; d++)  {
                    if(yuvMoveDepth.Y[j][i] == m_acLabel2Depth[d]) {
#ifdef POZNAN_TWO_COSTS
                        errors_left[d][pp] = 0;
                        errors_right[d][pp] = 0;
#else
                        errors[d][pp] = 0;
#endif
                    }
                }
            }
        }
    }

    //Free memory
    for (i=0; i<2; i++)
    {
        delete[] motion_v1[i];
        delete[] motion_v2[i];
    }
    delete[] motion_v1;
    delete[] motion_v2;

    for (i=0; i<height; i++)
    {
        delete[] PreviousY[i];
        delete[] CurrentY[i];
    }
    delete[] PreviousY;
    delete[] CurrentY;
    yuvMoveDepth.~CIYuv(); // release

    return 0;
} //update_error_cost_etri

//-----------------------------------------------------------------------------
void CEstimation::depth_estimation_post_processing_etri(DepthType **pDepth, CParameterDepthEstimation *cParameter, int **FirstFrame, int **SecondFrame)
{
    int j,i,pp, x, y;


    if (gc_mode != SEMI_D_REFRESH && gc_mode != SEMI_D_E_REFRESH) {
        // ETRI start
        //---------------- post processing -------------------------------
        int new_height = m_iHeight;
#ifdef SUB_PEL_PRECISION
        int new_width  = m_iWidth * cParameter->getPrecision();
#else
        int new_width = m_iWidth;
#endif
        int **tempB, **tempP;
        tempB = new int *[new_height];
        tempP = new int *[new_height];

        for (i=0; i<new_height; i++)
        {
            tempB[i] = new int [new_width];
            tempP[i] = new int [new_width];
        }

        // -------------- Get frame difference --------------
        int Threshold   = cParameter->getThresOfDepthDiff();
        for (j=0; j<new_height; j++)
            for (i=0; i<new_width; i++)
            {
                tempP[j][i] = abs(FirstFrame[j][i] - SecondFrame[j][i]);
            }

        for (j=0; j<new_height; j++)
            for (i=0; i<new_width; i++)
            {
                if (tempP[j][i] > Threshold) tempP[j][i] = 255;
                else tempP[j][i] = 0;
            }

        // -------------- convert tempP(moving pixel mep) into tempB(moving block map) --------------
        CIYuv<ImageType> yuvTemp;
        yuvTemp.resize(m_iHeight, m_iWidth*cParameter->getPrecision(), POZNAN_DEPTHMAP_CHROMA_FORMAT);  // Alloc memory
        int BSize;
        switch (cParameter->getMovingObjectsBSize())
        {
        case 0:
            BSize = 32;
            break;
        case 1:
            BSize = 64;
            break;
        case 2:
            BSize = 128;
        }
        MakeBlockMap(tempP, tempB, PrevTempB, new_height, new_width, cParameter->getPrecision(), BSize);

        for (j=0; j<new_height; j++)
        {
            for (i=0; i<new_width; i++)
            {
                yuvTemp.Y[j][i] = 0;
                if (tempB[j][i] == (MAX_LUMA-1)) yuvTemp.Y[j][i] = (MAX_LUMA-1);
            }
        }

        // --------------------- post processing depth map --------------------
        int blocksize;
        switch (cParameter->getMotionSearchBSize())
        {
        case 0:
            blocksize = 4;
            break;
        case 1:
            blocksize = 8;
            break;
        case 2:
            blocksize = 16;
        }

        int nPrecision;
        for (j=0; j<m_iHeight; j+=blocksize) {
            for (i=0; i<m_iWidth; i+=blocksize) {

                for (y=j; y<j+blocksize; y++)
                    for (x=i; x<i+blocksize; x++)
                    {
                        for (nPrecision=0; nPrecision<(int)cParameter->getPrecision(); nPrecision++)
                        {
                            if (yuvTemp.Y[y][x*cParameter->getPrecision()+nPrecision] == (MAX_LUMA-1))
                            {
                                //yuvDepth.Y[y][x] = yuvPreDepth.Y[y][x]; // Calculated Depth
                            }
                            else
                            {
                                //yuvDepth.Y[y][x] = yuvInitDepth.Y[y][x];
                                pDepth[y][x] = m_acLabel2Depth[labels_prev[y*m_iWidth+x]]; // Prev.Frame Depth
                            }
                        }
                    }
            }//i
        }//j

        median(m_iHeight, m_iWidth, 5, pDepth);

        //Free memory
        for (i=0; i<new_height; i++)
        {
            delete[] tempB[i];
            delete[] tempP[i];
        }
        delete[] tempB;
        delete[] tempP;
        yuvTemp.~CIYuv(); // release
        // ETRI end
    }

    for(j=pp=0; j<m_iHeight; j++)   {
        for(i=0; i<m_iWidth; i++,pp++) {
            labels_prev[pp] = (BYTE)depth2label[pDepth[j][i]];
        }
    }
}

//-----------------------------------------------------------------------------
void CEstimation::depth_estimation_post_processing(DepthType **pDepth, CParameterDepthEstimation *cParameter)
{
    int j,i,pp;

        // Copy Results and depth map
        for(j=pp=0; j<m_iHeight; j++)
        {
            for(i=0; i<m_iWidth; i++,pp++)
            {

                if (gc_mode == SEMI_D_E_REFRESH)
                    ipl_Drefresh->imageData[pp] = m_acLabel2Depth[labels[pp]];

                if ( (gc_mode == SEMI_TEMPORAL) && ipl_manstatic != NULL && ipl_manstatic->imageData[pp]
                    )
                    pDepth[j][i] = (DepthType)ipl_Drefresh->imageData[pp];
                else
                    pDepth[j][i] = m_acLabel2Depth[labels[pp]];
            }
        }

        median(m_iHeight, m_iWidth, 3, pDepth);

    for(j=pp=0; j<m_iHeight; j++)   {
        for(i=0; i<m_iWidth; i++,pp++) {
            labels_prev[pp] = (DepthType)depth2label[pDepth[j][i]];
        }
    }
} //CEstimation::depth_estimation_post_processing

//-----------------------------------------------------------------------------
//Graph cut optimisation for semi-automatic modes
void CEstimation::depth_estimation_by_graph_cut_semi(DepthType **pDepth, int iCycle, CIYuv<ImageType> *yuvCenter)
{
    int i, j, c, pp, source;
    int right, down;
    int cost_cur_temp, cost_cur_temp2, cost_right_temp, cost_down_temp;

    if (iCycle > 0)
        printf("Graph Cuts\n");
    else
        printf("Skipping Graph Cuts....\n");

    for(c=0; c<iCycle; c++)
    {
        for(source=0; source<m_iNumOfLabels; source+=1)
        {
            Graph *g = new Graph();
            int counter = 0;
            for(pp = 0; pp< m_iPicsize; pp++)
            {
                nodes[pp] = g -> add_node();

#ifdef POZNAN_TWO_COSTS
                if(labels[pp] == source)
                    g -> set_tweights(nodes[pp], MINSEL(errors_left[source][pp],errors_right[source][pp]), COST_INF);
                else
                    g -> set_tweights(nodes[pp], MINSEL(errors_left[source][pp],errors_right[source][pp]), MINSEL(errors_left[labels[pp]][pp],errors_right[labels[pp]][pp]));
#else
                if(labels[pp] == source)
                    g -> set_tweights(nodes[pp], errors[source][pp], COST_INF);
                else
                    g -> set_tweights(nodes[pp], errors[source][pp], errors[labels[pp]][pp]);
#endif
            }

            double Ls = 0.1;

            for(j = pp = 0; j < m_iHeight; j++)
            {
                for(i = 0; i < m_iWidth; i++, pp++)
                {
                    // set condition
                    right = pp+1;
                    down = pp+m_iWidth;

                    cost_cur_temp = m_aiEdgeCost[byte_abs[labels[pp]-source]];
                    cost_cur_temp2 = m_aiEdgeCost2[byte_abs[labels[pp]-source]];

                    // add auxiliary node
                    if(i != m_iWidth_minus1)
                    {
                        if( gc_mode == SEMI_D_E_REFRESH &&
                               (ipl_manedge->imageData[j*ipl_manedge->widthStep + i*ipl_manedge->nChannels] ||
                                ipl_manedge->imageData[j*ipl_manedge->widthStep + (i+1)*ipl_manedge->nChannels])
                           ) Ls = 0.1;
                        else Ls = 1;

                        if(labels[pp] != labels[right])
                        {
                            auxiliary[counter] = g -> add_node();
                            cost_right_temp = m_aiEdgeCost[byte_abs[labels[right]-source]];
                            if(m_iImageSegmentation == 1 && !refresh_frame() && segmlabel[j * m_iWidth + i] != segmlabel[j * m_iWidth + i+1])
                            {
                                int cost_right_temp2 = m_aiEdgeCost2[byte_abs[labels[right]-source]];
                                if(m_aiEdgeCost2[byte_abs[labels[pp] - labels[right]]]!=0){
                                  g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost2[byte_abs[labels[pp] - labels[right]]]);
                                  g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp2, cost_cur_temp2);
                                  g -> add_edge(auxiliary[counter], nodes[right], cost_right_temp2, cost_right_temp2);
                                }
                                else
                                  g -> add_edge(nodes[pp], nodes[right], cost_cur_temp2, cost_cur_temp2);
                            }
                            else
                            {
                                g -> set_tweights(auxiliary[counter], 0, Ls*m_aiEdgeCost[byte_abs[labels[pp] - labels[right]]]);
                                g -> add_edge(nodes[pp], auxiliary[counter], Ls*cost_cur_temp, Ls*cost_cur_temp);
                                g -> add_edge(auxiliary[counter], nodes[right], Ls*cost_right_temp, Ls*cost_right_temp);
                            }
                          counter++;
                        }
                        else
                        {
                            g -> add_edge(nodes[pp], nodes[right], Ls*cost_cur_temp, Ls*cost_cur_temp);
                        }
                    }

                    // add auxiliary node - DOWN
                    if(j != m_iHeight_minus1)
                    {
                        if( gc_mode == SEMI_D_E_REFRESH &&
                               (ipl_manedge->imageData[j*ipl_manedge->widthStep + i*ipl_manedge->nChannels] ||
                                ipl_manedge->imageData[(j+1)*ipl_manedge->widthStep + i*ipl_manedge->nChannels])
                           ) Ls = 0.1;
                        else Ls = 1;

                        if(labels[pp] != labels[down])
                        {
                            auxiliary[counter] = g -> add_node();
                            cost_down_temp = m_aiEdgeCost[byte_abs[labels[down]-source]];

                            if(m_iImageSegmentation == 1 && !refresh_frame() && segmlabel[j * m_iWidth + i] != segmlabel[(j+1) * m_iWidth + i])
                            {
                                int cost_down_temp2 = m_aiEdgeCost2[byte_abs[labels[down]-source]];
                                if(m_aiEdgeCost2[byte_abs[labels[pp] - labels[down]]]!=0){
                                  g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost2[byte_abs[labels[pp] - labels[down]]]);
                                  g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp2, cost_cur_temp2);
                                  g -> add_edge(auxiliary[counter], nodes[down], cost_down_temp2, cost_down_temp2);
                                }
                                else
                                  g -> add_edge(nodes[pp], nodes[down], cost_cur_temp2, cost_cur_temp2);
                            }
                            else
                            {
                                g -> set_tweights(auxiliary[counter], 0, Ls*m_aiEdgeCost[byte_abs[labels[pp] - labels[down]]]);
                                g -> add_edge(nodes[pp], auxiliary[counter], Ls*cost_cur_temp, Ls*cost_cur_temp);
                                g -> add_edge(auxiliary[counter], nodes[down], Ls*cost_down_temp, Ls*cost_down_temp);
                            }
                            counter++;
                        }
                        else
                        {
                            g -> add_edge(nodes[pp], nodes[down], Ls*cost_cur_temp, Ls*cost_cur_temp);
                        }
                    }//if(j != m_iHeight_minus1)
                } //for i
            }
            printf(".");
            fflush(stdout);

#ifdef GRAPH_CUT_LOG
            printf("cycle:%d, label=(%d), auxiliary:%d\n", c+1, source, counter);
#endif
            Graph::flowtype flow = g -> maxflow();

            for(pp = 0; pp < m_iPicsize; pp++)
            {
                if (g->what_segment(nodes[pp]) != Graph::SOURCE)
                {
                    labels[pp] = (DepthType) source;
                }
            }

            delete g;
        } // source

        printf("Next Cycle\n");

    } // cycle

    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++,pp++)
        {
            pDepth[j][i] = m_acLabel2Depth[labels[pp]];
        }
    }
} // depth_estimation_by_graph_cut_semi
//-----------------------------------------------------------------------------

#ifdef POZNAN_OCC

inline double CEstimation::data_term(int label, int pp)
{
  return (this->*func_data_term)(label,pp);
}

double CEstimation::data_term_disparity(int label, int pp)
{
	/*
	return MINSEL(errors_left[label][pp],errors_right[label][pp]);
	//*/

	/*
	return ((occ_l<=label&&error_l[label][pp]<255)?error_l[label][pp]/2:1)+((occ_r<=label&&error_r[label][pp]<255)?error_r[label][pp]/2:1);
	//*/

	//Gladko��  p q r |p+r-2q|+|p-r|

	/*
	return ((occ_l<=label&&error_l[label][pp]<255)?error_l[label][pp]:10);
	//*/

	/*
	return ((ipl_manoccl->imageData[pp]<=m_aiLabel2Disparity[0][label]&&errors_left[label][pp]<255)?errors_left[label][pp]:10);
	//*/

	/* //Occ Gdy occ jest liczone jako cienie
	if(ipl_manoccl->imageData[pp]<=m_aiLabel2Disparity[0][label]&&errors_left[label][pp]<255)
	{
		//Lewy nie zas�oniety
		if(ipl_manoccr->imageData[pp]<=m_aiLabel2Disparity[1][label]&&errors_right[label][pp]<255)
		{
			//Lewy i przey nie zas�oni?y
			return (errors_left[label][pp]+errors_right[label][pp]+1)>>1;
		}
		else
		{
			//Tylko Lewy niezas�oniety
			return errors_left[label][pp];
		}
	}
	else
	{
		//Lewy zas�oni?y
		if(ipl_manoccr->imageData[pp]<=m_aiLabel2Disparity[1][label]&&errors_right[label][pp]<255)
		{
			//Prawy nie zas�oniety
			return errors_right[label][pp];
		}
		else
		{
			//Oba zas�oniete
			return 10;
		}

	}
	//*/

	//* //Occ Gdy Occ jest zsyntezowan?map?g�ebi
	int xl  = pp % m_iWidth + m_aiLabel2Disparity[0][label]/m_iPrecision;
	int xr  = pp % m_iWidth - m_aiLabel2Disparity[1][label]/m_iPrecision;
	int ppl = pp +  m_aiLabel2Disparity[0][label]/m_iPrecision;
	int ppr = pp -  m_aiLabel2Disparity[1][label]/m_iPrecision;
#ifdef POZNAN_OCC_ALWAYS_DEPTH
	if(xl>=0&&xl<m_iWidth&&disparity_map_left[ppl]<=m_acLabel2Depth[label]&&errors_left[label][pp]<COST_MAX-1)
#else
	if(xl>=0&&xl<m_iWidth&&disparity_map_left[ppl]<=m_aiLabel2Disparity[0][label]&&errors_left[label][pp]<COST_MAX)
#endif
	{
		//Lewy nie zas�oniety
#ifdef POZNAN_OCC_ALWAYS_DEPTH
		if(xr>=0&&xr<m_iWidth&&disparity_map_right[ppr]<=m_acLabel2Depth[label]&&errors_right[label][pp]<COST_MAX-1)
#else
		if(xr>=0&&xr<m_iWidth&&disparity_map_right[ppr]<=m_aiLabel2Disparity[1][label]&&errors_right[label][pp]<COST_MAX)
#endif
		{
			//Lewy i Prawy nie zas�oni?y
#ifdef POZNAN_DOUBLE_COST
			return (errors_left[label][pp]+errors_right[label][pp])/2;
#else
			return (errors_left[label][pp]+errors_right[label][pp]+1)/2;
#endif
		}
		else
		{
			//Tylko Lewy niezas�oniety
			return errors_left[label][pp];
		}
	}
	else
	{
		//Lewy zas�oni?y
#ifdef POZNAN_OCC_ALWAYS_DEPTH
		if(xr>=0&&xr<m_iWidth&&disparity_map_right[ppr]<=m_acLabel2Depth[label]&&errors_right[label][pp]<COST_MAX-1)
#else
		if(xr>=0&&xr<m_iWidth&&disparity_map_right[ppr]<=m_aiLabel2Disparity[1][label]&&errors_right[label][pp]<COST_MAX)
#endif
		{
			//Prawy nie zas�oniety
			return errors_right[label][pp];
		}
		else
		{
			//Oba zas�oniete
			return 10;
		}

	}
	//*/



	/* //Occ Gdy Occ jest zsyntezowan?map?g�ebi i odejmowanie kosztu
	int x = pp % m_iWidth;
	int xl  = pp % m_iWidth + m_aiLabel2Disparity[0][label]/m_iPrecision;
	int xr  = pp % m_iWidth - m_aiLabel2Disparity[1][label]/m_iPrecision;
	int ppl = pp +  m_aiLabel2Disparity[0][label]/m_iPrecision;
	int ppr = pp -  m_aiLabel2Disparity[1][label]/m_iPrecision;
	if(0&&xl>=0&&xl<m_iWidth&&ipl_manoccl->imageData[ppl]<=m_aiLabel2Disparity[0][label]&&errors_left[label][pp]<255)
	{
		//Lewy nie zas�oniety
		if(xr>=0&&xr<m_iWidth&&ipl_manoccr->imageData[ppr]<=m_aiLabel2Disparity[1][label]&&errors_right[label][pp]<255)
		{
			int e = (errors_left[label][pp]+errors_right[label][pp]+1)>>1;
			int mel = 0;
			int mer = 0;
			//Lewy i przey nie zas�oni?y
			for(int i=1;i<m_aiLabel2Disparity[0][label]/m_iPrecision;i++)
			{
				//Pobierz punkt w lewo
				int xt = x + i;
				int ppt = pp + i;
				int pptt;
				if(xt>=0&&xt<m_iWidth)
				{
					//Sprawdz jak sie resyntezuje
					xt = x + i + m_aiLabel2Disparity[0][labels[ppt]]/m_iPrecision;
					pptt = pp + i + m_aiLabel2Disparity[0][labels[ppt]]/m_iPrecision;
					if(ppl==pptt)
					{
						//je�li wskazuj?na ten sam punkt 
						//e += 10 - errors_left[labels[ppt]][ppt]/2;
						mel = max(mel,errors_left[labels[ppt]][ppt]-10);

					}
				}

				//Pobierz punkt w prawo
				xt = x - i;
				ppt = pp - i;
				pptt = 0;
				if(xt>=0&&xt<m_iWidth)
				{
					//Sprawdz jak sie resyntezuje
					xt = x - i - m_aiLabel2Disparity[1][labels[ppt]]/m_iPrecision;
					pptt = pp - i - m_aiLabel2Disparity[1][labels[ppt]]/m_iPrecision;
					if(ppr==pptt)
					{
						//je�li wskazuj?na ten sam punkt 
						//e += 10 - errors_right[labels[ppt]][ppt]/2;
						mer = max(mer,errors_right[labels[ppt]][ppt]-10);

					}
				}
			}

			return e-(mer+mel+1)>>1;
		}
		else
		{
			int e = errors_left[label][pp];
			int me = 0;
			//Tylko Lewy niezas�oniety
			for(int i=1;i<m_aiLabel2Disparity[0][label]/m_iPrecision;i++)
			{
				//Pobierz punkt w lewo
				int xt = x + i;
				int ppt = pp + i;
				int pptt;
				if(xt>=0&&xt<m_iWidth)
				{
					//Sprawdz jak sie resyntezuje
					xt = x + i + m_aiLabel2Disparity[0][labels[ppt]]/m_iPrecision;
					pptt = pp + i + m_aiLabel2Disparity[0][labels[ppt]]/m_iPrecision;
					if(ppl==pptt&&m_aiLabel2Disparity[0][label]>=m_aiLabel2Disparity[0][labels[ppt]])
					{
						//je�li wskazuj?na ten sam punkt 
						//e += 10 - errors_left[labels[ppt]][ppt];
						me = max(me,errors_left[labels[ppt]][ppt]-10);

					}
				}
			}
			return e-me;
		}
	}
	else
	{
		//Lewy zas�oni?y
		if(xr>=0&&xr<m_iWidth&&ipl_manoccr->imageData[ppr]<=m_aiLabel2Disparity[1][label]&&errors_right[label][pp]<255)
		{
			int e = errors_right[label][pp];
			int me = -1000;
			//Prawy nie zas�oniety
			for(int i=1;i<m_aiLabel2Disparity[1][label]/m_iPrecision;i++)
			{
				//Pobierz punkt w prawo
				int xt = x - i;
				int ppt = pp - i;
				int pptt;
				if(xt>=0&&xt<m_iWidth)
				{
					//Sprawdz jak sie resyntezuje
					xt = x - i - m_aiLabel2Disparity[1][labels[ppt]]/m_iPrecision;
					pptt = pp - i - m_aiLabel2Disparity[1][labels[ppt]]/m_iPrecision;
					if(ppr==pptt)
					{
						//je�li wskazuj?na ten sam punkt 
						//e += 10 - errors_right[labels[ppt]][ppt];
						me = max(me,errors_right[labels[ppt]][ppt]-10);

					}
				}
			}
			return max(e-me,5);
		}
		else
		{
			//Oba zas�oniete
			return 10;
		}

	}
	//*/
	
}

#define pow2(x) ((x)*(x))

double CEstimation::data_term_depth(int label, int pp)
{
	//* //Occ Gdy Occ jest zsyntezowan?map?g�ebi
  int target_pixel_ul,target_pixel_vl,target_pixel_ur,target_pixel_vr;
  CvMat* pix[2];
  pix[0] = cvCreateMat(3, 1, CV_64F);
  pix[1] = cvCreateMat(3, 1, CV_64F);
  cvmSet(pix[0], 2, 0, 1.0);
  cvmSet(pix[0], 1, 0, pp/m_iWidth);
  cvmSet(pix[0], 0, 0, pp%m_iWidth);
  //cvmMul(m_matH_V2L[labels[pp]], pix[0], pix[1]);
  cvmMul(m_matH_V2L[label], pix[0], pix[1]);
  //target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
  target_pixel_ul = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
  target_pixel_vl = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
  //cvmMul(m_matH_V2R[labels[pp]], pix[0], pix[1]);
  cvmMul(m_matH_V2R[label], pix[0], pix[1]);
  //target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
  target_pixel_ur = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
  target_pixel_vr = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);
  cvReleaseMat(&pix[0]);
  cvReleaseMat(&pix[1]);
	if( target_pixel_ul >= 0 && target_pixel_ul < m_iWidth && target_pixel_vl >= 0 && target_pixel_vl < m_iHeight && disparity_map_left[target_pixel_vl*m_iWidth+target_pixel_ul]<=m_acLabel2Depth[label]&&errors_left[label][pp]<COST_MAX)
	{
		//Lewy nie zas�oniety
		if(target_pixel_ur >= 0 && target_pixel_ur < m_iWidth && target_pixel_vr >= 0 && target_pixel_vr < m_iHeight && disparity_map_right[target_pixel_vr*m_iWidth+target_pixel_ur]<=m_acLabel2Depth[label]&&errors_right[label][pp]<COST_MAX)
		{
			//Lewy i Prawy nie zas�oni?y
#ifdef POZNAN_DOUBLE_COST
#ifdef POZNAN_OCC_POW_MEAN
      return sqrt((pow2(errors_left[label][pp])+pow2(errors_right[label][pp]))/2);
#else
			return (errors_left[label][pp]+errors_right[label][pp])/2;
#endif      
#else
			return (errors_left[label][pp]+errors_right[label][pp]+1)/2;
#endif
		}
		else
		{
			//Tylko Lewy niezas�oniety
			return errors_left[label][pp];
		}
	}
	else
	{
		//Lewy zas�oni?y
		if(target_pixel_ur >= 0 && target_pixel_ur < m_iWidth && target_pixel_vr >= 0 && target_pixel_vr < m_iHeight && disparity_map_right[target_pixel_vr*m_iWidth+target_pixel_ur]<=m_acLabel2Depth[label]&&errors_right[label][pp]<COST_MAX)
		{
			//Prawy nie zas�oniety
			return errors_right[label][pp];
		}
		else
		{
			//Oba zas�oniete
			return 10;
		}

	}	
}


//Owieczka
void CEstimation::depth_estimation_by_graph_cut_occ(DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter)
{
  int i, j, c, pp, source,x,y;
  int right, down;
  int cost_cur_temp, cost_right_temp, cost_down_temp;

  memset(labels, 0, sizeof(DepthType)*m_iPicsize);

  for(c=0; c<iCycle; c++)
  {
    for(source=0; source<m_iNumOfLabels; source+=1)
      //for(source=m_iNumOfLabels-1; source>=0; source-=1)
    {
      //Occ
      estimation_occ();		

      Graph *g = new Graph();
      int counter = 0;

      for(pp = 0; pp< m_iPicsize; pp++)
      {
        y = pp / m_iWidth;
        x = pp % m_iWidth;
        nodes[pp] = g -> add_node();
        
        if(labels[pp] == source)
          //g -> set_tweights(nodes[pp], MINSEL(errors_left[source][pp],errors_right[source][pp]), SHRT_MAX);
            g -> set_tweights(nodes[pp], data_term(source,pp), COST_INF);
        //g -> set_tweights(nodes[pp], ((ipl_manoccl->imageData[y*m_iWidth+x]<=source&&errors_left[source][pp]<255)?errors_left[source][pp]/2:1)+((ipl_manoccr->imageData[y*m_iWidth+x]<=source&&errors_right[source][pp]<255)?errors_right[source][pp]/2:1), SHRT_MAX);
        else
          //g -> set_tweights(nodes[pp], MINSEL(errors_left[source][pp],errors_right[source][pp]), MINSEL(errors_left[labels[pp]][pp],errors_right[labels[pp]][pp]));
          g -> set_tweights(nodes[pp], data_term(source,pp), data_term(labels[pp],pp));
        //g -> set_tweights(nodes[pp], ((ipl_manoccl->imageData[y*m_iWidth+x]<=source&&errors_left[source][pp]<255)?errors_left[source][pp]/2:1)+((ipl_manoccr->imageData[y*m_iWidth+x]<=source&&errors_right[source][pp]<255)?errors_right[source][pp]/2:1), ((ipl_manoccl->imageData[y*m_iWidth+x]<=labels[pp]&&errors_left[labels[pp]][pp]<255)?errors_left[labels[pp]][pp]/2:1)+((ipl_manoccr->imageData[y*m_iWidth+x]<=labels[pp]&&errors_right[labels[pp]][pp]<255)?errors_right[labels[pp]][pp]/2:1));

      }

      for(j = pp = 0; j < m_iHeight; j++)
      {
        for(i = 0; i < m_iWidth; i++, pp++)
        {
          // set condition
          right = pp+1;
          down = pp+m_iWidth;

          cost_cur_temp = m_aiEdgeCost[byte_abs[labels[pp]-source]];

          // add auxiliary node
          if(i != m_iWidth_minus1)
          {
            if(labels[pp] != labels[right])
            {
              auxiliary[counter] = g -> add_node();
              cost_right_temp = m_aiEdgeCost[byte_abs[labels[right]-source]];
              g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[right]]]);
              g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp, cost_cur_temp);
              g -> add_edge(auxiliary[counter], nodes[right], cost_right_temp, cost_right_temp);
              counter++;
            }
            else
            {
              g -> add_edge(nodes[pp], nodes[right], cost_cur_temp, cost_cur_temp);
            }
          }

          if(j != m_iHeight_minus1)
          {
            if(labels[pp] != labels[down])
            {
              auxiliary[counter] = g -> add_node();
              cost_down_temp = m_aiEdgeCost[byte_abs[labels[down]-source]];
              g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[down]]]);
              g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp, cost_cur_temp);
              g -> add_edge(auxiliary[counter], nodes[down], cost_down_temp, cost_down_temp);
              counter++;
            }
            else
            {
              g -> add_edge(nodes[pp], nodes[down], cost_cur_temp, cost_cur_temp);
            }
          }
        }
      }
      printf(".");
      fflush(stdout);

#ifdef GRAPH_CUT_LOG
      printf("cycle:%d, label=(%d), auxiliary:%d\n", c+1, source, counter);
#endif

      Graph::flowtype flow = g -> maxflow();

      for(pp = 0; pp < m_iPicsize; pp++)
      {
        if (g->what_segment(nodes[pp]) != Graph::SOURCE)
        {
          labels[pp] = (DepthType) source;
        }
      }

      delete g;
    } // source

    printf("Next Cycle\n");

  } // cycle


  for(j=pp=0; j<m_iHeight; j++)
  {
    for(i=0; i<m_iWidth; i++,pp++)
    {
      // GIST start
      //          labels_prev[pp] = labels[pp];  //now in post-processing function
      // GIST end
      pDepth[j][i] = m_acLabel2Depth[labels[pp]];
    }
  }
}
#endif

void CEstimation::depth_estimation_by_graph_cut(DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter) //Changed by SZK
{
    int i, j, c, pp, source;
    int right, down;
    int cost_cur_temp, cost_right_temp, cost_down_temp;

    memset(labels, 0, m_iPicsize*sizeof(DepthType));

    for(c=0; c<iCycle; c++)
    {
        for(source=0; source<m_iNumOfLabels; source+=1)
        {
            Graph *g = new Graph();
            int counter = 0;

            for(pp = 0; pp< m_iPicsize; pp++)
            {
                nodes[pp] = g -> add_node();

                if(labels[pp] == source)
#ifdef POZNAN_TWO_COSTS
                    g -> set_tweights(nodes[pp], MINSEL(errors_left[source][pp],errors_right[source][pp]), COST_INF);
#else
                    g -> set_tweights(nodes[pp], errors[source][pp], COST_INF);
#endif
                else
#ifdef POZNAN_TWO_COSTS
                    g -> set_tweights(nodes[pp], MINSEL(errors_left[source][pp],errors_right[source][pp]), MINSEL(errors_left[labels[pp]][pp],errors_right[labels[pp]][pp]));
#else
                    g -> set_tweights(nodes[pp], errors[source][pp], errors[labels[pp]][pp]);
#endif
            }

            for(j = pp = 0; j < m_iHeight; j++)
            {
                for(i = 0; i < m_iWidth; i++, pp++)
                {

                    // set condition
                    right = pp+1;
                    down = pp+m_iWidth;

                    cost_cur_temp = m_aiEdgeCost[byte_abs[labels[pp]-source]];

                    // add auxiliary node
                    if(i != m_iWidth_minus1)
                    {
                        if(labels[pp] != labels[right])
                        {
                            auxiliary[counter] = g -> add_node();
                            cost_right_temp = m_aiEdgeCost[byte_abs[labels[right]-source]];
                            g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[right]]]);
                            g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp, cost_cur_temp);
                            g -> add_edge(auxiliary[counter], nodes[right], cost_right_temp, cost_right_temp);
                            counter++;
                        }
                        else
                        {
                            g -> add_edge(nodes[pp], nodes[right], cost_cur_temp, cost_cur_temp);
                        }
                    }

                    if(j != m_iHeight_minus1)
                    {
                        if(labels[pp] != labels[down])
                        {
                            auxiliary[counter] = g -> add_node();
                            cost_down_temp = m_aiEdgeCost[byte_abs[labels[down]-source]];
                            g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[down]]]);
                            g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp, cost_cur_temp);
                            g -> add_edge(auxiliary[counter], nodes[down], cost_down_temp, cost_down_temp);
                            counter++;
                        }
                        else
                        {
                            g -> add_edge(nodes[pp], nodes[down], cost_cur_temp, cost_cur_temp);
                        }
                    }
                }
            }
            printf(".");
            fflush(stdout);

#ifdef GRAPH_CUT_LOG
            printf("cycle:%d, label=(%d), auxiliary:%d\n", c+1, source, counter);
#endif

            Graph::flowtype flow = g -> maxflow();

            for(pp = 0; pp < m_iPicsize; pp++)
            {
                if (g->what_segment(nodes[pp]) != Graph::SOURCE)
                {
                    labels[pp] = (DepthType) source;
                }
            }

            delete g;
        } // source

        printf("Next Cycle\n");

    } // cycle


    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++,pp++)
        {
            // GIST start
//          labels_prev[pp] = labels[pp];  //now in post-processing function
            // GIST end
            pDepth[j][i] = m_acLabel2Depth[labels[pp]];
        }
    }
}


//Nagoya start
//In segmented image, adjacent two pixels are in the same segment, smoothing coefficient is "SmoothingCoefficient".
//Otherwise, smoothing coefficient is the product of "SmoothingCoefficient" and "SmoothingCoefficient2".
//
//Smoothing coefficient is as follows:
//adjacent two pixels are in the same segment ==> (SmoothingCoefficient)*(SmoothingCoefficient2)
//Otherwise ==> SmoothingCoefficient
//
void CEstimation::depth_estimation_by_graph_cut_segmentation(DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter)
{
    int i, j, c, pp, source;
    int right, down;
    int cost_cur_temp, cost_right_temp, cost_down_temp;
    int cost_right_temp2, cost_down_temp2;
    int cost_cur_temp2;

    memset(labels, 0, m_iPicsize*sizeof(DepthType));

    for(c=0; c<iCycle; c++)
    {
        for(source=0; source<m_iNumOfLabels; source+=1)
        {
            Graph *g = new Graph();
            int counter = 0;

            for(pp = 0; pp< m_iPicsize; pp++)
            {
                nodes[pp] = g -> add_node();

                if(labels[pp] == source)
#ifdef POZNAN_TWO_COSTS
                    g -> set_tweights(nodes[pp], MINSEL(errors_left[source][pp],errors_right[source][pp]), COST_INF);
#else
                    g -> set_tweights(nodes[pp], errors[source][pp], COST_INF);
#endif
                else
#ifdef POZNAN_TWO_COSTS
                    g -> set_tweights(nodes[pp], MINSEL(errors_left[source][pp],errors_right[source][pp]),MINSEL(errors_left[labels[pp]][pp],errors_right[labels[pp]][pp]));
#else
                    g -> set_tweights(nodes[pp], errors[source][pp], errors[labels[pp]][pp]);
#endif
            }

            for(j = pp = 0; j < m_iHeight; j++)
            {
                for(i = 0; i < m_iWidth; i++, pp++)
                {
                    // set condition
                    right = pp+1;
                    down = pp+m_iWidth;

                    cost_cur_temp = m_aiEdgeCost[byte_abs[labels[pp]-source]];
                    cost_cur_temp2 = m_aiEdgeCost2[byte_abs[labels[pp]-source]];

                    // add auxiliary node
                    if(i != m_iWidth_minus1)
                    {
                        if(labels[pp] != labels[right])
                        {
                            auxiliary[counter] = g -> add_node();
                            cost_right_temp = m_aiEdgeCost[byte_abs[labels[right]-source]];
                            cost_right_temp2 = m_aiEdgeCost2[byte_abs[labels[right]-source]];
                            if(segmlabel[j * m_iWidth + i] != segmlabel[j * m_iWidth + i+1]){
                                if(m_aiEdgeCost2[byte_abs[labels[pp] - labels[right]]]!=0){
                                    g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost2[byte_abs[labels[pp] - labels[right]]]);
                                    g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp2, cost_cur_temp2);
                                    g -> add_edge(auxiliary[counter], nodes[right], cost_right_temp2, cost_right_temp2);
                                }
                                else
                                    g -> add_edge(nodes[pp], nodes[right], cost_cur_temp2, cost_cur_temp2);
                            }
                            else
                            {
                                g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[right]]]);
                                g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp, cost_cur_temp);
                                g -> add_edge(auxiliary[counter], nodes[right], cost_right_temp, cost_right_temp);
                            }
                            counter++;
                        }
                        else
                        {
//                          if(segmlabel[j * m_iWidth + i] != segmlabel[j * m_iWidth + i+1])
//                              g -> add_edge(nodes[pp], nodes[right], cost_cur_temp2, cost_cur_temp2);
//                          else
                                g -> add_edge(nodes[pp], nodes[right], cost_cur_temp, cost_cur_temp);
                        }
                    }

                    if(j != m_iHeight_minus1)
                    {
                        if(labels[pp] != labels[down])
                        {
                            auxiliary[counter] = g -> add_node();
                            cost_down_temp = m_aiEdgeCost[byte_abs[labels[down]-source]];
                            cost_down_temp2 = m_aiEdgeCost2[byte_abs[labels[down]-source]];
                            if(segmlabel[j * m_iWidth + i] != segmlabel[(j+1) * m_iWidth + i]){
                                if(m_aiEdgeCost2[byte_abs[labels[pp] - labels[down]]]!=0){
                                    g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost2[byte_abs[labels[pp] - labels[down]]]);
                                    g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp2, cost_cur_temp2);
                                    g -> add_edge(auxiliary[counter], nodes[down], cost_down_temp2, cost_down_temp2);
                                }
                                else
                                    g -> add_edge(nodes[pp], nodes[down], cost_cur_temp2, cost_cur_temp2);
                            }
                            else
                            {
                                g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[down]]]);
                                g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp, cost_cur_temp);
                                g -> add_edge(auxiliary[counter], nodes[down], cost_down_temp, cost_down_temp);
                            }
                            counter++;
                        }
                        else
                        {
//                          if(segmlabel[j * m_iWidth + i] != segmlabel[(j+1) * m_iWidth + i])
//                              g -> add_edge(nodes[pp], nodes[down], cost_cur_temp2, cost_cur_temp2);
//                          else
                                g -> add_edge(nodes[pp], nodes[down], cost_cur_temp, cost_cur_temp);
                        }
                    }
                }
            }
            printf(".");
            fflush(stdout);

#ifdef GRAPH_CUT_LOG
            printf("cycle:%d, label=(%d), auxiliary:%d\n", c+1, source, counter);
#endif

            Graph::flowtype flow = g -> maxflow();

            for(pp = 0; pp < m_iPicsize; pp++)
            {
                if (g->what_segment(nodes[pp]) != Graph::SOURCE)
                {
                    labels[pp] = (DepthType) source;
                }
            }

            delete g;
        } // source

        printf("Next Cycle\n");

    } // cycle

    //IplImage *dp1 = cvCreateImage(cvSize(m_iWidth, m_iHeight), IPL_DEPTH_8U, 1);
    //for(j=pp=0; j < m_iHeight; j++)
    //{
    //  for(i=0; i < m_iWidth; i++,pp++)
    //  {
    //      dp1->imageData[pp] = m_acLabel2Depth[labels[pp]];
    //  }
    //}
    //cvSaveImage("DepthGC.bmp", dp1);
    //cvReleaseImage (&dp1);

    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++,pp++)
        {
            // GIST start
//          labels_prev[pp] = labels[pp];
            // GIST end
            pDepth[j][i] = m_acLabel2Depth[labels[pp]];
        }
    }
}
//Nagoya end

#ifdef POZNAN_OCC
void CEstimation::depth_estimation_by_graph_cut_segmentation_occ(DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter)
{
    int i, j, c, pp, source;
    int right, down;
    int cost_cur_temp, cost_right_temp, cost_down_temp;
    int cost_right_temp2, cost_down_temp2;
    int cost_cur_temp2;

    memset(labels, 0, m_iPicsize*sizeof(DepthType));

    for(c=0; c<iCycle; c++)
    {
        for(source=0; source<m_iNumOfLabels; source+=1)
        {
            estimation_occ();		

            Graph *g = new Graph();
            int counter = 0;

            for(pp = 0; pp< m_iPicsize; pp++)
            {
                nodes[pp] = g -> add_node();

                if(labels[pp] == source)
                    g -> set_tweights(nodes[pp], data_term(source,pp), COST_INF);
                else
                    g -> set_tweights(nodes[pp], data_term(source,pp),data_term(labels[pp],pp));
            }

            for(j = pp = 0; j < m_iHeight; j++)
            {
                for(i = 0; i < m_iWidth; i++, pp++)
                {
                    // set condition
                    right = pp+1;
                    down = pp+m_iWidth;

                    cost_cur_temp = m_aiEdgeCost[byte_abs[labels[pp]-source]];
                    cost_cur_temp2 = m_aiEdgeCost2[byte_abs[labels[pp]-source]];

                    // add auxiliary node
                    if(i != m_iWidth_minus1)
                    {
                        if(labels[pp] != labels[right])
                        {
                            auxiliary[counter] = g -> add_node();
                            cost_right_temp = m_aiEdgeCost[byte_abs[labels[right]-source]];
                            cost_right_temp2 = m_aiEdgeCost2[byte_abs[labels[right]-source]];
                            if(segmlabel[j * m_iWidth + i] != segmlabel[j * m_iWidth + i+1]){
                                if(m_aiEdgeCost2[byte_abs[labels[pp] - labels[right]]]!=0){
                                    g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost2[byte_abs[labels[pp] - labels[right]]]);
                                    g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp2, cost_cur_temp2);
                                    g -> add_edge(auxiliary[counter], nodes[right], cost_right_temp2, cost_right_temp2);
                                }
                                else
                                    g -> add_edge(nodes[pp], nodes[right], cost_cur_temp2, cost_cur_temp2);
                            }
                            else
                            {
                                g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[right]]]);
                                g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp, cost_cur_temp);
                                g -> add_edge(auxiliary[counter], nodes[right], cost_right_temp, cost_right_temp);
                            }
                            counter++;
                        }
                        else
                        {
//                          if(segmlabel[j * m_iWidth + i] != segmlabel[j * m_iWidth + i+1])
//                              g -> add_edge(nodes[pp], nodes[right], cost_cur_temp2, cost_cur_temp2);
//                          else
                                g -> add_edge(nodes[pp], nodes[right], cost_cur_temp, cost_cur_temp);
                        }
                    }

                    if(j != m_iHeight_minus1)
                    {
                        if(labels[pp] != labels[down])
                        {
                            auxiliary[counter] = g -> add_node();
                            cost_down_temp = m_aiEdgeCost[byte_abs[labels[down]-source]];
                            cost_down_temp2 = m_aiEdgeCost2[byte_abs[labels[down]-source]];
                            if(segmlabel[j * m_iWidth + i] != segmlabel[(j+1) * m_iWidth + i]){
                                if(m_aiEdgeCost2[byte_abs[labels[pp] - labels[down]]]!=0){
                                    g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost2[byte_abs[labels[pp] - labels[down]]]);
                                    g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp2, cost_cur_temp2);
                                    g -> add_edge(auxiliary[counter], nodes[down], cost_down_temp2, cost_down_temp2);
                                }
                                else
                                    g -> add_edge(nodes[pp], nodes[down], cost_cur_temp2, cost_cur_temp2);
                            }
                            else
                            {
                                g -> set_tweights(auxiliary[counter], 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[down]]]);
                                g -> add_edge(nodes[pp], auxiliary[counter], cost_cur_temp, cost_cur_temp);
                                g -> add_edge(auxiliary[counter], nodes[down], cost_down_temp, cost_down_temp);
                            }
                            counter++;
                        }
                        else
                        {
//                          if(segmlabel[j * m_iWidth + i] != segmlabel[(j+1) * m_iWidth + i])
//                              g -> add_edge(nodes[pp], nodes[down], cost_cur_temp2, cost_cur_temp2);
//                          else
                                g -> add_edge(nodes[pp], nodes[down], cost_cur_temp, cost_cur_temp);
                        }
                    }
                }
            }
            printf(".");
            fflush(stdout);

#ifdef GRAPH_CUT_LOG
            printf("cycle:%d, label=(%d), auxiliary:%d\n", c+1, source, counter);
#endif

            Graph::flowtype flow = g -> maxflow();

            for(pp = 0; pp < m_iPicsize; pp++)
            {
                if (g->what_segment(nodes[pp]) != Graph::SOURCE)
                {
                    labels[pp] = (DepthType) source;
                }
            }

            delete g;
        } // source

        printf("Next Cycle\n");

    } // cycle

    //IplImage *dp1 = cvCreateImage(cvSize(m_iWidth, m_iHeight), IPL_DEPTH_8U, 1);
    //for(j=pp=0; j < m_iHeight; j++)
    //{
    //  for(i=0; i < m_iWidth; i++,pp++)
    //  {
    //      dp1->imageData[pp] = m_acLabel2Depth[labels[pp]];
    //  }
    //}
    //cvSaveImage("DepthGC.bmp", dp1);
    //cvReleaseImage (&dp1);

    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++,pp++)
        {
            // GIST start
//          labels_prev[pp] = labels[pp];
            // GIST end
            pDepth[j][i] = m_acLabel2Depth[labels[pp]];
        }
    }
}
//Nagoya end
#endif

#ifdef SEOULTECH_CUDA_SUPPORT
static void copy_error(CostType **errors, int* errors_cudacuts, int width, int height, int nLabels, double coeff)
{
    int n = 0, i = 0;
    for (int c = 0; c < nLabels; c++) {
		n = c;
		for (i = 0; i < width * height; i++) {
            errors_cudacuts[n] = coeff*errors[c][i];
			n += nLabels;
		}
	}
}

static void generate_smoothness_cost(int* smoothCostArray, int nLabels, double coeff)
{
	for(int i = 0; i < nLabels; i++)
	{
		for(int j = 0; j < nLabels; j++)
		{
			smoothCostArray[i*nLabels + j] = coeff*abs(i-j);
		}
	}
}

void CEstimation::depth_estimation_by_graph_cut_cuda(DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter, double datacoeff, double smoothcoeff, int is_stochatic)
{
    memset(labels, 0, m_iPicsize*sizeof(DepthType));

    // Setup Dataterm error and smoothness function for CudaCuts
    int* errors_cudacuts = (int*)malloc(sizeof(int)* m_iWidth * m_iHeight * m_iNumOfLabels);
    int* smooth_cost = (int*)malloc(sizeof(int)*m_iNumOfLabels * m_iNumOfLabels);
    copy_error(errors, errors_cudacuts, m_iWidth, m_iHeight, m_iNumOfLabels, datacoeff);
    generate_smoothness_cost(smooth_cost, m_iNumOfLabels, smoothcoeff);
    
    // Setup CudaCuts Memory
    CudaCuts cuts(m_iWidth, m_iHeight, m_iNumOfLabels, errors_cudacuts, smooth_cost);
    
    // Run CudaCuts
	std::vector<int> label_vec(m_iNumOfLabels);
	int n=0;
	generate(label_vec.begin(), label_vec.end(), [&n] { return n++;});
    cuts.run(label_vec, is_stochatic);

    int pp, row, col;
    for(row=pp=0; row<m_iHeight; row++)
    {
        for(col=0; col<m_iWidth; col++,pp++)
        {
            labels[pp] = (DepthType)cuts.pixelLabel[row*cuts.width + col];

            // GIST start
            // labels_prev[pp] = labels[pp];  //now in post-processing function
            // GIST end
            pDepth[row][col] = m_acLabel2Depth[labels[pp]];
        }
    }
}

void CEstimation::depth_estimation_by_graph_cut_no_auxnode(DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter) //Changed by SZK
{
    int i, j, c, pp, source;
    int right, down;
    int cost_cur_temp, cost_right_temp1, cost_right_temp2, cost_down_temp1, cost_down_temp2;

    memset(labels, 0, m_iPicsize*sizeof(DepthType));

    for(c=0; c<iCycle; c++)
    {
        for(source=0; source<m_iNumOfLabels; source+=1)
        {
            Graph *g = new Graph();

            for(pp = 0; pp< m_iPicsize; pp++)
            {
                nodes[pp] = g -> add_node();
                g -> set_tweights(nodes[pp], errors[source][pp], errors[labels[pp]][pp]);
            }

            for(j = pp = 0; j < m_iHeight; j++)
            {
                for(i = 0; i < m_iWidth; i++, pp++)
                {

                    // set condition
                    right = pp+1;
                    down = pp+m_iWidth;

                    cost_cur_temp = m_aiEdgeCost[byte_abs[labels[pp]-source]];

                    // no auxiliary node
                    if(i != m_iWidth_minus1)
                    {
                        if(labels[pp] != labels[right])
                        {
                            cost_right_temp1 = m_aiEdgeCost[byte_abs[labels[right]-source]];
                            cost_right_temp2 = m_aiEdgeCost[byte_abs[labels[pp]-source]];
                            g -> add_edge(nodes[pp], nodes[right], cost_right_temp1, cost_right_temp2);
                        }
                        else
                        {
                            g -> add_edge(nodes[pp], nodes[right], cost_cur_temp, cost_cur_temp);
                        }
                    }

                    if(j != m_iHeight_minus1)
                    {
                        if(labels[pp] != labels[down])
                        {
                            cost_down_temp1 = m_aiEdgeCost[byte_abs[labels[down]-source]];
                            cost_down_temp2 = m_aiEdgeCost[byte_abs[labels[pp]-source]];
                            g -> add_edge(nodes[pp], nodes[right], cost_down_temp1, cost_down_temp2);
                        }
                        else
                        {
                            g -> add_edge(nodes[pp], nodes[down], cost_cur_temp, cost_cur_temp);
                        }
                    }
                }
            }
            printf(".");
            fflush(stdout);

#ifdef GRAPH_CUT_LOG
            printf("cycle:%d, label=(%d), auxiliary:%d\n", c+1, source, counter);
#endif

            Graph::flowtype flow = g -> maxflow();

            for(pp = 0; pp < m_iPicsize; pp++)
            {
                if (g->what_segment(nodes[pp]) != Graph::SOURCE)
                {
                    labels[pp] = (DepthType) source;
                }
            }

            delete g;
        } // source

        printf("Next Cycle\n");

    } // cycle


    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++,pp++)
        {
            // GIST start
//          labels_prev[pp] = labels[pp];  //now in post-processing function
            // GIST end
            pDepth[j][i] = m_acLabel2Depth[labels[pp]];
        }
    }
}

static void push_relabel_tweight(PushRelabel::Graph &g, int node, int source, int sink, int source_cap, int sink_cap)
{
    int tr_cap = abs(source_cap - sink_cap);
    if(source_cap < sink)
    {
        g.addEdge(node, sink, tr_cap);
    }
    else    
    {
        g.addEdge(source, node, tr_cap);
    }
}

void CEstimation::depth_estimation_by_graph_cut_push_relabel(DepthType **pDepth, int iCycle, BYTE ***srcSEGM, CIYuv<ImageType> *yuvCenter)
{
    int i, j, c, pp, source;
    int right, down;
    int cost_cur_temp, cost_right_temp, cost_down_temp;

    int* nodes = (int*)malloc(m_iPicsize*sizeof(int));
    int* auxiliary = (int*)malloc((m_iPicsize<<1)*sizeof(int));

    memset(labels, 0, m_iPicsize*sizeof(DepthType));

    for(c=0; c<iCycle; c++)
    {
        for(source=0; source<m_iNumOfLabels; source+=1)
        {
            PushRelabel::Graph g;
            int source_node = g.addNode();
            int sink_node = g.addNode();

            int counter = 0;

            for(pp = 0; pp< m_iPicsize; pp++)
            {
                nodes[pp] = g.addNode();

                if(labels[pp] == source)
                    push_relabel_tweight(g, nodes[pp], source_node, sink_node, errors[source][pp], COST_INF);
                else
                    push_relabel_tweight(g, nodes[pp], source_node, sink_node, errors[source][pp], errors[labels[pp]][pp]);
            }

            for(j = pp = 0; j < m_iHeight; j++)
            {
                for(i = 0; i < m_iWidth; i++, pp++)
                {

                    // set condition
                    right = pp+1;
                    down = pp+m_iWidth;

                    cost_cur_temp = m_aiEdgeCost[byte_abs[labels[pp]-source]];

                    // add auxiliary node
                    if(i != m_iWidth_minus1)
                    {
                        if(labels[pp] != labels[right])
                        {
                            auxiliary[counter] = g.addNode();
                            cost_right_temp = m_aiEdgeCost[byte_abs[labels[right]-source]];
                            push_relabel_tweight(g, auxiliary[counter], source_node, sink_node, 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[right]]]);

                            g.addEdge(nodes[pp], auxiliary[counter], cost_cur_temp);
                            g.addEdge(auxiliary[counter], nodes[pp], cost_cur_temp);
                            g.addEdge(nodes[right], auxiliary[counter], cost_right_temp);
                            g.addEdge(auxiliary[counter], nodes[right], cost_right_temp);
                            counter++;
                        }
                        else
                        {
                            g.addEdge(nodes[pp], nodes[right], cost_cur_temp);
                            g.addEdge(nodes[right], nodes[pp], cost_cur_temp);
                        }
                    }

                    if(j != m_iHeight_minus1)
                    {
                        if(labels[pp] != labels[down])
                        {
                            auxiliary[counter] = g.addNode();
                            cost_down_temp = m_aiEdgeCost[byte_abs[labels[down]-source]];
                            push_relabel_tweight(g, auxiliary[counter], source_node, sink_node, 0, m_aiEdgeCost[byte_abs[labels[pp] - labels[down]]]);

                            g.addEdge(nodes[pp], auxiliary[counter], cost_cur_temp);
                            g.addEdge(auxiliary[counter], nodes[pp], cost_cur_temp);
                            g.addEdge(nodes[down], auxiliary[counter], cost_down_temp);
                            g.addEdge(auxiliary[counter], nodes[down], cost_down_temp);
                            counter++;
                        }
                        else
                        {
                            g.addEdge(nodes[pp], nodes[down], cost_cur_temp);
                            g.addEdge(nodes[down], nodes[pp], cost_cur_temp);
                        }
                    }
                }
            }

            cout << "graph construction done" << endl;
            cout << "Maximum flow is " << g.getMaxFlow(source_node, sink_node) << std::endl;
            // mincut 
            g.getMinCut(source_node);

            printf(".");
            fflush(stdout);

#ifdef GRAPH_CUT_LOG
            printf("cycle:%d, label=(%d), auxiliary:%d\n", c+1, source, counter);
#endif

            printf("max flow done");

            for(pp = 0; pp < m_iPicsize; pp++)
            {
                if(g.isVisited(nodes[pp]) != 1)
                {
                    labels[pp] = (DepthType) source;
                }
            }
        } // source
        printf("Next Cycle\n");

    } // cycle


    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++,pp++)
        {
            // GIST start
//          labels_prev[pp] = labels[pp];  //now in post-processing function
            // GIST end
            pDepth[j][i] = m_acLabel2Depth[labels[pp]];
        }
    }
}

#endif

//Nagoya start
void CEstimation::center_image_segmentation(BYTE ***srcSEGM)
{
    int i, j;
    int h, w;

    printf("Center image segmentation...\n");

    IplImage *src_img_segm = cvCreateImage(cvSize(m_iWidth, m_iHeight), IPL_DEPTH_8U, 3);
    IplImage *dst_img_segm = cvCreateImage(cvSize(m_iWidth, m_iHeight), IPL_DEPTH_8U, 3);

    cvZero(src_img_segm);
    cvZero(dst_img_segm);
    CvMat *clusters;
    CvMat *points;
    CvMat *color_mat = cvCreateMat (m_iMaxCluster, 1, CV_32FC3);
    CvMat *count_mat = cvCreateMat (m_iMaxCluster, 1, CV_32SC1);

    // Read Center Image  (OpenCV �޸𸮷� ����)
    for (j=0; j<m_iHeight; j++)
        for (i=0; i<m_iWidth; i++)
        {
            src_img_segm->imageData[(j * m_iWidth + i)*3]     = srcSEGM[0][j][i]; // e.g. B
            src_img_segm->imageData[(j * m_iWidth + i)*3 + 1] = srcSEGM[1][j][i]; //      G
            src_img_segm->imageData[(j * m_iWidth + i)*3 + 2] = srcSEGM[2][j][i]; //      R
        }

    dst_img_segm = cvCloneImage (src_img_segm);  // �� �̸� ����?
    clusters     = cvCreateMat  (m_iWidth * m_iHeight, 1, CV_32SC1);  // �ش� ũ����Ʈ id?
    points       = cvCreateMat  (m_iWidth * m_iHeight, 1, CV_32FC3);  // (R/G/B)

    int segm_level;
        CvRect roi;
        CvMemStorage *storage = 0;
        CvSeq *comp = 0;


    switch (m_iSegmentationMethod){
        case 1:
        segm_level=3; //Phyramid level
        roi.x = roi.y = 0;
        roi.width = src_img_segm->width & -(1 << segm_level);
        roi.height = src_img_segm->height & -(1 << segm_level);
        cvSetImageROI (src_img_segm, roi);
        dst_img_segm = cvCloneImage (src_img_segm);
        cvPyrMeanShiftFiltering (src_img_segm, dst_img_segm, 5.0, 5.0, segm_level, cvTermCriteria (CV_TERMCRIT_ITER + CV_TERMCRIT_EPS, 5, 1.0));
        break;

        case 2:
        segm_level=2; //Phyramid level
        double threshold1, threshold2;
        roi.x = roi.y = 0;
        roi.width = src_img_segm->width & -(1 << segm_level);
        roi.height = src_img_segm->height & -(1 << segm_level);
        cvSetImageROI (src_img_segm, roi);
        dst_img_segm = cvCloneImage (src_img_segm);
        storage = cvCreateMemStorage (0);
        threshold1 = 10.0;
        threshold2 = 10.0;
        cvPyrSegmentation(src_img_segm, dst_img_segm, storage, &comp, segm_level, threshold1, threshold2);
        break;

        case 3:
        // Input pixel value into a matrix 
		// @TODO OpenCV �� 1�������� ������ �Է��� �޴� ���� ����� �һ��
        cvSmooth(src_img_segm, src_img_segm, CV_MEDIAN, 3);
        for (j = 0; j < m_iHeight; j++)
        {
            for (i = 0; i < m_iWidth; i++)
            {
                points->data.fl[(j * m_iWidth + i) * 3]     = (uchar) src_img_segm->imageData[(j * m_iWidth + i) * 3];
                points->data.fl[(j * m_iWidth + i) * 3 + 1] = (uchar) src_img_segm->imageData[(j * m_iWidth + i) * 3 + 1];
                points->data.fl[(j * m_iWidth + i) * 3 + 2] = (uchar) src_img_segm->imageData[(j * m_iWidth + i) * 3 + 2];
            }
        }

        // Clustering (OpenCV ���)
        cvKMeans2 (points, m_iMaxCluster, clusters, cvTermCriteria (CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 100, 1.0));

        // Compute average for each cluster
        cvSetZero (color_mat);
        cvSetZero (count_mat);

        for (i = 0; i < m_iPicsize; i++)
        {
            int idx = clusters->data.i[i];
            int v = ++count_mat->data.i[idx];

            color_mat->data.fl[idx * 3]     = color_mat->data.fl[idx * 3]   * (v - 1) / v + points->data.fl[i * 3 + 0] / v;
            color_mat->data.fl[idx * 3 + 1] = color_mat->data.fl[idx * 3 + 1] * (v - 1) / v + points->data.fl[i * 3 + 1] / v;
            color_mat->data.fl[idx * 3 + 2] = color_mat->data.fl[idx * 3 + 2] * (v - 1) / v + points->data.fl[i * 3 + 2] / v;
        }

        // Draw color for each cluster
        for (i = 0; i < m_iPicsize; i++)
        {
            int idx = clusters->data.i[i];
            dst_img_segm->imageData[i * 3]      = (char) color_mat->data.fl[idx * 3];
            dst_img_segm->imageData[i * 3 + 1]  = (char) color_mat->data.fl[idx * 3 + 1];
            dst_img_segm->imageData[i * 3 + 2]  = (char) color_mat->data.fl[idx * 3 + 2];
        }

        break;

        default:
            fprintf(stderr, "Unknown number on SegmentationMethod [%d]\n", m_iSegmentationMethod);
        break;
    }

    for (j=0; j < m_iHeight; j++)
        for (i=0; i < m_iWidth; i++)
        {
            srcSEGM[0][j][i] = dst_img_segm->imageData[(j * m_iWidth + i)*3];
            srcSEGM[1][j][i] = dst_img_segm->imageData[(j * m_iWidth + i)*3 + 1];
            srcSEGM[2][j][i] = dst_img_segm->imageData[(j * m_iWidth + i)*3 + 2];
        }

    cvSaveImage("LowColorImage.bmp", dst_img_segm);  // segmentation ����� �̹����� ����?
	printf("The segmentation is done.\n");

	// ���⼭ ���ʹ� Labeling �ϴ� ������ ���ϴµ�. �� �̰� ���� �̷��� ���� ������?
	// OpenCV ���� ���������� �̹� ���� �� ������.
	printf("Labeling has started\n");

    int num_segm;

    unsigned int *headofsegm_x = (unsigned int *)malloc(m_iWidth * m_iHeight * sizeof(unsigned int));
    unsigned int *headofsegm_y = (unsigned int *)malloc(m_iWidth * m_iHeight * sizeof(unsigned int));
    unsigned int *segm_head_pixel = (unsigned int *)malloc(m_iWidth * m_iHeight * sizeof(unsigned int)) ;


    for(h=0; h < m_iHeight; h++)
        for (w=0; w < m_iWidth; w++)
        {
            segmlabel[h * m_iWidth + w]=0; //Label number of segment
            segm_head_pixel[h * m_iWidth + w] = 0;
            headofsegm_x[h * m_iWidth + w] = 0;
            headofsegm_y[h * m_iWidth + w] = 0;
        }

    num_segm=0; //The number of segment

    int color_thresh = 1;
    int block_h = m_iHeight;
    int block_w = m_iWidth;
    int tmp;

    //Count the number of pixels in each segment

    num_segm++; //num_segm = 1
    segmlabel[0]            =   num_segm; //Label 1 is assidned to pixel (0, 0).
    segm_head_pixel[0]      =   1; //Pixel (0, 0) is considered as head of segment.
    headofsegm_x[num_segm]  =   0;
    headofsegm_y[num_segm]  =   0;

    for (w = 1; w < m_iWidth; w++){ //0th row of image (h=0)
        //Compara the colors of adjacent pixels
        //srcSEGM[0][h][w]: Red of pixel (w, h)
        //srcSEGM[1][h][w]: Green of pixel (w, h)
        //srcSEGM[2][h][w]: Blue of pixel (w, h)
        //Is the difference between colors of pixels (w, 0) and (w-1, 0) smaller than threshold?
        if(abs(srcSEGM[0][0][w] - srcSEGM[0][headofsegm_y[segmlabel[w-1]]][headofsegm_x[segmlabel[w-1]]])<=color_thresh
            && abs(srcSEGM[1][0][w] - srcSEGM[1][headofsegm_y[segmlabel[w-1]]][headofsegm_x[segmlabel[w-1]]])<=color_thresh
            && abs(srcSEGM[2][0][w] - srcSEGM[2][headofsegm_y[segmlabel[w-1]]][headofsegm_x[segmlabel[w-1]]])<=color_thresh
        ){
            segmlabel[w]    =   segmlabel[w-1]; //If the difference of two colors is smaller than threshold, label of pixel (w, 0) is same as that of pixel(w-1, 0).
            segm_head_pixel[w]  =   0;
        } else { //If the difference of two colors is larger than threshold, new label is assigned to pixel (w, 0).
            num_segm++; //increment segment number
            segmlabel[w] =  num_segm; //assign new label number to new segment
            segm_head_pixel[w]  =   1; //Pixel (0, w) is considered as head of segment.
            headofsegm_x[num_segm]  =   w;
            headofsegm_y[num_segm]  =   0;
        }
    }

    //Not 0th row of image (h>=1)
    for(h= 1; h < m_iHeight; h++){
        //Compare the colors of pixels (0, h) and (0, h-1).
        //Is the difference of colors smaller than threshold?
        if( abs(srcSEGM[0][h][0] - srcSEGM[0][headofsegm_y[segmlabel[(h-1)*m_iWidth]]][headofsegm_x[segmlabel[(h-1)*m_iWidth]]])<=color_thresh
            && abs(srcSEGM[1][h][0] - srcSEGM[1][headofsegm_y[segmlabel[(h-1)*m_iWidth]]][headofsegm_x[segmlabel[(h-1)*m_iWidth]]])<=color_thresh
            && abs(srcSEGM[2][h][0] - srcSEGM[2][headofsegm_y[segmlabel[(h-1)*m_iWidth]]][headofsegm_x[segmlabel[(h-1)*m_iWidth]]])<=color_thresh
        ){
            segmlabel[h * m_iWidth] = segmlabel[(h - 1) * m_iWidth]; //If the difference of two colors is smaller than threshold, label of pixel (0, h) is same as that of pixel(0, h-1).
            segm_head_pixel[h * m_iWidth] = 0;
        } else {
            num_segm++;
            segmlabel[h * m_iWidth]=num_segm; //If the difference of two colors is larger than threshold, new label number is assigned to pixel (0, h).
            segm_head_pixel[h * m_iWidth] = 1;
            headofsegm_x[num_segm]=0;
            headofsegm_y[num_segm]=h;
        }

        //1 <= w < width
        for (w=1; w < m_iWidth; w++){
            //Compare the colors of pixels (w, h) and (w-1, h).
            if(abs(srcSEGM[0][h][w] - srcSEGM[0][headofsegm_y[segmlabel[h*m_iWidth+w-1]]][headofsegm_x[segmlabel[h*m_iWidth+w-1]]])<=color_thresh
                && abs(srcSEGM[1][h][w] - srcSEGM[1][headofsegm_y[segmlabel[h*m_iWidth+w-1]]][headofsegm_x[segmlabel[h*m_iWidth+w-1]]])<=color_thresh
                && abs(srcSEGM[2][h][w] - srcSEGM[2][headofsegm_y[segmlabel[h*m_iWidth+w-1]]][headofsegm_x[segmlabel[h*m_iWidth+w-1]]])<=color_thresh
            ){
                segmlabel[h * m_iWidth + w] = segmlabel[h * m_iWidth + w-1]; //If the difference of two colors is smaller than threshold, label of pixel (w, h) is same as that of pixel(w-1, h).
                segm_head_pixel[h * m_iWidth + w] = 0;

                //Compare the colors of pixels (w, h) and (w, h-1).
                if(abs(srcSEGM[0][h][w] - srcSEGM[0][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]])<=color_thresh
                    && abs(srcSEGM[1][h][w] - srcSEGM[1][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]])<=color_thresh
                    && abs(srcSEGM[2][h][w] - srcSEGM[2][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]])<=color_thresh
//                  && abs(srcSEGM[0][headofsegm_y[segmlabel[h*m_iWidth+w]]][headofsegm_x[segmlabel[h*m_iWidth+w]]] - srcSEGM[0][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]])<=color_thresh
//                  && abs(srcSEGM[1][headofsegm_y[segmlabel[h*m_iWidth+w]]][headofsegm_x[segmlabel[h*m_iWidth+w]]] - srcSEGM[1][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]])<=color_thresh
//                  && abs(srcSEGM[2][headofsegm_y[segmlabel[h*m_iWidth+w]]][headofsegm_x[segmlabel[h*m_iWidth+w]]] - srcSEGM[2][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]])<=color_thresh
                ){
                    if(segmlabel[h * m_iWidth + w] < segmlabel[(h-1) * m_iWidth + w])
                    {
                        segm_head_pixel[headofsegm_y[segmlabel[(h-1) * m_iWidth + w]] * m_iWidth + headofsegm_x[segmlabel[(h-1) * m_iWidth + w]]] = 0;
                        tmp = segmlabel[(h-1) * m_iWidth + w];
                        for(int kh=0; kh<=h*m_iWidth+w-1; kh++)
                        {
                            if(segmlabel[kh] == tmp)
                                segmlabel[kh] = segmlabel[h * m_iWidth + w];
                        }
                    }
                    else if (segmlabel[h * m_iWidth + w] > segmlabel[(h-1) * m_iWidth + w])
                    {
                        segm_head_pixel[headofsegm_y[segmlabel[h * m_iWidth + w]] * m_iWidth + headofsegm_x[segmlabel[h * m_iWidth + w]]] = 0;
                        tmp = segmlabel[h * m_iWidth + w];
                        segmlabel[h * m_iWidth + w] = segmlabel[(h-1) * m_iWidth + w];
                        for(int kh = 0; kh <= h*m_iWidth+w-1; kh++)
                        {
                            if(segmlabel[kh] == tmp)
                                segmlabel[kh] = segmlabel[(h-1) * m_iWidth + w];
                        }
                    }
                }
            } else if ( abs(srcSEGM[0][h][w] - srcSEGM[0][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]]) <= color_thresh
                        && abs(srcSEGM[1][h][w] - srcSEGM[1][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]]) <=  color_thresh
                        && abs(srcSEGM[2][h][w] - srcSEGM[2][headofsegm_y[segmlabel[(h-1)*m_iWidth+w]]][headofsegm_x[segmlabel[(h-1)*m_iWidth+w]]]) <=  color_thresh
            ){
                    segmlabel[h * m_iWidth + w] = segmlabel[(h-1) * m_iWidth + w];
                    segm_head_pixel[h * m_iWidth + w]=0;
            } else {
                    num_segm++;
                    segmlabel[h * m_iWidth + w]=num_segm;
                    headofsegm_x[num_segm]=w;
                    headofsegm_y[num_segm]=h;
                    segm_head_pixel[h * m_iWidth + w] = 1;
            }
        } // for w
        
		//printf(".");
		printf("\r%4d/%4d", h, m_iHeight);
        fflush(stdout);
    } //for h

    m_iNumOfSegms=0;

    for (i=0; i<m_iWidth * m_iHeight; i++)
        m_iNumOfSegms += segm_head_pixel[i];
    printf("\nNumber of segments = %d \n", m_iNumOfSegms);

    int head_segm=0;

    printf("Renumbering start\n");

    //Reassignment of segment number
    for (j=0; j<m_iHeight; j++)
    {
        for (i=0; i<m_iWidth; i++)
        {
            if(segm_head_pixel[j * m_iWidth + i]==1)
            {
                head_segm++;
                tmp = segmlabel[j * m_iWidth + i];
                for(int k= j * m_iWidth + i; k < m_iWidth * m_iHeight; k++)
//              for(int k= 0; k < m_iWidth * m_iHeight; k++)
                {
                    if(segmlabel[k] == tmp)
                        segmlabel[k] = head_segm;
                }
            }
        }
    }

    printf("Number of segments = %d, renumbering OK\n", head_segm);

    IplImage *color_assignment = cvCreateImage(cvSize(m_iWidth, m_iHeight), IPL_DEPTH_8U, 3);
    for(h=0; h < m_iHeight; h++){
        for (w=0; w < m_iWidth; w++){
            color_assignment->imageData[(h * m_iWidth + w)*3]   = segmlabel[h * m_iWidth + w] >> 16;
            color_assignment->imageData[(h * m_iWidth + w)*3 + 1] = (segmlabel[h * m_iWidth + w] - (segmlabel[h * m_iWidth + w] >> 16)) >> 8;
            color_assignment->imageData[(h * m_iWidth + w)*3 + 2] = segmlabel[h * m_iWidth + w] % (2<<8);
        }
    }
    cvSaveImage("ColorAssignedImage.bmp", color_assignment);  // 
    cvReleaseImage(&color_assignment);

    cvReleaseImage(&src_img_segm);
    cvReleaseImage(&dst_img_segm);
    cvReleaseMat(&clusters);
    cvReleaseMat(&points);
    cvReleaseMat(&color_mat);
    cvReleaseMat(&count_mat);

    free(segm_head_pixel);
    free(headofsegm_x);
    free(headofsegm_y);

    cvReleaseMemStorage(&storage);
    return;
}
//Nagoya end

//Nagoya start
void CEstimation::plane_fitting(DepthType **pDepth, BYTE ***srcSEGM, int num_segm)
{
    int i, j;

    printf("num of segms = %d\n", num_segm);
    num_segm=(int)(2 * num_segm);

    if( (reliability_matrix = (int *)malloc(num_segm * sizeof(int)))==NULL)
    {
        fprintf(stderr, "Can't allocate enough memory\n");
    }

    memset(reliability_matrix, 0, num_segm);

// memory allocation

    if( (matrix       = (long double **)malloc(num_segm*sizeof(long double *))) == NULL ||
        (Dmatrix      = (long double **)malloc(num_segm*sizeof(long double *))) == NULL ||
        (inv_matrix   = (long double **)malloc(num_segm*sizeof(long double *))) == NULL ||
        (prod_matrix  = (long double **)malloc(num_segm*sizeof(long double *))) == NULL ||
        (check_matrix = (long double **)malloc(num_segm*sizeof(long double *))) == NULL )
    {
        fprintf(stderr, "Can't allocate enough memory\n");
    }

    if( (matrix[0]       = (long double *)malloc(num_segm * 9 * sizeof(long double))) == NULL ||
        (Dmatrix[0]      = (long double *)malloc(num_segm * 3 * sizeof(long double))) == NULL ||
        (inv_matrix[0]   = (long double *)malloc(num_segm * 9 * sizeof(long double))) == NULL ||
        (prod_matrix[0]  = (long double *)malloc(num_segm * 3 * sizeof(long double))) == NULL ||
        (check_matrix[0] = (long double *)malloc(num_segm * 9 * sizeof(long double))) == NULL
    )
    {
        fprintf(stderr, "Can't allocate enough memory\n");
    }

    int pos;

    for(i=1, pos = 9; i < num_segm; i++, pos += 9)
    {
        matrix[i] = &matrix[0][pos];
        inv_matrix[i] = &inv_matrix[0][pos];
        check_matrix[i] = &check_matrix[0][pos];
    }
    for(i=1, pos = 3; i < num_segm; i++, pos += 3)
    {
        Dmatrix[i] = &Dmatrix[0][pos];
        prod_matrix[i] = &prod_matrix[0][pos];
    }

    memset(matrix[0],       0.0, num_segm * 9);
    memset(Dmatrix[0],      0.0, num_segm * 3);
    memset(inv_matrix[0],   0.0, num_segm * 9);
    memset(prod_matrix[0],  0.0, num_segm * 3);
    memset(check_matrix[0], 0.0, num_segm * 9);

    int nrow = 3; //The number of rows of matrix
    int ncol = 3; //The number of columns of matrix

    CvMat *src_matrix, *dst_matrix, *mul_matrix;
    CvMat *src_vec;
    double det;

    for(j = 0; j < m_iHeight; j++)
    {
        for (i = 0; i < m_iWidth; i++)
        {
            matrix[segmlabel[j * m_iWidth + i] - 1][0] += i * i;
            matrix[segmlabel[j * m_iWidth + i] - 1][1] += i * j;
            matrix[segmlabel[j * m_iWidth + i] - 1][2] += i;
            matrix[segmlabel[j * m_iWidth + i] - 1][3] += i * j;
            matrix[segmlabel[j * m_iWidth + i] - 1][4] += j * j;
            matrix[segmlabel[j * m_iWidth + i] - 1][5] += j;
            matrix[segmlabel[j * m_iWidth + i] - 1][6] += i;
            matrix[segmlabel[j * m_iWidth + i] - 1][7] += j;
            matrix[segmlabel[j * m_iWidth + i] - 1][8] += 1.0;

            Dmatrix[segmlabel[j * m_iWidth + i] - 1][0] += i * (labels[j * m_iWidth + i] * 1.0);
            Dmatrix[segmlabel[j * m_iWidth + i] - 1][1] += j * (labels[j * m_iWidth + i] * 1.0);
            Dmatrix[segmlabel[j * m_iWidth + i] - 1][2] += (labels[j * m_iWidth + i] * 1.0);
        }
    }

    // (1) Memory allocation for matrices
    src_matrix = cvCreateMat (nrow, ncol, CV_64FC1);
    dst_matrix = cvCreateMat (ncol, nrow, CV_64FC1);
    mul_matrix = cvCreateMat (nrow, nrow, CV_64FC1);
    src_vec    = cvCreateMat (nrow, 1, CV_64FC1);

    // (2) Substitute values into matrices
    for (int k = 0; k < num_segm; k++)
    {
            cvZero(src_matrix);
            cvZero(dst_matrix);
            cvZero(mul_matrix);
            cvZero(src_vec);

            for (i = 0; i < src_matrix->rows; i++) {
                for (j = 0; j < src_matrix->cols; j++) {
                    cvmSet (src_matrix, i, j, matrix[k][3*i + j]);
                }
            }
            // (3) Compute inverse matrix of src_matrix
            // Double type is not enough for computing inverse matrix.
            // So, long double is used.
            det = cvInvert (src_matrix, dst_matrix, CV_SVD);

            // (4) Display the determinant of matrix
//          printf ("det(src_matrix)=%f\n", det);

            // (5) Display inverse matrix
//          printf ("dst_matrix\n");
//          for (int ii = 0; ii < dst_matrix->rows; ii++)
//          {
//              for (int jj = 0; jj < dst_matrix->cols; jj++)
//              {
//                  printf ("%f \t", cvmGet (dst_matrix, ii, jj));
//              }
//              printf ("\n");
//          }

            // (6) Compute the product of matrix and its inverse matrix
            cvMatMul (src_matrix, dst_matrix, mul_matrix);
//          printf ("mul_matrix\n");
//          for (int ii = 0; ii < mul_matrix->rows; ii++)
//          {
//              for (int jj = 0; jj < mul_matrix->cols; jj++)
//              {
//                  printf ("%f \t", cvmGet (mul_matrix, ii, jj));
//              }
//              printf ("\n");
//          }

            //If mul_matrix is 3x3 identity matrix, a plane is fitted to the segment and depth values in the segment are recalculated.
            //Other wise, a plane is not fitted to the segment and depth values are not recalculated.
            if(fabs(fabs(mul_matrix->data.db[0]) - 1.0) <= 0.0001 && fabs(mul_matrix->data.db[1]) <= 0.0001 && fabs(mul_matrix->data.db[2]) <= 0.0001
                && fabs(mul_matrix->data.db[3]) <= 0.0001 && fabs(fabs(mul_matrix->data.db[4]) - 1.0) <= 0.0001 && fabs(mul_matrix->data.db[5]) <= 0.0001
                && fabs(mul_matrix->data.db[6]) <= 0.0001 && fabs(mul_matrix->data.db[7]) <= 0.0001 && fabs(fabs(mul_matrix->data.db[8]) - 1.0) <= 0.0001
            )
                reliability_matrix[k] = 1;
            else
                reliability_matrix[k] = 0;

            for (i = 0; i < src_matrix->rows; i++)
            {
                for (j = 0; j < src_matrix->cols; j++)
                {
                    inv_matrix[k][3*i + j] = cvmGet(dst_matrix, i, j);
                }
            }

            prod_matrix[k][0] = inv_matrix[k][0] * Dmatrix[k][0] + inv_matrix[k][1] * Dmatrix[k][1] + inv_matrix[k][2] * Dmatrix[k][2];
            prod_matrix[k][1] = inv_matrix[k][3] * Dmatrix[k][0] + inv_matrix[k][4] * Dmatrix[k][1] + inv_matrix[k][5] * Dmatrix[k][2];
            prod_matrix[k][2] = inv_matrix[k][6] * Dmatrix[k][0] + inv_matrix[k][7] * Dmatrix[k][1] + inv_matrix[k][8] * Dmatrix[k][2];
    } //for k


    int tmp_label;
    for(j = 0; j < m_iHeight; j++)
    {
        for (i = 0; i < m_iWidth; i++)
        {
            if(reliability_matrix[segmlabel[j * m_iWidth + i] - 1] == 1)
            {
                tmp_label = (int)(prod_matrix[segmlabel[j * m_iWidth + i] - 1][0] * i + prod_matrix[segmlabel[j * m_iWidth + i] - 1][1] * j + prod_matrix[segmlabel[j * m_iWidth + i] - 1][2] + 0.5);
                if(tmp_label < 0)
                    labels[j * m_iWidth + i] = (DepthType)(0);
                //else if (tmp_label>=(MAX_DEPTH-1))
                //    labels[j * m_iWidth + i] = (DepthType)(MAX_DEPTH-1);
                else if (tmp_label>=(m_iNumOfLabels-1))                        //KDDI R&D Labs., Inc. Bugfix m32570
                    labels[j * m_iWidth + i] = (DepthType)(m_iNumOfLabels-1);
                else
                    labels[j * m_iWidth + i] = (DepthType)(tmp_label);
            }
            else
            {
                tmp_label = (int)labels[j * m_iWidth + i];
                labels[j * m_iWidth + i] = (DepthType)(tmp_label);
            }
        }
    }
    cvReleaseMat (&src_matrix);
    cvReleaseMat (&dst_matrix);
    cvReleaseMat (&mul_matrix);
    cvReleaseMat (&src_vec);

    int pp;
    for(j=pp=0; j<m_iHeight; j++)
    {
        for(i=0; i<m_iWidth; i++,pp++)
        {
            // GIST start
//          labels_prev[pp] = labels[pp]; // now done in post-processing function
            // GIST end
            pDepth[j][i] = m_acLabel2Depth[labels[pp]];
        }
    }

    free(matrix[0]);
    free(matrix);
    free(Dmatrix[0]);
    free(Dmatrix);
    free(inv_matrix[0]);
    free(inv_matrix);
    free(prod_matrix[0]);
    free(prod_matrix);
    free(check_matrix[0]);
    free(check_matrix);

    free(reliability_matrix);

}
//Nagoya end
//Owieczka
void CEstimation::store_errors(CIYuv<ImageType> *yuvErrors,FILE *fpo_e)
{
    int i, j, d, pp;
    for(d=0; d<m_iNumOfLabels; d++)
    {
        for(j=pp=0; j<m_iHeight; j++)
        {
             for(i=0; i<m_iWidth; i++, pp++)
             {
#ifdef POZNAN_TWO_COSTS
                  yuvErrors->Y[j][i] = MINSEL(errors_left[d][pp],errors_right[d][pp]);
#else
                  yuvErrors->Y[j][i] = errors[d][pp];
#endif
             }
        }
        yuvErrors->writeOneFrame(fpo_e);
    }
}

#ifdef POZNAN_OCC
//Owieczka Occ Map Podzieli?przez precision
/*
void CEstimation::estimation_occ()
{
  int y,x,pp;
  
#ifdef POZNAN_OCC_VERBOSE
  IplImage *ImgDisparity = cvCreateImage(cvSize(m_iWidth,m_iHeight), IPL_DEPTH_8U, 1);

  //median(m_iHeight, m_iWidth, 3, pDepth);

  for(y=0,pp=0; y<m_iHeight; y++)   
  {
    for(x=0; x<m_iWidth; x++,pp++) 
    {
      ImgDisparity->imageData[y*ImgDisparity->widthStep+x] = m_aiLabel2Disparity[0][labels[pp]]; //Powinien by?drugi taki blok dla drugiego obrazu (mo�e po porstu odwo�ywa?si?do labels a nie do imgdisp
    }
  }
  
  cvSaveImage("disp.bmp",ImgDisparity);

  SAFE_RELEASE_IMAGE(ImgDisparity)
#endif

  int depth = 0;
  int dir = 1;
  int step = 0;
  int j;

  for(pp=0,y=0; y<m_iHeight; y++)   
  {
    depth = 0;
    x = (dir>0)?0:m_iWidth-1;
    pp += (dir>0)?0:m_iWidth-1;
    step = dir;
    while((x<m_iWidth)&&(x>=0))
    {
      disparity_map_left[pp] = MAX_DEPTH;
      DepthType disp = m_aiLabel2Disparity[0][labels[pp]];
      if(depth<disp)
      {
        depth = disp;
        if((disp!=0)&&(disp*dir/m_iPrecision+x>=0)&&(disp*dir/m_iPrecision+x<m_iWidth))
        {
          disparity_map_left[pp] = 0;
        }
      }
      disparity_map_left[pp] = depth;
      depth -= m_iPrecision;
      x += step;
      pp += step;
    }
    pp += (dir>0)?0:m_iWidth;
  }

  //cvSmooth(ipl_manoccl, ipl_manoccl, CV_MEDIAN, 3);
  //cvDilate(ipl_manoccl,ipl_manoccl,NULL,10);
#ifdef POZNAN_OCC_VERBOSE
  IplImage* ipl_manoccl = cvCreateImage(cvSize(m_iWidth,m_iHeight), IPL_DEPTH_8U, 1);
  
  for(y=0,pp=0; y<m_iHeight; y++)   
  {
    for(x=0; x<m_iWidth; x++,pp++) 
    {
      ipl_manoccl->imageData[y*ipl_manoccl->widthStep+x] = disparity_map_left[pp]; 
    }
  }
  
  cvSaveImage("occl.bmp",ipl_manoccl);

  SAFE_RELEASE_IMAGE(ipl_manoccl)
#endif

  dir = -1;

  for(pp=0,y=0; y<m_iHeight; y++)   
  {
    //printf("y=%d\n",y);
    depth = 0;
    x = (dir>0)?0:m_iWidth-1;
    pp += (dir>0)?0:m_iWidth-1;
    step = dir;
    while((x<m_iWidth)&&(x>=0))
    {
      disparity_map_right[pp] = MAX_DEPTH;
      DepthType disp = m_aiLabel2Disparity[1][labels[pp]];
      if(depth<disp)
      {
        depth = disp;
        if((disp!=0)&&(disp*dir/m_iPrecision+x>=0)&&(disp*dir/m_iPrecision+x<m_iWidth))
        {
          disparity_map_right[pp] = 0;
        }
      }
      disparity_map_right[pp] = depth;
      depth -= m_iPrecision;
      x += step;
      pp += step;
    }
    pp += (dir>0)?0:m_iWidth;
  }
  //cvSmooth(ipl_manoccr, ipl_manoccr, CV_MEDIAN, 3);
  //cvDilate(ipl_manoccr,ipl_manoccr,NULL,10);
#ifdef POZNAN_OCC_VERBOSE
  IplImage* ipl_manoccr = cvCreateImage(cvSize(m_iWidth,m_iHeight), IPL_DEPTH_8U, 1);
  
  for(y=0,pp=0; y<m_iHeight; y++)   
  {
    for(x=0; x<m_iWidth; x++,pp++) 
    {
      ipl_manoccr->imageData[y*ipl_manoccr->widthStep+x] = disparity_map_right[pp]; 
    }
  }
  
  cvSaveImage("occr.bmp",ipl_manoccr);

  SAFE_RELEASE_IMAGE(ipl_manoccr)
#endif

} //*/

inline void CEstimation::estimation_occ()
{
  (this->*func_estimation_occ)();
}

void CEstimation::estimation_occ_disparity()
{
  static int itr = 0;
  itr++;
	int y,x,pp;

  for(pp=0,y=0; y<m_iHeight; y++)   
	{
		for(x=0;x<m_iWidth;x++,pp++)
		{
      disparity_map_left[pp] = 0;
      disparity_map_right[pp] = 0;
    }
  }

	int depth = 0;
	int dir = 1;
	int step = 0;
	int j;
	int a;

	for(pp=0,y=0; y<m_iHeight; y++)   
	{
		for(x=0;x<m_iWidth;x++,pp++)
		{
      a = m_aiLabel2Disparity[0][labels[pp]]*dir/m_iPrecision+x;
			if(a>=0&&a<m_iWidth)
			{
				//ipl_manoccl->imageData[y*m_iWidth+a] = max(ipl_manoccl->imageData[y*m_iWidth+a],ImgDisparity->imageData[y*m_iWidth+x]);
#ifdef POZNAN_OCC_ALWAYS_DEPTH
				disparity_map_left[y*m_iWidth+a] = max(disparity_map_left[y*m_iWidth+a],(DepthType)m_acLabel2Depth[labels[pp]]);
#else
				disparity_map_left[y*m_iWidth+a] = max(disparity_map_left[y*m_iWidth+a],(DepthType)m_aiLabel2Disparity[0][labels[pp]]);
#endif
			}
		}
  }
	//cvSmooth(ipl_manoccl, ipl_manoccl, CV_MEDIAN, 3);
	//cvDilate(ipl_manoccl,ipl_manoccl,NULL,3);
	//cvErode(ipl_manoccl,ipl_manoccl,NULL,1);
	//cvDilate(ipl_manoccl,ipl_manoccl,NULL,1);

#ifdef POZNAN_OCC_VERBOSE
  CIYuv<DepthType> yuvOccLeft;
  if(yuvOccLeft.resize(m_iHeight, m_iWidth, POZNAN_DEPTHMAP_CHROMA_FORMAT))
  {
    for(y=0,pp=0; y<m_iHeight; y++)   
    {
      for(x=0; x<m_iWidth; x++,pp++) 
       {
        yuvOccLeft.Y[y][x] = disparity_map_left[pp]; 
      }
    }
    char name[255];
    static int frameno=0;
    FILE* fpo;
#ifdef POZNAN_16BIT_DEPTH
    sprintf(name,"occl%dx%d_16bps_cf%03d.yuv",m_iWidth,m_iHeight,POZNAN_DEPTHMAP_CHROMA_FORMAT);
#else
    sprintf(name,"occl%dx%d_cf%03d.yuv",m_iWidth,m_iHeight,POZNAN_DEPTHMAP_CHROMA_FORMAT);
#endif

    if(frameno==0)
    {
      if((fpo=fopen(name, "wb+"))!=NULL)
      {
        fclose(fpo);
      }
      frameno++;
    }    
    if((fpo=fopen(name, "ab"))!=NULL)
    {
      yuvOccLeft.writeOneFrame(fpo);
      fclose(fpo);
    }
  }
  /*
  IplImage* ipl_manoccl = cvCreateImage(cvSize(m_iWidth,m_iHeight), IPL_DEPTH_8U, 1);
  
  for(y=0,pp=0; y<m_iHeight; y++)   
  {
    for(x=0; x<m_iWidth; x++,pp++) 
    {
      ipl_manoccl->imageData[y*ipl_manoccl->widthStep+x] = disparity_map_left[pp]; 
    }
  }
  
  char name[255];
  sprintf(name,"occl%04d.bmp",itr);
  cvSaveImage(name,ipl_manoccl);

  {
  char a = 128;
  sprintf(name,"occl%dx%d.yuv",m_iWidth,m_iHeight);
  FILE* fi = fopen(name,"a");
  fwrite(ipl_manoccl->imageData,1,ipl_manoccl->widthStep*ipl_manoccl->height,fi);
  for(y=0;y<m_iHeight*m_iWidth/2;y++)
    fwrite(&a,1,1,fi);
  fclose(fi);
  }
  SAFE_RELEASE_IMAGE(ipl_manoccl)
  //*/
#endif


	dir = -1;
	
	for(pp=0,y=0; y<m_iHeight; y++)   
	{
		for(x=0;x<m_iWidth;x++,pp++)
		{
			//a = ImgDisparity->imageData[y*m_iWidth+x]*dir/m_iPrecision+x;
      a = m_aiLabel2Disparity[1][labels[pp]]*dir/m_iPrecision+x;
			if(a>=0&&a<m_iWidth)
			{
				//ipl_manoccr->imageData[y*m_iWidth+a] = max(ipl_manoccr->imageData[y*m_iWidth+a],ImgDisparity->imageData[y*m_iWidth+x]);
#ifdef POZNAN_OCC_ALWAYS_DEPTH
				disparity_map_right[y*m_iWidth+a] = max(disparity_map_right[y*m_iWidth+a],(DepthType)m_acLabel2Depth[labels[pp]]);
#else
				disparity_map_right[y*m_iWidth+a] = max(disparity_map_right[y*m_iWidth+a],(DepthType)m_aiLabel2Disparity[1][labels[pp]]);
#endif
			}
		}
  }
	//cvSmooth(ipl_manoccr, ipl_manoccr, CV_MEDIAN, 3);
	//cvDilate(ipl_manoccr,ipl_manoccr,NULL,3);
	//cvErode(ipl_manoccr,ipl_manoccr,NULL,1);
	//cvDilate(ipl_manoccr,ipl_manoccr,NULL,1);

#ifdef POZNAN_OCC_VERBOSE
  CIYuv<DepthType> yuvOccRight;
  if(yuvOccRight.resize(m_iHeight,m_iWidth,POZNAN_DEPTHMAP_CHROMA_FORMAT))
  {
    for(y=0,pp=0; y<m_iHeight; y++)   
    {
      for(x=0; x<m_iWidth; x++,pp++) 
       {
        yuvOccRight.Y[y][x] = disparity_map_right[pp]; 
      }
    }
    char name[255];
    static int frameno=0;
    FILE* fpo;
#ifdef POZNAN_16BIT_DEPTH
    sprintf(name,"occr%dx%d_16bps_cf%03d.yuv",m_iWidth,m_iHeight,POZNAN_DEPTHMAP_CHROMA_FORMAT);
#else
    sprintf(name,"occr%dx%d_cf%03d.yuv",m_iWidth,m_iHeight,POZNAN_DEPTHMAP_CHROMA_FORMAT);
#endif
    //if((fpo=fopen(name, "wb+"))!=NULL)
    if(frameno==0)
    {
      if((fpo=fopen(name, "wb+"))!=NULL)
      {
        fclose(fpo);
      }
      frameno++;
    }   
    if((fpo=fopen(name, "ab"))!=NULL)
    {
      yuvOccRight.writeOneFrame(fpo);
      fclose(fpo);
    }
  }
  /*
  IplImage* ipl_manoccr = cvCreateImage(cvSize(m_iWidth,m_iHeight), IPL_DEPTH_8U, 1);
  
  for(y=0,pp=0; y<m_iHeight; y++)   
  {
    for(x=0; x<m_iWidth; x++,pp++) 
    {
      ipl_manoccr->imageData[y*ipl_manoccr->widthStep+x] = disparity_map_right[pp]; 
    }
  }
  
  sprintf(name,"occr%04d.bmp",itr);
  cvSaveImage(name,ipl_manoccr);

  {
  char a = 128;
  sprintf(name,"occr%dx%d.yuv",m_iWidth,m_iHeight);
  FILE* fi = fopen(name,"a");
  fwrite(ipl_manoccr->imageData,1,ipl_manoccr->widthStep*ipl_manoccr->height,fi);
  for(y=0;y<m_iHeight*m_iWidth/2;y++)
    fwrite(&a,1,1,fi);
  fclose(fi);
  }
  SAFE_RELEASE_IMAGE(ipl_manoccr)
  //*/
#endif
} 
//*/

void CEstimation::estimation_occ_depth()
{
  static int itr = 0;
  itr++;
	int y,x,pp;
  CvMat* pix[2];
  pix[0] = cvCreateMat(3, 1, CV_64F);
  pix[1] = cvCreateMat(3, 1, CV_64F);
  cvmSet(pix[0], 2, 0, 1.0);

  for(pp=0,y=0; y<m_iHeight; y++)   
	{
		for(x=0;x<m_iWidth;x++,pp++)
		{
      disparity_map_left[pp] = 0;
      disparity_map_right[pp] = 0;
    }
  }

	int depth = 0;
	int dir = 1;
	int step = 0;
	int j;
	int target_pixel_u, target_pixel_v;

	for(pp=0,y=0; y<m_iHeight; y++)   
	{
    cvmSet(pix[0], 1, 0, y);
		for(x=0;x<m_iWidth;x++,pp++)
		{
      cvmSet(pix[0], 0, 0, x);
      cvmMul(m_matH_V2L[labels[pp]], pix[0], pix[1]);

      //target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
      target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
      target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);

      if(target_pixel_u >= 0 && target_pixel_u < m_iWidth && target_pixel_v >= 0 && target_pixel_v < m_iHeight)
      {
        //ipl_manoccl->imageData[y*m_iWidth+a] = max(ipl_manoccl->imageData[y*m_iWidth+a],ImgDisparity->imageData[y*m_iWidth+x]);
        disparity_map_left[target_pixel_v*m_iWidth+target_pixel_u] = max(disparity_map_left[target_pixel_v*m_iWidth+target_pixel_u],(DepthType)m_acLabel2Depth[labels[pp]]);
      }
    }
  }
  //cvSmooth(ipl_manoccl, ipl_manoccl, CV_MEDIAN, 3);
  //cvDilate(ipl_manoccl,ipl_manoccl,NULL,3);
  //cvErode(ipl_manoccl,ipl_manoccl,NULL,1);
  //cvDilate(ipl_manoccl,ipl_manoccl,NULL,1);

#ifdef POZNAN_OCC_VERBOSE
  CIYuv<DepthType> yuvOccLeft;
  if(yuvOccLeft.resize(m_iHeight, m_iWidth, POZNAN_DEPTHMAP_CHROMA_FORMAT))
  {
    for(y=0,pp=0; y<m_iHeight; y++)   
    {
      for(x=0; x<m_iWidth; x++,pp++) 
       {
        yuvOccLeft.Y[y][x] = disparity_map_left[pp]; 
      }
    }
    char name[255];
    static int frameno=0;
    FILE* fpo;
#ifdef POZNAN_16BIT_DEPTH
    sprintf(name,"occl%dx%d_16bps_cf%03d.yuv",m_iWidth,m_iHeight,POZNAN_DEPTHMAP_CHROMA_FORMAT);
#else
    sprintf(name,"occl%dx%d_cf%03d.yuv",m_iWidth,m_iHeight,POZNAN_DEPTHMAP_CHROMA_FORMAT);
#endif

    if(frameno==0)
    {
      if((fpo=fopen(name, "wb+"))!=NULL)
      {
        fclose(fpo);
      }
      frameno++;
    }    
    if((fpo=fopen(name, "ab"))!=NULL)
    {
      yuvOccLeft.writeOneFrame(fpo);
      fclose(fpo);
    }
  }
  /*
  IplImage* ipl_manoccl = cvCreateImage(cvSize(m_iWidth,m_iHeight), IPL_DEPTH_8U, 1);
  
  for(y=0,pp=0; y<m_iHeight; y++)   
  {
    for(x=0; x<m_iWidth; x++,pp++) 
    {
      ipl_manoccl->imageData[y*ipl_manoccl->widthStep+x] = disparity_map_left[pp]; 
    }
  }
  
  char name[255];
  sprintf(name,"occl%04d.bmp",itr);
  cvSaveImage(name,ipl_manoccl);
  SAFE_RELEASE_IMAGE(ipl_manoccl)
  //*/
#endif


	dir = -1;
	
	for(pp=0,y=0; y<m_iHeight; y++)   
	{
    cvmSet(pix[0], 1, 0, y);
		for(x=0;x<m_iWidth;x++,pp++)
		{
      cvmSet(pix[0], 0, 0, x);
      cvmMul(m_matH_V2R[labels[pp]], pix[0], pix[1]);
      //target_pixel_u = (int)(cvmGet(pix[1], 0, 0)*m_iPrecision/cvmGet(pix[1], 2, 0) + 0.5);
      target_pixel_u = (int)(cvmGet(pix[1], 0, 0)/cvmGet(pix[1], 2, 0) + 0.5);
      target_pixel_v = (int)(cvmGet(pix[1], 1, 0)/cvmGet(pix[1], 2, 0) + 0.5);

      if(target_pixel_u >= 0 && target_pixel_u < m_iWidth && target_pixel_v >= 0 && target_pixel_v < m_iHeight)
      {
        //ipl_manoccr->imageData[y*m_iWidth+a] = max(ipl_manoccr->imageData[y*m_iWidth+a],ImgDisparity->imageData[y*m_iWidth+x]);
        disparity_map_right[target_pixel_v*m_iWidth+target_pixel_u] = max(disparity_map_right[target_pixel_v*m_iWidth+target_pixel_u],(DepthType)m_acLabel2Depth[labels[pp]]);
      }
    }
  }
  //cvSmooth(ipl_manoccr, ipl_manoccr, CV_MEDIAN, 3);
  //cvDilate(ipl_manoccr,ipl_manoccr,NULL,3);
  //cvErode(ipl_manoccr,ipl_manoccr,NULL,1);
  //cvDilate(ipl_manoccr,ipl_manoccr,NULL,1);

#ifdef POZNAN_OCC_VERBOSE
  CIYuv<DepthType> yuvOccRight;
  if(yuvOccRight.resize(m_iHeight,m_iWidth,POZNAN_DEPTHMAP_CHROMA_FORMAT))
  {
    for(y=0,pp=0; y<m_iHeight; y++)   
    {
      for(x=0; x<m_iWidth; x++,pp++) 
       {
        yuvOccRight.Y[y][x] = disparity_map_right[pp]; 
      }
    }
    char name[255];
    static int frameno=0;
    FILE* fpo;
#ifdef POZNAN_16BIT_DEPTH
    sprintf(name,"occr%dx%d_16bps_cf%03d.yuv",m_iWidth,m_iHeight,POZNAN_DEPTHMAP_CHROMA_FORMAT);
#else
    sprintf(name,"occr%dx%d_cf%03d.yuv",m_iWidth,m_iHeight,POZNAN_DEPTHMAP_CHROMA_FORMAT);
#endif
    //if((fpo=fopen(name, "wb+"))!=NULL)
    if(frameno==0)
    {
      if((fpo=fopen(name, "wb+"))!=NULL)
      {
        fclose(fpo);
      }
      frameno++;
    }   
    if((fpo=fopen(name, "ab"))!=NULL)
    {
      yuvOccRight.writeOneFrame(fpo);
      fclose(fpo);
    }
  }

  /*
  IplImage* ipl_manoccr = cvCreateImage(cvSize(m_iWidth,m_iHeight), IPL_DEPTH_8U, 1);
  
  for(y=0,pp=0; y<m_iHeight; y++)   
  {
    for(x=0; x<m_iWidth; x++,pp++) 
    {
      ipl_manoccr->imageData[y*ipl_manoccr->widthStep+x] = disparity_map_right[pp]; 
    }
  }
  
  sprintf(name,"occr%04d.bmp",itr);
  cvSaveImage(name,ipl_manoccr);
  SAFE_RELEASE_IMAGE(ipl_manoccr)
  //*/
#endif
  cvReleaseMat(&pix[0]);
  cvReleaseMat(&pix[1]);
} 
//*/

#endif

#ifdef POZNAN_TWOVIEW_SUPPORT
void CEstimation::clear_error(CostType **error)
{
  int i, j, d, pp;

  for(j=pp=0; j<m_iHeight; j++)
  {
    for(i=0; i<m_iWidth; i++, pp++)
    {
      for(d=0; d<m_iNumOfLabels; d++)
      {
        error[pp][d] = COST_MAX;
      }
    }
  }
}
#endif