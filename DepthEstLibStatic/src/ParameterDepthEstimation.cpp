//======================================created by Ying Chen =====================================
//===============Tampere University of Technology (TUT)/Nokia Research Center (NRC)===============

#include "ParameterDepthEstimation.h"

#include <string>

#ifndef MSYS_WIN32
#define stricmp strcasecmp
#endif

#define equal(a,b)  (!stricmp((a),(b)))


using namespace std;

CParameterDepthEstimation::CParameterDepthEstimation()
: m_uiDepthType           ( 0     )
, m_uiSourceWidth         ( 0     )
, m_uiSourceHeight        ( 0     )
, m_uiNumberOfFrames      ( 0     )
, m_uiStartFrame          ( 0     ) //DT
, m_cLeftCameraName       (       )
, m_cm_uiCenterCameraName (       )
, m_cRightCameraName      (       )
// ETRI start
, m_iThresOfDepthDiff     (       )
, m_iMovingObjectsBSize   (       )
, m_iMotionSearchBSize    (       )
// ETRI end
, m_iMinValDisSearchRange (       )
, m_iMaxValDisSearchRange (       )
, m_iMinValDisRange       (       )
, m_iMaxValDisRange       (       )
, m_dSmoothCoeff          ( 0.0   )
, m_dSmoothCoeff2         ( 0.0   ) //Nagoya
, m_cFileLeftView         (       )
, m_cFileCenterView       (       )
, m_cFileRightView        (       )
, m_cFileOutputDepth      (       )
, m_cFileCamPara          (       )
// GIST start
, m_iTemporalEnhancement  (       )
, m_dThreshold            (       )
// GIST end
// Poznan start
, m_dSoftDistanceCoeff    ( 0.0   ) //Owieczka
, m_dSoftColorCoeff       ( 0.0   )
, m_uiSoftBlockWidth      ( 0     )
, m_uiSoftBlockHeight     ( 0     )
, m_cFileOutputErrors     (       )
#ifdef POZNAN_OCC
, m_iOcc                  ( 0     )
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
, m_iDirections           ( 3     )  
#endif
// Poznan end
{
}



CParameterDepthEstimation::~CParameterDepthEstimation()
{

    release();
}




Int
CParameterDepthEstimation::xInit( Int argc, Char**  argv )
{
    if ( argc < 2 ) return -1;

    std::string cFilename = argv[1] ;

    setup();

    if(xReadFromFile( cFilename )!=1) return -1;
    if(xReadFromCommandLine(argc, argv)!=1) return -1;
#if 1
    xPrintParam();
#endif
    release();

    return 1;
}


Int
CParameterDepthEstimation::xPrintUsage( Char **argv )
{
  printf("\n supported options:\n\n");
  printf(" Parameter File Name\n\n");
  printf("\n");
  return 1;
}

UInt
CParameterDepthEstimation::setup()
{
    UInt uiParLnCount=0;

    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("DepthType",                             & m_uiDepthType,                   0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("SourceWidth"   ,                        & m_uiSourceWidth ,                1280 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("SourceHeight",                          & m_uiSourceHeight,                960 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("TotalNumberOfFrames",                   & m_uiNumberOfFrames,              300 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("StartFrame",                            & m_uiStartFrame,                  0 ); // DT
    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("LeftCameraName",                        & m_cLeftCameraName,               "param_dog36" );
    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("CenterCameraName",                      & m_cm_uiCenterCameraName,         "param_dog38" );
    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("RightCameraName",                       & m_cRightCameraName,              "param_dog40" );
    // ETRI start
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("ThresholdOfDepthDifference",            & m_iThresOfDepthDiff,             10 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("MovingObjectsBSize",                    & m_iMovingObjectsBSize,           0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("MotionSearchBSize",                     & m_iMotionSearchBSize,            0 );
    // ETRI end
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("MinimumValueOfDisparitySearchRange",    & m_iMinValDisSearchRange,         0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("MaximumValueOfDisparitySearchRange",    & m_iMaxValDisSearchRange,         20 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("MinimumValueOfDisparityRange",          & m_iMinValDisRange,               0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("MaximumValueOfDisparityRange",          & m_iMaxValDisRange,               20 );
#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("NearestDepthValue",                     & m_dNearestDepthValue,            -1.0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("FarthestDepthValue",                    & m_dFarthestDepthValue,           -10.0);
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("NearestSearchDepthValue",               & m_dNearestSearchDepthValue,      -1.0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("FarthestSearchDepthValue",              & m_dFarthestSearchDepthValue,     -10.0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("NumberOfDepthSteps",                    & m_iNumberOfDepthSteps,           10 );
#endif
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("SmoothingCoefficient",                  & m_dSmoothCoeff,                  4.0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("SmoothingCoefficient2",                 & m_dSmoothCoeff2,                 1.0 ); //Nagoya
    //m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileLeftViewImage",                     & m_cFileLeftView,                 "dog036.yuv" );
    //m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileCenterViewImage",                   & m_cFileCenterView,               "dog038.yuv");
	//m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileRightViewImage",                    & m_cFileRightView,                "dog040.yuv" );
	m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileLeftViewImage",                     & m_cFileLeftView,                 "" ); // no default
	m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileCenterViewImage",                   & m_cFileCenterView,               "");  // no default
	m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileRightViewImage",                    & m_cFileRightView,                "" ); // no default 

    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileOutputDepthMapImage",               & m_cFileOutputDepth,              "depth_dog038.yuv" );
    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileCameraParameter",                   & m_cFileCamPara,                  "cam_param_dog.txt" );

    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("BaselineBasis"   ,                      & m_uiBaselineBasis ,              1 );

#ifdef SUB_PEL_PRECISION
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("Precision"   ,                          & m_uiPrecision ,                  4 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("SearchLevel"   ,                        & m_dSearchLevel ,                 2 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("Filter"   ,                             & m_uiFilter ,                     1 );
#ifdef SUB_PEL_VERTICAL_PRECISION
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("VerticalPrecision"   ,                  & m_uiVerticalPrecision ,          1 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("VerticalFilter"   ,                     & m_uiVerticalFilter ,             1 );
#endif
#endif
#ifdef ACCURATE_CORRESPONDENCE
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("MatchingMethod"   ,                     & m_uiMatchingMethod ,             2 );
    //Nagoya start
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("MatchingBlock"   ,                      & m_uiMatchingBlock ,             1 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("ImageSegmentation"   ,                  & m_uiImageSegmentation ,         0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("SegmentationMethod"   ,                 & m_uiSegmentationMethod ,        3 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("MaxCluster"   ,                         & m_uiMaxCluster ,               32 );
    //Nagoya end
#endif

    // GIST start
    m_pCfgLines[uiParLnCount++] = new ConfigLineInt ("TemporalEnhancement"   ,                & m_iTemporalEnhancement ,         0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("Threshold"   ,                          & m_dThreshold ,                   0.00 );
    // GIST end

    //Nagoya Semi-automatic start
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("DepthEstimationMode",                   & m_uiDEmode,                      0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileCenterManual",                      & m_cFileManual,                   "dog038.yuv");
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("TemporalWeight",                        & m_dTemporalWeight,               0.8 );

    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("RefDepthCameraName",                    & m_cRefDepthCameraName,           "CamParam36");
    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("RefDepthFile",                          & m_cRefDepthFile,                 "ref_depth.yuv");
    //Nagoya Semi-automatic end
    //Poznan start - Owieczka
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("SoftDistanceCoeff"   ,                  & m_dSoftDistanceCoeff ,           20.0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineDbl ("SoftColorCoeff"   ,                     & m_dSoftColorCoeff ,              20.0 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("SoftBlockWidth"   ,                     & m_uiSoftBlockWidth ,             3 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("SoftBlockHeight"  ,                     & m_uiSoftBlockHeight ,            3 );
    m_pCfgLines[uiParLnCount++] = new ConfigLineStr ("FileOutputErrors" ,                     & m_cFileOutputErrors , "errors.yuv" );    
#ifdef POZNAN_OCC
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("Occlusion" ,                            & m_iOcc ,                         0 );    
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
    m_pCfgLines[uiParLnCount++] = new ConfigLineUInt("MatchDirection" ,                            & m_iDirections ,                  3 );    
#endif    
    //Poznan end
    m_pCfgLines[uiParLnCount] = NULL;

    return uiParLnCount;
}
