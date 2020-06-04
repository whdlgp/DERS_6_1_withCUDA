//======================================created by Ying Chen =====================================
//===============Tampere University of Technology (TUT)/Nokia Research Center (NRC)===============

#if !defined AFX_PARAMETER_DEPTH_ESTIMATION_H
#define AFX_PARAMETER_DEPTH_ESTIMATION_H

#include "iostream"
#include "string"
#include "version.h"
#include "cstring"

#include "ParameterBase.h"

using namespace std;


class CParameterDepthEstimation : public ParameterBase
{

public:
    CParameterDepthEstimation               ();
    virtual ~CParameterDepthEstimation      ();

    Int                 xInit                       ( Int argc, Char** argv);

    Void                setDepthType                ( UInt ui )         { m_uiDepthType           = ui; }
    Void                setSourceWidth              ( UInt ui )         { m_uiSourceWidth         = ui; }
    Void                setSourceHeight             ( UInt ui )         { m_uiSourceHeight        = ui; }
    Void                setStartFrame               ( UInt ui )         { m_uiStartFrame          = ui; }
    Void                setNumberOfFrames           ( UInt ui )         { m_uiNumberOfFrames      = ui; }

    Void                setLeftCameraName           ( std::string s )   { m_cLeftCameraName       = s;  }
    Void                setCenterCameraName         ( std::string s )   { m_cm_uiCenterCameraName = s;  }
    Void                setRightCameraName          ( std::string s )   { m_cRightCameraName      = s;  }
    // ETRI start
    Void                setThresOfDepthDiff         ( Int i )           { m_iThresOfDepthDiff     = i;  }
    Void                setMovingObjectsBSize       ( Int i )           { m_iMovingObjectsBSize   = i;  }
    Void                setMotionSearchBSize        ( Int i )           { m_iMotionSearchBSize    = i;  }
    // ETRI end
    Void                setMinValDisSearchRange     ( Int i )           { m_iMinValDisSearchRange = i;  }
    Void                setMaxValDisSearchRange     ( Int i )           { m_iMaxValDisSearchRange = i;  }
    Void                setMinValDisRange           ( Int i )           { m_iMinValDisRange       = i;  }
    Void                setMaxValDisRange           ( Int i )           { m_iMaxValDisRange       = i;  }
    Void                setSmoothCoeff              ( Double d)         { m_dSmoothCoeff          = d;  }
    Void                setSmoothCoeff2             ( Double d)         { m_dSmoothCoeff2         = d;  } //Nagoya
    Void                setFileLeftView             ( std::string s )   { m_cFileLeftView         = s;  }
    Void                setFileCenterView           ( std::string s )   { m_cFileCenterView       = s;  }
    Void                setFileRightView            ( std::string s )   { m_cFileRightView        = s;  }
    Void                setFileOutputDepth          ( std::string s )   { m_cFileOutputDepth      = s;  }
    Void                setFileCamPara              ( std::string s )   { m_cFileCamPara          = s;  }

    Void                setBaselineBasis            ( UInt ui )         { m_uiBaselineBasis       = ui; }

#ifdef SUB_PEL_PRECISION
    Void                setPrecision                ( UInt ui )         { m_uiPrecision           = ui; }
    Void                setSearchLevel              ( Double d )        { m_dSearchLevel          = d; }
    Void                setFilter                   ( UInt ui )         { m_uiFilter              = ui; }
#endif
#ifdef ACCURATE_CORRESPONDENCE
    Void                setMatchingMethod           ( UInt ui )         { m_uiMatchingMethod      = ui; }
    //Nagoya start
    Void                setMatchingBlock            ( UInt ui )         { m_uiMatchingBlock         = ui; }
    Void                setImageSegmentation        ( UInt ui )         { m_uiImageSegmentation     = ui; }
    Void                setSegmentationmethod       ( UInt ui )         { m_uiSegmentationMethod    = ui; }
    Void                setMaxCluster               ( UInt ui )         { m_uiMaxCluster    = ui; }
    //Nagoya end
#endif
    // GIST start
    Void                setTemporalEnhancement      ( Int i   )         { m_iTemporalEnhancement = i; }
    Void                setThreshold                ( Double d)         { m_dThreshold           = d; }
    // GIST end

    //Poznan start
    Void                setSoftDistanceCoeff        ( Double d)         { m_dSoftDistanceCoeff    = d;  } //Owieczka
    Void                setSoftColorCoeff           ( Double d)         { m_dSoftColorCoeff       = d;  }
    Void                setSoftBlockWidth           ( UInt ui )         { m_uiSoftBlockWidth      = ui; }
    Void                setSoftBlockHeight          ( UInt ui )         { m_uiSoftBlockHeight     = ui; }
    Void                setFileOutputErrors         ( std::string s )   { m_cFileOutputErrors     = s;  }
    //Poznan end

    UInt                getDepthType              ()         { return m_uiDepthType; }
    UInt                getSourceWidth            ()         { return m_uiSourceWidth; }
    UInt                getSourceHeight           ()         { return m_uiSourceHeight; }
    UInt                getNumberOfFrames         ()         { return m_uiNumberOfFrames; }
    UInt                getStartFrame             ()         { return m_uiStartFrame; }
    const std::string   getLeftCameraName         ()         { return m_cLeftCameraName; }
    const std::string   getCenterCameraName       ()         { return m_cm_uiCenterCameraName; }
    const std::string   getRightCameraName        ()         { return m_cRightCameraName; }
    // ETRI start
    Int                 getThresOfDepthDiff       ()         { return m_iThresOfDepthDiff; }
    Int                 getMovingObjectsBSize     ()         { return m_iMovingObjectsBSize; }
    Int                 getMotionSearchBSize      ()         { return m_iMotionSearchBSize; }
    // ETRI end
    Int                 getMinValDisSearchRange   ()         { return m_iMinValDisSearchRange; }
    Int                 getMaxValDisSearchRange   ()         { return m_iMaxValDisSearchRange; }
    Int                 getMinValDisRange         ()         { return m_iMinValDisRange; }
    Int                 getMaxValDisRange         ()         { return m_iMaxValDisRange; }
#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
    Double              getNearestDepthValue      ()         { return m_dNearestDepthValue; };
    Double              getFarthestDepthValue     ()         { return m_dFarthestDepthValue; };
    Double              getNearestSearchDepthValue      ()   { return m_dNearestSearchDepthValue; };
    Double              getFarthestSearchDepthValue     ()   { return m_dFarthestSearchDepthValue; };
    Int                 getNumberOfDepthSteps     ()         { return m_iNumberOfDepthSteps; };
#endif
    Double              getSmoothCoeff            ()         { return m_dSmoothCoeff; }
    Double              getSmoothCoeff2           ()         { return m_dSmoothCoeff2; } //Nagoya
    const std::string   getFileLeftView           ()         { return m_cFileLeftView; }
    const std::string   getFileCenterView         ()         { return m_cFileCenterView; }
    const std::string   getFileRightView          ()         { return m_cFileRightView; }
    const std::string   getFileOutputDepth        ()         { return m_cFileOutputDepth; }
    const std::string   getFileCamPara            ()         { return m_cFileCamPara; }

    UInt                getBaselineBasis          ()         { return m_uiBaselineBasis; }

#ifdef SUB_PEL_PRECISION
    UInt                getPrecision              ()         { return m_uiPrecision; }
    Double              getSearchLevel            ()         { return m_dSearchLevel; }
    UInt                getFilter                 ()         { return m_uiFilter; }
#ifdef SUB_PEL_VERTICAL_PRECISION
    UInt                getVerticalPrecision      ()         { return m_uiVerticalPrecision; }
    UInt                getVerticalFilter         ()         { return m_uiVerticalFilter; }
#endif
#endif
#ifdef ACCURATE_CORRESPONDENCE
    UInt                getMatchingMethod         ()         { return m_uiMatchingMethod; }
    //Nagoya start
    UInt                getMatchingBlock          ()          {return m_uiMatchingBlock; }
    UInt                getImageSegmentation      ()          {return m_uiImageSegmentation; }
    UInt                getSegmentationMethod     ()          {return m_uiSegmentationMethod; }
    UInt                getMaxCluster             ()          {return m_uiMaxCluster; }
    //Nagoya end
#endif

    //Poznan start
    Double              getSoftDistanceCoeff      ()         { return m_dSoftDistanceCoeff; } //Owieczka
    Double              getSoftColorCoeff         ()         { return m_dSoftColorCoeff; }
    UInt                getSoftBlockWidth         ()         { return m_uiSoftBlockWidth; }
    UInt                getSoftBlockHeight        ()         { return m_uiSoftBlockHeight; }
    const std::string   getFileOutputErrors       ()         { return m_cFileOutputErrors; }
#ifdef POZNAN_OCC
    UInt                getOcc                    ()         { return m_iOcc; }
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
    UInt                getDirections             ()         { return m_iDirections; }
#endif
    //Poznan end

    // GIST start
    UInt                getTemporalEnhancement    ()         { return m_iTemporalEnhancement; }
    Double              getThreshold              ()         { return m_dThreshold; }
    // GIST end

    //Nagoya start
    UInt                getDEmode                 ()         {return m_uiDEmode; }
    const std::string   getFileCenterManual       ()         { return m_cFileManual; }
    Double              getTemporalWeight         ()         { return m_dTemporalWeight; }

    const std::string   getRefDepthCameraName     ()         { return m_cRefDepthCameraName; }
    const std::string   getRefDepthFile           ()         { return m_cRefDepthFile; }
    //Nagoya end

//  Void                setResult           ( Int     iResult )   { m_iResult = iResult;  }

    Bool                equals( const Char* str1, const Char* str2, UInt nLetter ) { return 0 == ::strncmp( str1, str2, nLetter); }

private:
    UInt                                setup                               ();

protected:
    Int             xPrintUsage         ( Char**  argv );

protected:
    UInt            m_uiDepthType  ;
    UInt            m_uiSourceWidth ;
    UInt            m_uiSourceHeight  ;
    UInt            m_uiNumberOfFrames  ;
    UInt            m_uiStartFrame; // DT

    std::string     m_cLeftCameraName ;
    std::string     m_cm_uiCenterCameraName ;
    std::string     m_cRightCameraName  ;
    // ETRI start
    Int             m_iThresOfDepthDiff;
    Int             m_iMovingObjectsBSize;
    Int             m_iMotionSearchBSize;
    // ETRI end
    Int             m_iMinValDisSearchRange ;
    Int             m_iMaxValDisSearchRange ;
    Int             m_iMinValDisRange ;
    Int             m_iMaxValDisRange ;
#ifdef POZNAN_ZNEAR_ZFAR_SEARCH_RANGE
    Double          m_dNearestDepthValue;
    Double          m_dFarthestDepthValue;
    Double          m_dNearestSearchDepthValue;
    Double          m_dFarthestSearchDepthValue;
    Int             m_iNumberOfDepthSteps;
#endif
    Double          m_dSmoothCoeff  ;
    Double          m_dSmoothCoeff2 ; //Nagoya
    std::string     m_cFileLeftView ;
    std::string     m_cFileCenterView ;
    std::string     m_cFileRightView  ;
    std::string     m_cFileOutputDepth  ;
    std::string     m_cFileCamPara  ;

    UInt            m_uiBaselineBasis ;

#ifdef SUB_PEL_PRECISION
    UInt            m_uiPrecision ;
    Double          m_dSearchLevel ;
    UInt            m_uiFilter;
#ifdef SUB_PEL_VERTICAL_PRECISION
    UInt            m_uiVerticalPrecision;
    UInt            m_uiVerticalFilter;
#endif
#endif
#ifdef ACCURATE_CORRESPONDENCE
    UInt            m_uiMatchingMethod ;
    //Nagoya start
    UInt            m_uiMatchingBlock;
    UInt            m_uiImageSegmentation;
    UInt            m_uiSegmentationMethod;
    UInt            m_uiMaxCluster;
    //Nagoya end
#endif

    //Poznan start - Soft Segment Owieczka
    Double          m_dSoftDistanceCoeff;
    Double          m_dSoftColorCoeff;
    UInt            m_uiSoftBlockWidth;
    UInt            m_uiSoftBlockHeight;
    std::string     m_cFileOutputErrors;
#ifdef POZNAN_OCC
    UInt            m_iOcc;
#endif
#ifdef POZNAN_TWOVIEW_SUPPORT
    UInt            m_iDirections;
#endif
    //Poznan end

    // GIST start
    Int             m_iTemporalEnhancement;
    Double          m_dThreshold;
    // GIST end

    //Nagoya start
    std::string     m_cFileManual;
    UInt            m_uiDEmode;
    Double          m_dTemporalWeight;
    std::string     m_cRefDepthCameraName;
    std::string     m_cRefDepthFile;
    //Nagoya end

};

#endif
