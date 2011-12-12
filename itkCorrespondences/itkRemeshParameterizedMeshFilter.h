#ifndef _itkRemeshParameterizedMeshFilter_h
#define _itkRemeshParameterizedMeshFilter_h

#include "itkMeshToMeshFilter.h"


namespace itk
{
  /** \class RemeshParameterizedMeshFilter
  *   \brief Class for mapping a mesh in parameter space into object space.
  * 
  * Basically a remeshing of the parameterized input mesh with vertices at positions indicated by the 
  * specified landmarks.
  * All point data from TParameterization is interpolated to TOutputMesh, so the used point types of
  * these template parameters have to be the same!
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, German Cancer Research Center, Heidelberg, Germany.
  */
  template <class TParameterization, class TLandmarkMesh, class TOutputMesh>
  class RemeshParameterizedMeshFilter : public MeshToMeshFilter
    <TParameterization, TOutputMesh>
  {

  public:

    /** Standard typedefs. */
    typedef RemeshParameterizedMeshFilter              Self;
    typedef MeshToMeshFilter<TParameterization, TOutputMesh> Superclass;
    typedef SmartPointer<Self>                        Pointer;
    typedef SmartPointer<const Self>                  ConstPointer;

    /** Convenient typedefs. */
    typedef TLandmarkMesh                               LandmarkMeshType;
    typedef TParameterization                           ParameterizationType;
    typedef TOutputMesh                                 OutputMeshType;
    typedef typename LandmarkMeshType::Pointer          LandmarkPointer;
    typedef typename LandmarkMeshType::PointType        LandmarkPointType;
    typedef typename ParameterizationType::PointType    ParameterizationPointType;
    typedef typename ParameterizationType::Pointer      ParameterizationPointer;
    typedef typename ParameterizationType::IndexType    IndexType;
    typedef typename ParameterizationType::IndexedTriangleMeshType MeshType;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(RemeshParameterizedMeshFilter, MeshToMeshFilter);

    void SetLandmarks( LandmarkMeshType *landmarks );

    /** Copy all cache information from another RemeshParameterizedMeshFilter.
    * This makes sense if the input is a reparameterized version of another mesh
    * which has its own remeshing filter with the same landmarks.
    * If the landmarks are different, the method does not do anything.
    */
    void CloneCache( Pointer source );
    
    void SetPointDataDimension (int n) ;
    
  protected:

    RemeshParameterizedMeshFilter();

    ~RemeshParameterizedMeshFilter();

    void PrintSelf( std::ostream& os, Indent indent ) const;

    virtual void GenerateData();

    /** Copies the cell information from landmarks to output and initializes the required point containers. */
    void InitializeOutput( LandmarkMeshType *landmarks );

    /** Sets one point of the output mesh. */
    virtual void CalculateOutputValues( unsigned int pointId, IndexType faceId, double barycentricP, double barycentricQ );

    
  private:

    ParameterizationPointer m_Parameterization;
    LandmarkPointer         m_Landmarks;
    typename ParameterizationType::PointsContainer    *m_InPoints;
    typename ParameterizationType::PointDataContainer *m_InPointData;
    typename OutputMeshType::PointsContainer          *m_OutPoints;
    typename OutputMeshType::PointDataContainer       *m_OutPointData;
    std::vector<IndexType>  m_LastFaceHit;
    bool                    m_NewLandmarksSet;
    bool                    m_FaceHitsValid;
    bool                    m_CalculatePointData;
    
    int                     m_PointDataDimension ;
  };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRemeshParameterizedMeshFilter.txx"
#endif

#endif
