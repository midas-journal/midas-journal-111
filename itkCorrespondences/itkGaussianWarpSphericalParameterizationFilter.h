#ifndef _itkGaussianWarpSphericalParameterizationFilter_h
#define _itkGaussianWarpSphericalParameterizationFilter_h

#include "itkMeshToMeshFilter.h"
#include "itkSphericalParameterizedTriangleMesh.h"


namespace itk
{

  /** \class GaussianWarpSphericalParameterizationFilter
  *   \brief This class locally modifies the parameterization of an 
  *          itk::SphericalParameterizedTriangleMesh. 
  *
  * Within a local environment around a control point, all points of the 
  * spherical map (the parameterization) are moved in a specific direction. 
  * The magnitude of the movement is controlled by a Gaussian envelope 
  * function, its variance determines the area that is modified. 
  * The filter offers presets for different variances that can be activated 
  * by the SetLevelOfDetail() method. The number of available control points 
  * also depends on the level of detail (since fewer local environments fit 
  * on the sphere if each one is larger) and can be queried by 
  * GetNumberOfControlPoints().
  * When the level of detail is set, the control point that is used for the 
  * warp has to be specified using SetActiveControlPoint(). Subsequently, 
  * the direction and maximum magnitude of the movement has to be set using 
  * SetDirection(). 
  * Note that the filter only modifies the area around one control point at 
  * a time.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
            German Cancer Research Center, Heidelberg, Germany.
  */
  template <class TParameterization>
  class GaussianWarpSphericalParameterizationFilter : public MeshToMeshFilter<TParameterization, TParameterization>
  {

  public:

    /** Standard class typedefs. */
    typedef GaussianWarpSphericalParameterizationFilter   Self;
    typedef MeshToMeshFilter
      <TParameterization, TParameterization>    Superclass;
    typedef SmartPointer<Self>                  Pointer;
    typedef SmartPointer<const Self>            ConstPointer;

    /** Convenient typedefs. */
    typedef TParameterization                               ParameterizationType;
    typedef typename ParameterizationType::Pointer          ParameterizationPointer;
    typedef typename ParameterizationType::SphericalMapType SphericalMapType;
    typedef typename ParameterizationType::PointType        PointType;
    typedef typename ParameterizationType::VectorType       VectorType;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianWarpSphericalParameterizationFilter, MeshToMeshFilter);

    void SetInput( ParameterizationType *input )
    {
      m_NewInputSet = true;
      Superclass::SetInput( input );
    }

    /** Sets the level of detail for the modification.
    * Valid values range from 0 up to (including) the maximum level of detail.
    */
    void SetLevelOfDetail( unsigned int level );

    /** Gets the current level of detail. */
    itkGetConstMacro( LevelOfDetail, unsigned int );

    /** Return the maximum possible level of detail. */
    unsigned int GetMaximumLevelOfDetail()
    {
      return 3;
    }

    /** Gets the current number of control points 
    * (which depends on the level of detail.) 
    */
    itkGetConstMacro( NumberOfControlPoints, unsigned int );

    /** Sets the control point that will be used for the warp. */
    itkSetMacro( ActiveControlPoint, unsigned int );

    /** Sets the direction for the warp. 
    * The direction is specified by a three-dimensional vector: 
    * The first element is the change in radians for the azimuth angle, 
    * the second is the change in radians for the polar angle. 
    * The third element is not used.
    */
    itkSetMacro( Direction, VectorType );
    

  protected:

    GaussianWarpSphericalParameterizationFilter();

    ~GaussianWarpSphericalParameterizationFilter();

    void PrintSelf( std::ostream& os, Indent indent ) const;

    virtual void GenerateData();

    /** Copies the cell information from landmarks to output and initializes the required point containers. */
    void InitializeOutput( ParameterizationType *input );


  private:

    bool                    m_NewInputSet;
    unsigned int            m_LevelOfDetail;
    unsigned int            m_NumberOfControlPoints;
    unsigned int            m_ActiveControlPoint;
    std::vector<PointType>  m_ControlPoints;
    double                  m_KernelStandardDeviation;
    VectorType              m_Direction;

  };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianWarpSphericalParameterizationFilter.txx"
#endif

#endif

