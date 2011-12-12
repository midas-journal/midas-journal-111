#ifndef _itkRotateSphericalParameterizationFilter_h
#define _itkRotateSphericalParameterizationFilter_h

#include "itkMeshToMeshFilter.h"
#include "itkSphericalParameterizedTriangleMesh.h"
#include "itkAffineTransform.h"


namespace itk
{

  /** \class RotateSphericalParameterizationFilter
  *   \brief This class rotates the parameterization of  an 
  *          itk::SphericalParameterizedTriangleMesh.

  * A rotation of the parameterization results in rotating all mapped 
  * points/landmarks around the mesh.
  * The function SetTransform() is used to specify the desired rotation 
  * via an itk::AffineTransform.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
            German Cancer Research Center, Heidelberg, Germany.
  */
  template <class TParameterization>
  class RotateSphericalParameterizationFilter : public MeshToMeshFilter<TParameterization, TParameterization>
  {

  public:

    /** Standard class typedefs. */
    typedef RotateSphericalParameterizationFilter   Self;
    typedef MeshToMeshFilter
      <TParameterization, TParameterization>        Superclass;
    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;

    /** Convenient typedefs. */
    typedef TParameterization                               ParameterizationType;
    typedef typename ParameterizationType::Pointer          ParameterizationPointer;
    typedef typename ParameterizationType::SphericalMapType SphericalMapType;
    typedef typename ParameterizationType::PointType        PointType;
    typedef typename ParameterizationType::VectorType       VectorType;
    typedef typename PointType::CoordRepType                CoordRepType;
    typedef AffineTransform<CoordRepType, 3>                TransformType;
    typedef typename TransformType::Pointer                 TransformPointer;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(RotateSphericalParameterizationFilter, MeshToMeshFilter);

    /** Sets the input parameterization. */
    void SetInput( ParameterizationType *input )
    {
      m_NewInputSet = true;
      Superclass::SetInput( input );
    }

    /** Get/Set the used transformation.
    * The set transformation has to be a rotation, i.e. no shearing, 
    * translation or scaling, otherwise the result will not be valid.
    */
    itkSetObjectMacro( Transform, TransformType );
    itkGetConstObjectMacro( Transform, TransformType );


  protected:

    RotateSphericalParameterizationFilter();

    ~RotateSphericalParameterizationFilter();

    void PrintSelf( std::ostream& os, Indent indent ) const;

    virtual void GenerateData();

    /** Copies the cell information from landmarks to output and initializes the required point containers. */
    void InitializeOutput( ParameterizationType *input );


  private:

    bool              m_NewInputSet;
    TransformPointer  m_Transform;

  };
  
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRotateSphericalParameterizationFilter.txx"
#endif

#endif
