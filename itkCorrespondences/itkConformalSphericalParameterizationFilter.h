#ifndef _itkConformalSphericalParameterizationFilter_h
#define _itkConformalSphericalParameterizationFilter_h

#include "itkMeshToMeshFilter.h"
#include "itkSphericalParameterizedTriangleMesh.h"

namespace itk
{

  /** \class ConformalSphericalParameterizationFilter
  *   \brief The class generates a conformal parameterization for meshes of 
  *          spherical topology. 
  *
  * It expects an itk::IndexedTriangleMesh as input and delivers an 
  * itk::SphericalParameterizedTriangleMesh as output. 
  * Internally, the first step of the algorithm is to compute a Gauss map of 
  * the mesh, where all points are mapped to positions specified by their 
  * normal vectors. The protected method ComputeGaussMap() calculates all face 
  * normals coherently (i.e. pointing in the correct direction) and averages 
  * the face normals around each point.
  * Subsequently, the string energy of the mesh is minimized in two steps, 
  * first optimizing the barycentric energy, then the conformal energy. 
  * In some rare cases, these optimizations do not converge and the filter 
  * will not return from the Update() function. 
  *
  * \par REFERENCES
  * \par 
  * [1] Gu X., Wang Y., et. al, Genus zero surface conformal mapping and its application to brain surface mapping,
  *     Proc. IPMI 2003, pp. 172-184, 2003.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, German Cancer Research Center, Heidelberg, Germany.
  */
  template <class TIndexedTriangleMesh, class TParameterizedTriangleMesh>
  class ConformalSphericalParameterizationFilter : public MeshToMeshFilter<TIndexedTriangleMesh, TParameterizedTriangleMesh>
  {

  public:

    /** Standard class typedefs. */
    typedef ConformalSphericalParameterizationFilter  Self;
    typedef MeshToMeshFilter<TIndexedTriangleMesh, 
              TParameterizedTriangleMesh>             Superclass;
    typedef SmartPointer<Self>                        Pointer;
    typedef SmartPointer<const Self>                  ConstPointer;

    /** Convenient typedefs. */
    typedef DataObject::Pointer                             DataObjectPointer;
    typedef TIndexedTriangleMesh                            InputMeshType;
    typedef TParameterizedTriangleMesh                      ParameterizationType;
    typedef typename ParameterizationType::Pointer          ParameterizationPointer;
    typedef typename ParameterizationType::SphericalMapType SphericalMapType;
    typedef typename InputMeshType::PointType               PointType;
    typedef typename InputMeshType::VectorType              VectorType;
    typedef typename InputMeshType::Pointer                 IndexedTriangleMeshPointer;


    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ConformalSphericalParameterizationFilter, MeshToMeshFilter);

    /** Get/Set the step length for the minimization of the Tutte energy. */
    itkSetMacro( TutteStepLength, double );
    itkGetMacro( TutteStepLength, double );
    /** Get/Set the convergence value for the minimization of the Tutte energy. */
    itkSetMacro( TutteConvergence, double );
    itkGetMacro( TutteConvergence, double );
    /** Get/Set the step length for the minimization of the harmonic energy. */
    itkSetMacro( HarmonicStepLength, double );
    itkGetMacro( HarmonicStepLength, double );
    /** Get/Set the convergence value for the minimization of the harmonic energy. */
    itkSetMacro( HarmonicConvergence, double );
    itkGetMacro( HarmonicConvergence, double );


  protected:

    ConformalSphericalParameterizationFilter();

    ~ConformalSphericalParameterizationFilter();

    /** Generate requested data. */
    virtual void GenerateData( void );

    /** Computes a Gauss map of the input mesh as initial start. */
    void ComputeGaussMap();

    /** Minimizes the string energy of the mesh, using the values in 
    * m_StringConstant. The iterative process moves all points by m_StepLength 
    * until the energy decrease is less than m_Convergence.
    */
    void MinimizeEnergy();

    /** Translates the spherical map to an average position of zero. */
    void CenterSphere();

    /** Projects all points to the unit sphere. */
    void RemapToSphere();

    /** Returns the current string energy of the mesh, using the values in 
    * m_StringConstant.
    */
    double GetEnergy();

    /** Calculates the piecewice Laplacian at the given point. */
    VectorType GetPiecewiseLaplacian( unsigned int pointId );

    /** Calculates the absolute derivative at the given point. */
    VectorType GetAbsoluteDerivative( unsigned int pointId );

    /** Calculates the cotangens of an angle opposing the given edge.
    * localFaceId is 0 or 1 and specifies which of the angles is calculated.
    */
    double GetCotangensAgainstEdge( unsigned int edgeId, int localFaceId );


  private:

    SphericalMapType                *m_Sphere;
    IndexedTriangleMeshPointer      m_SourceMesh;
    std::vector<double>             m_StringConstant;
    double                          m_TutteStepLength;
    double                          m_TutteConvergence;
    double                          m_HarmonicStepLength;
    double                          m_HarmonicConvergence;
    double                          m_StepLength;
    double                          m_Convergence;
      

  };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConformalSphericalParameterizationFilter.txx"
#endif

#endif
