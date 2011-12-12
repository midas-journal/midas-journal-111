#ifndef _itkStatisticalShapeModel3DCalculator_h
#define _itkStatisticalShapeModel3DCalculator_h

#include "itkProcessObject.h"
#include "itkParameterizedTriangleMesh.h"
#include "itkRemeshParameterizedMeshFilter.h"
#include "itkProcrustesAlign3DMeshFilter.h"
#include "itkGaussianWarpSphericalParameterizationFilter.h"
#include "itkRotateSphericalParameterizationFilter.h"
#include "itkMultiThreader.h"
#include "vnl/vnl_matrix.h"



namespace itk
{

  template<class TShapeModelCalculator>
  class ShapeModelCalculatorCostFunction;


  /** \class StatisticalShapeModel3DCalculator
  *   \brief Finds correspondences across a set of 2-manifold triangular 
  *          training meshes. 
  *
  * After calling SetNumberOfInputs() to pass the number of meshes in the set, 
  * SetInput() can be used to specify the parameterized input meshes. Since the 
  * parameterizations will be changed in the course of filter execution, the 
  * initial parameterizations should not be that critical, but conformal 
  * parameterizations (as generated with the 
  * itk::ConformalSphericalParameterizationFilter) have shown to deliver good 
  * results. 
  * In addition to the input meshes, the user has to specify a landmark mesh in
  * parameter space that will be used to create the corresponding shapes. For 
  * the used spherical parameterizations, this means the desired number of 
  * vertices for the model should be spread equally on the unit sphere. The 
  * necessary coordinates can be obtained e.g. by subdividing one of the
  * platonic solids and project all points to the sphere.
  * The search for optimal point correspondences is guided by a cost function 
  * that estimates the quality of a shape model. SetCostFunction() has to be 
  * called to specify the function to be used. The best results have been 
  * obtained with cost functions based on the Minimum Description Length of the
  * shape model.
  * A call of Update() triggers the optimization process, which can take 
  * several hours or even days to finish, depending on the used hardware and 
  * number and complexity of sample meshes. The algorithm should scale linear 
  * with number of samples, number of vertices on each sample and number of 
  * landmarks. A rule of thumb is to decimate the sample meshes to reach a 
  * similar number of vertices as the used number of landmarks.
  * The optimization uses a gradient descent algorithm with fixed step size 
  * that can be set using SetParameterizationWarpStepLength() and 
  * SetParameterizationStartStepLength() for optimization of local 
  * correspondence and rotation, respectively. Convergence is determined by 
  * comparing the value of the cost function every 50 iterations, if the 
  * difference is less than the value set with SetConvergence(), the algorithm
  * raises its level of detail and stops if the maximum level is reached.
  * After that, the resulting point correspondences can be queried using three 
  * methods: GetOutputMesh() returns the corresponding points (landmarks) in 
  * the original coordinate system of each mesh, GetOutputAlignedMesh() returns
  * the landmarks in the aligned coordinate system used to build the model and 
  * GetOutputParameterization() returns the optimal parameterization for each 
  * input.
  * By deriving from itk::StatisticalShapeModel3DGenerator and overwriting 
  * GetNumberOfComponents() and InitializeMatrix(), it is possible to use 
  * arbitrary features for the optimization of correspondence.
  *
  * * \par REFERENCES
  * \par 
  * [1] Heimann T., Wolf I., Williams T., Meinzer H.-P.:
  *     3D Active Shape Models Using Gradient Descent Optimisation of 
  *     Description Length. In: Proc. IPMI 2005, pp. 566-577.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
  *         German Cancer Research Center, Heidelberg, Germany.
  */
  template<class TLandmarks, class TParameterization, class TOutputMeshes>
  class StatisticalShapeModel3DCalculator : public Object
  {

    friend class ShapeModelCalculatorCostFunction<StatisticalShapeModel3DCalculator>;


  public:

    /** Standard typedefs. */
    typedef StatisticalShapeModel3DCalculator      Self;
    typedef Object                            Superclass;
    typedef SmartPointer<Self>                Pointer;
    typedef SmartPointer<const Self>          ConstPointer;

    /** Convenient typedefs. */
    typedef TLandmarks                                  LandmarkType;
    typedef TParameterization                           ParameterizationType;
    typedef TOutputMeshes                               OutputMeshType;
    typedef typename LandmarkType::Pointer              LandmarkPointer;
    typedef typename OutputMeshType::Pointer            OutputMeshPointer;
    typedef typename ParameterizationType::VectorType   VectorType;
    typedef typename ParameterizationType::Pointer      ParameterizationPointer;
    typedef std::vector<ParameterizationPointer>        ParameterizationPointerArray;
    typedef RotateSphericalParameterizationFilter
      <ParameterizationType>                            RotateParameterizationFilterType;
    typedef typename RotateParameterizationFilterType::
                                               Pointer  RotateParameterizationFilterPointer;
    typedef std::vector
      <RotateParameterizationFilterPointer>             RotateParameterizationFilterPointerArray;
    typedef GaussianWarpSphericalParameterizationFilter
      <ParameterizationType>                            ReparameterizationFilterType;
    typedef typename ReparameterizationFilterType::
                                               Pointer  ReparameterizationFilterPointer;
    typedef std::vector<ReparameterizationFilterPointer>    ReparameterizationFilterPointerArray;
    typedef RemeshParameterizedMeshFilter
      <ParameterizationType, LandmarkType, OutputMeshType>  RemeshingType;
    typedef typename RemeshingType::Pointer                 RemeshingPointer;
    typedef std::vector<RemeshingPointer>               RemeshingPointerArray;
    typedef ProcrustesAlign3DMeshFilter
      <OutputMeshType, OutputMeshType>                  ProcrustesAlignFilterType;
    typedef typename OutputMeshType::CoordRepType       CoordRepType;
    typedef typename ProcrustesAlignFilterType::TransformType    TransformType;
    typedef typename TransformType::Pointer             TransformPointer;
    typedef TransformMeshFilter
      <OutputMeshType, OutputMeshType, TransformType>   TransformMeshType;
    typedef typename TransformMeshType::Pointer         TransformMeshPointer;
    typedef std::vector<TransformMeshPointer>           TransformMeshArray;
    typedef vnl_matrix<double>                          MatrixType;
    typedef std::vector<MatrixType>                     MatrixArray;
    typedef ShapeModelCalculatorCostFunction<Self>      CostFunctionType;
    
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Standard part of every itk Object. */
    itkTypeMacro(StatisticalShapeModel3DCalculator, Object);

    /** Sets the number of meshes to be analyzed. */
    void SetNumberOfInputs( unsigned int number );

    /** Returns the number of analyzed meshes. */
    itkGetConstMacro( NumberOfInputs, unsigned int );

    /** Sets the specified parameterized input mesh.
    * Copies the mesh information  when the method is called, i.e. the mesh
    * and its parameterization have to be initialized before.
    */
    void SetInput( unsigned int idx, ParameterizationPointer parameterization );

    /** Returns the optimal landmarks for the specified input,
    * using the original coordinate system.
    * Use only after calling Update().
    */
    OutputMeshPointer GetOutputMesh( unsigned int idx )
    {
      return m_ShapeGenerators[idx]->GetOutput();
    }

    /** Returns the optimal landmarks for the specified input,
    * using the internal, normalized coordinate system.
    * Use only after calling Update().
    */
    OutputMeshPointer GetOutputAlignedMesh( unsigned int idx )
    {
      return m_Transforms[idx]->GetOutput();
    }

    /** Returns the optimal parameterization for the specified input.
    * Use only after calling Update().
    */
    ParameterizationPointer GetOutputParameterization( unsigned int idx )
    {
      return m_Parameterizations[idx];
    }

    OutputMeshPointer GetResampledOutputMesh( unsigned int idx ) ;
    
    /** Determines if the algorithm should also optimize the start offsets of 
    * the parameterizations, which practically leads to optimizing the rotation
    * of the sample meshes. If all samples have the same orientation, this can 
    * be turned off for a slightly faster optimization. The default is on.
    */
    itkSetMacro( OptimizeParameterizationStart, bool );
    itkGetMacro( OptimizeParameterizationStart, bool );
    itkBooleanMacro( OptimizeParameterizationStart );

    /** Set/Get the convergence value. */
    itkSetMacro( Convergence, double );
    itkGetConstMacro( Convergence, double );

    /** Set/Get the step length for optimization of local parameterization. */
    itkSetMacro( ParameterizationWarpStepLength, double );
    itkGetConstMacro( ParameterizationWarpStepLength, double );
    
    /** Set/Get the step length for optimization of parameterization start. */
    itkSetMacro( ParameterizationStartStepLength, double );
    itkGetConstMacro( ParameterizationStartStepLength, double );
    
    /** Sets the used landmarks in parameter space. */
    virtual void SetLandmarks( LandmarkPointer landmarks );
    
    virtual void SetOutputLandmarks ( LandmarkPointer landmarks ) ;

    /** Returns the number of used landmarks. */
    unsigned int GetNumberOfLandmarks() const
    {
      return m_Landmarks->GetNumberOfPoints();
    }

    /** Sets the cost function used for optimization. */
    void SetCostFunction( CostFunctionType *costFunction )
    {
      m_CostFunction = costFunction;
      m_CostFunction->SetModel( this );
      this->Modified();
    }
    
    /** Triggers the optimization process. */
    void Update()
    {
      this->GenerateData();
    }
   
    /** Returns the state of the automatic alignment. */
    bool GetUseAutomaticAlignment() const
    {
      return m_ProcrustesFilter->GetAlignRotation();
    }

    /** Determines if the algorithm aligns the different input meshes 
    * automatically for best fit. This means a generalized procrustes 
    * alignment plus an individual scaling to tangent size.
    * Default is on. If switched off, the shapes are only translated
    * to zero origin (which is necessary for the optimization).
    */
    void SetUseAutomaticAlignment( bool on )
    {
      m_ProcrustesFilter->SetAlignRotation( on );
      m_ProcrustesFilter->SetAlignScale( on );
      m_UseTangentScale = on;
    }

    itkBooleanMacro( UseAutomaticAlignment );
    

  protected:

    StatisticalShapeModel3DCalculator();
    ~StatisticalShapeModel3DCalculator();

    void PrintSelf( std::ostream& os, Indent indent ) const;

    virtual void GenerateData();

    void OptimizeParameterizationStarts();

    void OptimizeParameterizations();

    void BuildModel();

    void RotateParameterizationsAndLandmarks();

    double GetTangentSize( unsigned int idx ) const;

    /** Uses finite differences to compute the gradients of column sampleIdx 
    * in the AlignedData matrix with respect to movements of the control points
    * on the parameterization of sampleIdx.
    * Each column of gradients corresponds to moving a certain control point in
    * a certain direction. 
    */
    void ComputeReparameterizationGradients( unsigned int sampleIdx, MatrixType *gradients ) const;

    /** Starts multiple threads to calculate the gradients for reparameterizing each sample. */
    void ComputeAllReparameterizationGradients( MatrixArray &contourGradient ) const;

    /** Calculates the gradients for rotating the parameterizations. */
    void ComputeAllParameterizationStartGradients( MatrixArray &contourGradient ) const;

    /** Initializes column sampleIdx of matrix, given the values in mesh. 
    * The x-coordinates are stored in matrix[3*pointId+0][sampleIdx],
    * y-coordinates in matrix[3*pointId+1][sampleIdx] and z-coordinates in
    * matrix[3*pointId+2][sampleIdx].
    */
    virtual void InitializeMatrix( MatrixType *matrix, unsigned int sampleIdx, OutputMeshPointer mesh ) const;

    /** Returns the number of components each point of the mesh has in the 
    * data matrix used for the PCA.
    * This are x,y and z coordinates, i.e. 3 components.
    */
    virtual unsigned int GetNumberOfComponents() const
    {
      return 3;
    }

    void DebugHelper() const;

    //  private:

    unsigned int                                m_NumberOfInputs;
    LandmarkPointer                             m_Landmarks;
    LandmarkPointer                             m_OutputLandmarks ;
    TransformMeshPointer                        m_LandmarkRotator;
    ParameterizationPointerArray                m_Parameterizations;
    RotateParameterizationFilterPointerArray    m_ParameterizationRotator;
    /** Filters to change the current parameterizations */
    ReparameterizationFilterPointerArray        m_ReparameterizationFilters;
    RemeshingPointerArray                       m_ShapeGenerators;
    typename ProcrustesAlignFilterType::Pointer m_ProcrustesFilter;
    TransformMeshArray                          m_Transforms;
    MatrixType                                  m_AlignedData;
    RemeshingPointerArray                       m_GradientShapeGenerators;
    TransformMeshArray                          m_GradientTransforms;
    bool                                        m_OptimizeParameterizationStart;
    bool                                        m_UseTangentScale;
    double                                      m_ParameterizationWarpStepLength;
    double                                      m_ParameterizationStartStepLength;
    double                                      m_Convergence;
    unsigned int                                m_CurrentNumberOfIterations;
    std::vector<double>                         m_SingularValues;
    std::vector<double>                         m_EigenValues;
    MatrixType                                  m_SVDMatrixU;
    MatrixType                                  m_SVDMatrixV;
    typename MultiThreader::Pointer             m_MultiThreader;
    unsigned int                                m_NumberOfThreads;

    CostFunctionType                            *m_CostFunction;
    
    OutputMeshPointer                           m_ResampleResult ;    

    /** Thread-Data Structure   */
    struct ThreadDataStruct
    {
      const Self      *modelCalculator;
      unsigned int    index;
      MatrixType      *gradients;
    };

    static ITK_THREAD_RETURN_TYPE ComputeReparameterizationGradientsThreaderCallback( void *arg );

  };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStatisticalShapeModel3DCalculator.txx"
#endif

#endif
