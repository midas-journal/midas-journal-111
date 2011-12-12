#ifndef _itkStatisticalShapeModel3DCalculatorWithFeatures_h
#define _itkStatisticalShapeModel3DCalculatorWithFeatures_h

#include "itkStatisticalShapeModel3DCalculator.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkMeanCalculator.h"


#define NUM_COMPONENTS 2

namespace itk
{
  template<class TShapeModelCalculator>
  class ShapeModelCalculatorCostFunction;


  /** \class StatisticalShapeModel3DCalculatorWithFeatures
  *   \brief Finds correspondences across a set of 2-manifold triangular 
  *          training meshes. 
  *
  * This class derives from StatisticalShapeModel3DCalculator to make it capable of using
  * arbitrary features for correspondence optimization. 
  * The major changes are the overwritten GetNumberOfComponents() and 
  * InitializeMatrix() methods, which use the point data associated with the mesh, 
  * rather than the vertex coordinates. 
  * The additional ComputeFeatureMeans() method computes the feature means once and 
  * stores it in the new class member m_FeatureMeans, so that it can be reused every 
  * time it is needed in the InitializeMatrix() method.
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
  class StatisticalShapeModel3DCalculatorWithFeatures : 
    public StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  {

    friend class ShapeModelCalculatorCostFunction<StatisticalShapeModel3DCalculatorWithFeatures>;


  public:

    /** Standard typedefs. */
    typedef StatisticalShapeModel3DCalculatorWithFeatures   Self;
    typedef StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes> Superclass;
    typedef SmartPointer<Self>                              Pointer;
    typedef SmartPointer<const Self>                        ConstPointer;

    typedef ShapeModelCalculatorCostFunction<Self>          CostFunctionType;
    typedef typename Superclass::LandmarkType               LandmarkType;
    typedef typename Superclass::LandmarkPointer            LandmarkPointer;
    typedef typename Superclass::MatrixType                 MatrixType;
    typedef typename Superclass::OutputMeshPointer          OutputMeshPointer;
    typedef typename LandmarkType::PixelType                FeatureVectorType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Standard part of every itk Object. */
    itkTypeMacro(StatisticalShapeModel3DCalculatorWithFeatures, Object);

    /** Sets the cost function used for optimization. */
    void SetCostFunction( CostFunctionType * costFunction )
    {
      Superclass::SetCostFunction (( typename Superclass::CostFunctionType *) costFunction) ;
    }

    virtual void SetLandmarks( LandmarkPointer landmarks ) ;

    
  protected:

    StatisticalShapeModel3DCalculatorWithFeatures();
    ~StatisticalShapeModel3DCalculatorWithFeatures();

    void PrintSelf( std::ostream& os, Indent indent ) const;

    /** Initializes column sampleIdx of matrix, given the values in mesh. 
    * The x-coordinates are stored in matrix[3*pointId+0][sampleIdx],
    * y-coordinates in matrix[3*pointId+1][sampleIdx] and z-coordinates in
    * matrix[3*pointId+2][sampleIdx].
    */
    virtual void InitializeMatrix( MatrixType *matrix, unsigned int sampleIdx, OutputMeshPointer mesh ) const ;

    /** Returns the number of components each point of the mesh has in the 
    * data matrix used for the PCA.
    * These are Koenderick's C and S metrics, i.e. 2 components.
    */
    virtual unsigned int GetNumberOfComponents() const
    {
      return NUM_COMPONENTS ;
    }
    
    
    mutable FeatureVectorType m_FeatureMeans ;

    virtual void ComputeFeatureMeans () const ;  
  };  

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStatisticalShapeModel3DCalculatorWithFeatures.txx"
#endif

#endif
