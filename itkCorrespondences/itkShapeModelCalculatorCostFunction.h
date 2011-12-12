#ifndef _itkShapeModelCalculatorCostFunction_h
#define _itkShapeModelCalculatorCostFunction_h

#include "itkObject.h"
#include "itkStatisticalShapeModel3DCalculator.h"


namespace itk
{

  
  /** \class ShapeModelCalculatorCostFunction
  *   \brief Baseclass for cost functions to be used with the 
  *          StatisticalShapeModel3DCalculator.
  * 
  * The StatisticalShapeModel3DCalculator class depends on external cost 
  * functions to guide the model building process.
  * Virtual functions that have to be implemented in subclasses are GetValue(), 
  * GetGradient() and PrepareGradients(). While GetValue() returns the current 
  * costs for the model, GetGradient() returns the gradient in a certain 
  * direction. How this direction is interpreted is specified by the 
  * implementation of PrepareGradients().
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
            German Cancer Research Center, Heidelberg, Germany.
  */
  template<class TShapeModelCalculator>
  class ShapeModelCalculatorCostFunction : public Object
  {

  public:

    /** Standard class typedefs. */
    typedef ShapeModelCalculatorCostFunction  Self;
    typedef Object                            Superclass;
    typedef SmartPointer<Self>                Pointer;
    typedef SmartPointer<const Self>          ConstPointer;

    /** Convenient typedefs. */
    typedef TShapeModelCalculator             ShapeModelCalculatorType;
    typedef double                            MeasureType;
    typedef typename ShapeModelCalculatorType
      ::VectorType                            VectorType;
    typedef typename ShapeModelCalculatorType
      ::ReparameterizationFilterType          ReparameterizationFilterType;
    typedef typename ShapeModelCalculatorType
      ::RemeshingType                         RemeshingFilterType;
    typedef typename ShapeModelCalculatorType
      ::OutputMeshType                        OutputMeshType;
    typedef vnl_matrix<double>                MatrixType;
    typedef std::vector<MatrixType>           MatrixArray;
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(ShapeModelCalculatorCostFunction, Object);

    /** Returns the costs for the current model. */
    virtual MeasureType GetValue() const = 0;

    /** Returns the gradient for the current model. */
    virtual MeasureType GetGradient( unsigned int sampleIdx, unsigned int componentIdx) const = 0;

    /** Returns the gradient for a given control point of the current model. */
    virtual VectorType GetGradientVector( unsigned int sampleIdx, unsigned int controlPntIdx ) const;
        
    /** Prepares the cost gradients by specifying the gradients of all coordinates. */
    virtual void PrepareGradients( MatrixArray &dataGradients ) = 0;

    /** Sets the shape model that is used for all calculations. */
    itkSetConstObjectMacro( Model, ShapeModelCalculatorType );


  protected:

    ShapeModelCalculatorCostFunction();

    ~ShapeModelCalculatorCostFunction();

    double GetEigenValue( unsigned int idx ) const
    {
      return m_Model->m_EigenValues[idx];
    }

    double GetSingularValue( unsigned int idx ) const
    {
      return m_Model->m_SingularValues[idx];
    }

    const MatrixType* GetSVDMatrixU() const
    {
      return &(m_Model->m_SVDMatrixU);
    }

    const MatrixType* GetSVDMatrixV() const
    {
      return &(m_Model->m_SVDMatrixV);
    }
  
    void ComputeAllReparameterizationGradients( MatrixArray &contourGradient ) const
    {
      m_Model->ComputeAllReparameterizationGradients( contourGradient );
    }

    void ComputeAllParameterizationStartGradients( MatrixArray &contourGradient ) const
    {
      m_Model->ComputeAllParameterizationStartGradients( contourGradient );
    }

    /** The shape model calculator delivering the values for the cost function. */
    const ShapeModelCalculatorType      *m_Model;

  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShapeModelCalculatorCostFunction.txx"
#endif

#endif
