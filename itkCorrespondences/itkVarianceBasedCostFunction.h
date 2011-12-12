#ifndef _itkVarianceBasedCostFunction_h
#define _itkVarianceBasedCostFunction_h

#include "itkShapeModelCalculatorCostFunction.h"
#include "itkMultiThreader.h"


namespace itk
{

  /** \class VarianceBasedCostFunction
  * \brief Implements the PrepareGradients() method for all cost functions
  *        based on the variance, i.e. the eigenvalues of the model. 
  *
  * Virtual functions that have to be implemented in subclasses are GetValue()
  * and GetGradient(). Both can make use of the eigenvalue gradients stored in 
  * the EVGradients matrix.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
            German Cancer Research Center, Heidelberg, Germany.
  */
  template<class TShapeModelCalculator>
  class VarianceBasedCostFunction : public ShapeModelCalculatorCostFunction<TShapeModelCalculator>
  {

  public:

    /** Standard class typedefs. */
    typedef VarianceBasedCostFunction           Self;
    typedef ShapeModelCalculatorCostFunction
      <TShapeModelCalculator>                   Superclass;
    typedef SmartPointer<Self>                  Pointer;
    typedef SmartPointer<const Self>            ConstPointer;

    /** Convenient typedefs. */
    typedef typename Superclass::ShapeModelCalculatorType ShapeModelCalculatorType;
    typedef typename Superclass::MeasureType              MeasureType;
    typedef typename Superclass::VectorType               VectorType;
    typedef typename Superclass::MatrixType               MatrixType;
    typedef typename Superclass::MatrixArray              MatrixArray;

    /** Run-time type information (and related methods). */
    itkTypeMacro(VarianceBasedCostFunction, ShapeModelCalculatorCostFunction);

    /** Prepares the cost gradients by specifying the gradients of all coordinates.
    * Stores the results in the m_EVGradients matrix.
    */
    virtual void PrepareGradients( MatrixArray &dataGradients );


  protected:

    VarianceBasedCostFunction();
    ~VarianceBasedCostFunction();
    
    void ComputeEigenValueGradients( unsigned int ev, MatrixType *evGradients, MatrixArray *dataGradients ) const;

    static ITK_THREAD_RETURN_TYPE ComputeEigenValueGradientsThreaderCallback( void *arg );
  
    double GetMatrixColumnDotProduct( const MatrixType *m1, unsigned int col1, const MatrixType *m2, unsigned int col2 ) const;


    MatrixArray   m_EVGradients;


    typename MultiThreader::Pointer     m_MultiThreader;
    unsigned int                        m_NumberOfThreads;

    /** Thread-Data Structure   */
    struct ThreadDataStruct
    {
      const Self      *costFunction;
      unsigned int    index;
      MatrixType      *gradients;
      MatrixArray     *cgArray;
    };


  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVarianceBasedCostFunction.txx"
#endif

#endif
