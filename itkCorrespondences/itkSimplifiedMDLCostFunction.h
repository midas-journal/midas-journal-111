#ifndef _itkSimplifiedMDLCostFunction_h
#define _itkSimplifiedMDLCostFunction_h

#include "itkVarianceBasedCostFunction.h"


namespace itk
{

  /** \class BoundedCovarianceMagnitudeCostFunction
  * \brief Implements a costfunction using the sum of the log values of the
  * bounded Eigenvalues of the Covariance matrix.
  *
  * The gradients are calculated based on the objective function of Thodberg [2],
  * while some constant terms are added to the function values to bring 
  * it closer to a value similar to Minimum Description Length by Davies[1].
  * Before using this cost function, the 
  * variance threshold that determines the amount of noise in the image data 
  * has to be set. The easiest way to set this threshold is to use 
  * SetVarianceCutForModelRadius() and pass the average voxel radius of the 
  * sample shapes as argument.
  *
  * \par REFERENCES
  * \par 
  * [1] Davies, R.H., Twining, C.J., Cootes, T.F., Waterton, J.C., Taylor C.J.:
  *     3D Statistical Shape Models Using Direct Optimisation of Description
  *     Length. In: Proc. ECCV 2002, Part 3, pp. 3-20.
  * \par
  * [2] Thodberg H.H.: Minimum Description Length Shape And Appearance Models.
  *     In: Proc. IPMI 2003, pp. 51-62.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
            German Cancer Research Center, Heidelberg, Germany.
  */
  template<class TShapeModelCalculator>
  class SimplifiedMDLCostFunction : public VarianceBasedCostFunction<TShapeModelCalculator>
  {

  public:

    /** Standard class typedefs. */
    typedef SimplifiedMDLCostFunction           Self;
    typedef VarianceBasedCostFunction
      <TShapeModelCalculator>                   Superclass;
    typedef SmartPointer<Self>                  Pointer;
    typedef SmartPointer<const Self>            ConstPointer;

    /** Convenient typedefs. */
    typedef typename Superclass::ShapeModelCalculatorType ShapeModelCalculatorType;
    typedef typename Superclass::MeasureType              MeasureType;
    typedef typename Superclass::VectorType               VectorType;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SimplifiedMDLCostFunction, VarianceBasedCostFunction);

    /** Returns the costs for the current model. */
    virtual MeasureType GetValue() const;

    /** Returns the gradient for the current model. */
    virtual MeasureType GetGradient( unsigned int sampleIdx, unsigned int componentIdx) const;

    /** Get/Set the variance threshold for the cost function. */
    itkSetMacro( VarianceCut, double );
    itkGetConstMacro( VarianceCut, double );

    /** Sets the variance threshold for the cost function depending on model 
    * radius (in voxels). 
    */
    void SetVarianceCutForModelRadius( double radius );
    

  protected:

    SimplifiedMDLCostFunction();

    ~SimplifiedMDLCostFunction();

    void PrintSelf( std::ostream& os, Indent indent ) const;


  private:

    double m_VarianceCut;

  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSimplifiedMDLCostFunction.txx"
#endif

#endif
