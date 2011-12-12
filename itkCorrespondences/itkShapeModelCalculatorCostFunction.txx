#ifndef _itkShapeModelCalculatorCostFunction_txx
#define _itkShapeModelCalculatorCostFunction_txx

#include "itkShapeModelCalculatorCostFunction.h"


namespace itk
{

  template <class TShapeModelCalculator>
  ShapeModelCalculatorCostFunction<TShapeModelCalculator>
  ::ShapeModelCalculatorCostFunction()
  {
    m_Model = 0;
  }


  template <class TShapeModelCalculator>
  ShapeModelCalculatorCostFunction<TShapeModelCalculator>
  ::~ShapeModelCalculatorCostFunction()
  {
  }


  template <class TShapeModelCalculator>
  typename ShapeModelCalculatorCostFunction<TShapeModelCalculator>::VectorType
  ShapeModelCalculatorCostFunction<TShapeModelCalculator>
  ::GetGradientVector( unsigned int sampleIdx, unsigned int controlPntIdx ) const
  {
    VectorType moveDir;
    moveDir[0] = GetGradient( sampleIdx, 2*controlPntIdx );
    moveDir[1] = GetGradient( sampleIdx, 2*controlPntIdx+1 );
    moveDir[2] = 0.0;
    return moveDir;
  }

}

#endif
