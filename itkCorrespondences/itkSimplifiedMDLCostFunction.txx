#ifndef _itkSimplifiedMDLCostFunction_txx
#define _itkSimplifiedMDLCostFunction_txx

#include "itkSimplifiedMDLCostFunction.h"


namespace itk
{

  template<class TShapeModelCalculator>
  SimplifiedMDLCostFunction<TShapeModelCalculator>
  ::SimplifiedMDLCostFunction()
  {
    m_VarianceCut = 1.0;
  }


  template<class TShapeModelCalculator>
  SimplifiedMDLCostFunction<TShapeModelCalculator>
  ::~SimplifiedMDLCostFunction()
  {
  }


  template <class TShapeModelCalculator>
  void
  SimplifiedMDLCostFunction<TShapeModelCalculator>
  ::PrintSelf( std::ostream& os, Indent indent ) const
  {
    Superclass::PrintSelf( os, indent );
  }


  template <class TShapeModelCalculator>
  void
  SimplifiedMDLCostFunction<TShapeModelCalculator>
  ::SetVarianceCutForModelRadius( double radius )
  {
    double stdDevCut = 0.3 / radius;
    this->SetVarianceCut( stdDevCut*stdDevCut );
  }


  template <class TShapeModelCalculator>
  typename SimplifiedMDLCostFunction<TShapeModelCalculator>::MeasureType
  SimplifiedMDLCostFunction<TShapeModelCalculator>
  ::GetValue() const
  {
    double result = 0;
    for (unsigned int i=0; i<this->m_Model->GetNumberOfInputs()-1; i++) {
      double ev = this->GetEigenValue( i );
      if (ev >= m_VarianceCut) {
        result += 1.0 + log( ev / m_VarianceCut );
      }
      else {
        result += ev / m_VarianceCut;
      }		
    }
    return result;
  }


  template <class TShapeModelCalculator>
  typename SimplifiedMDLCostFunction<TShapeModelCalculator>::MeasureType
  SimplifiedMDLCostFunction<TShapeModelCalculator>
  ::GetGradient( unsigned int sampleIdx, unsigned int componentIdx) const
  {
    unsigned int numSamples = this->m_Model->GetNumberOfInputs();
    MeasureType gradient[2] = {0, 0};
    MeasureType divValue;
    
    for (unsigned int ev=0; ev<numSamples-1; ev++) 
    {
      divValue = this->GetEigenValue( ev ) ;
      if ( m_VarianceCut > divValue )
      {
        divValue = m_VarianceCut ;
      }
      for (int g=0; g<2; g++) 
      {
          gradient[g] += (this->m_EVGradients[ev])[sampleIdx][2*componentIdx+g] / divValue;
      }
    }
    
    if (gradient[0] < 0) { return -gradient[0]; }
    else 
    {
      if (gradient[1] < 0) { return gradient[1]; }
      else                 { return 0.0;         }
    }
  }

}

#endif
