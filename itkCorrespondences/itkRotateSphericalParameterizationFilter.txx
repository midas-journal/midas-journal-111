#ifndef _itkRotateSphericalParameterizationFilter_txx
#define _itkRotateSphericalParameterizationFilter_txx

#include "itkRotateSphericalParameterizationFilter.h"


namespace itk
{

  template <class TParameterization>
  RotateSphericalParameterizationFilter<TParameterization>
  ::RotateSphericalParameterizationFilter()
  {
    this->ReleaseDataBeforeUpdateFlagOff();
    m_NewInputSet = false;
    m_Transform = TransformType::New();
  }


  template <class TParameterization>
  RotateSphericalParameterizationFilter<TParameterization>
  ::~RotateSphericalParameterizationFilter()
  {
  }


  template <class TParameterization>
  void
  RotateSphericalParameterizationFilter<TParameterization>
  ::PrintSelf( std::ostream& os, Indent indent ) const
  {
    Superclass::PrintSelf( os, indent );
  }

  
  template <class TParameterization>
  void
  RotateSphericalParameterizationFilter<TParameterization>
  ::InitializeOutput( ParameterizationType *input )
  {
    ParameterizationPointer output = this->GetOutput();
    (*output) = input;
    // initialize parameterization
    output->InitializeSphericalMap();
  }


  template <class TParameterization>
  void
  RotateSphericalParameterizationFilter<TParameterization>
  ::GenerateData()
  {
    if (this->GetNumberOfInputs() < 1) 
    { 
      itkExceptionMacro( "Input mesh not set." ); 
    }

    ParameterizationType *input = const_cast<ParameterizationType*>(this->GetInput()); 

    if (m_NewInputSet)
    {
      this->InitializeOutput( input );
      m_NewInputSet = false;
    }
    
    SphericalMapType *inputMap = input->GetSphericalMap();
    SphericalMapType *outputMap = this->GetOutput()->GetSphericalMap();
    for (unsigned int i=0; i<outputMap->size(); i++) 
    {
      (*outputMap)[i] = m_Transform->TransformPoint( (*inputMap)[i] );
    }
    this->GetOutput()->SetBufferedRegion( 0 );
  }

}

#endif
