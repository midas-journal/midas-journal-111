#ifndef _itkGaussianWarpSphericalParameterizationFilter_txx
#define _itkGaussianWarpSphericalParameterizationFilter_txx

#include "itkGaussianWarpSphericalParameterizationFilter.h"


namespace itk
{

  template <class TParameterization>
  GaussianWarpSphericalParameterizationFilter<TParameterization>
  ::GaussianWarpSphericalParameterizationFilter()
  {
    this->ReleaseDataBeforeUpdateFlagOff();
    m_NewInputSet = false;
    m_ActiveControlPoint = 0;
    m_Direction[0] = 0;
    m_Direction[1] = 0;
    this->SetLevelOfDetail( 0 );
  }


  template <class TParameterization>
  GaussianWarpSphericalParameterizationFilter<TParameterization>
  ::~GaussianWarpSphericalParameterizationFilter()
  {
  }


  template <class TParameterization>
  void
  GaussianWarpSphericalParameterizationFilter<TParameterization>
  ::PrintSelf( std::ostream& os, Indent indent ) const
  {
    Superclass::PrintSelf(os, indent);
  }


  template <class TParameterization>
  void
  GaussianWarpSphericalParameterizationFilter<TParameterization>
  ::SetLevelOfDetail( unsigned int level )
  {
    if (level > this->GetMaximumLevelOfDetail())
    {
      itkExceptionMacro( "Invalid level of detail." );
    }
    m_LevelOfDetail = level;
    PointType cpnt;
    if (level==0 || level==1) 
    {
      if (level==0) { m_KernelStandardDeviation = 0.4; }
      else          { m_KernelStandardDeviation = 0.3; }
      itkDebugMacro( "Using 4 control points." );
      m_NumberOfControlPoints = 4;
      m_ControlPoints.resize( m_NumberOfControlPoints );
      cpnt[0]= 1;  cpnt[1]= 0;  cpnt[2]=0;  m_ControlPoints[0]=cpnt;
      cpnt[0]= 0;  cpnt[1]= 1;  cpnt[2]=0;  m_ControlPoints[1]=cpnt;
      cpnt[0]=-1;  cpnt[1]= 0;  cpnt[2]=0;  m_ControlPoints[2]=cpnt;
      cpnt[0]= 0;  cpnt[1]=-1;  cpnt[2]=0;  m_ControlPoints[3]=cpnt;
    }
    else if (level==2) 
    {
      m_KernelStandardDeviation = 0.2;
      itkDebugMacro( "Using 12 control points." );
      m_NumberOfControlPoints = 12;
      m_ControlPoints.resize( m_NumberOfControlPoints );
      double a=0, b=0;
      a = 0;
      b = 6.2831853 * 25/360;
      for (int i=0; i<6; i++) 
      {
        cpnt[0] = cos( a ) * cos( b );
        cpnt[1] = sin( a ) * cos( b );
        cpnt[2] = sin( b );
        m_ControlPoints[i] = cpnt;
        a += 6.2831853 / 6;
      }
      a = 6.2831853 / 12;
      b = -6.2831853 * 25/360;
      for (int i=6; i<12; i++) 
      {
        cpnt[0] = cos( a ) * cos( b );
        cpnt[1] = sin( a ) * cos( b );
        cpnt[2] = sin( b );
        m_ControlPoints[i] = cpnt;
        a += 6.2831853 / 6;
      }
    }
    else if (level==3) 
    {
      m_KernelStandardDeviation = 0.1;
      itkDebugMacro( "Using 24 control points.\n" );
      m_NumberOfControlPoints = 24;
      m_ControlPoints.resize( m_NumberOfControlPoints );
      double a=0, b=0;
      for (int i=0; i<8; i++) 
      {
        cpnt[0] = cos( a ) * cos( b );
        cpnt[1] = sin( a ) * cos( b );
        cpnt[2] = sin( b );
        m_ControlPoints[i] = cpnt;
        a += 6.2831853 / 8;
      }
      a = 6.2831853 / 16;
      b = 6.2831853 * 35/360;
      for (int i=8; i<16; i++) 
      {
        cpnt[0] = cos( a ) * cos( b );
        cpnt[1] = sin( a ) * cos( b );
        cpnt[2] = sin( b );
        m_ControlPoints[i] = cpnt;
        a += 6.2831853 / 8;
      }
      a = 6.2831853 / 16;
      b = -6.2831853 * 35/360;
      for (int i=16; i<24; i++) 
      {
        cpnt[0] = cos( a ) * cos( b );
        cpnt[1] = sin( a ) * cos( b );
        cpnt[2] = sin( b );
        m_ControlPoints[i] = cpnt;
        a += 6.2831853 / 8;
      }
    }
  }


  template <class TParameterization>
  void
  GaussianWarpSphericalParameterizationFilter<TParameterization>
  ::GenerateData()
  {
    if (this->GetNumberOfInputs() < 1) 
    { 
      itkExceptionMacro( "Input mesh not set." ); 
    }

    ParameterizationType *input = const_cast<ParameterizationType*>(this->GetInput()); 
    ParameterizationType *output = this->GetOutput();

    if (m_NewInputSet)
    {
      this->InitializeOutput( input );
      m_NewInputSet = false;
    }
    double dir[2];
    // *2pi because change of 1 should be one rotation
    dir[0] = m_Direction[0] * 6.2831853;
    dir[1] = m_Direction[1] * 6.2831853;

    double distThresh = 9 * m_KernelStandardDeviation * m_KernelStandardDeviation;    // 3 standard deviations
    double expDiv = 2 * m_KernelStandardDeviation * m_KernelStandardDeviation;
    double zeroVal = exp( - distThresh / expDiv );
    PointType cPos = m_ControlPoints[m_ActiveControlPoint]; 
    PointType vPos;
    double distance, newAzimuth, newPolar, sinPolar;

    SphericalMapType *inputMap = input->GetSphericalMap();
    SphericalMapType *outputMap = output->GetSphericalMap();
    for (unsigned int i=0; i<outputMap->size(); i++) 
    {
      vPos = (*inputMap)[i];
      double sqrDist = (cPos[0]-vPos[0])*(cPos[0]-vPos[0]) + 
                       (cPos[1]-vPos[1])*(cPos[1]-vPos[1]) + 
                       (cPos[2]-vPos[2])*(cPos[2]-vPos[2]);
      if (sqrDist < distThresh) 
      {
        distance = exp( - sqrDist / expDiv ) - zeroVal;
        // subtract values because mapping vertices move in opposite directions as landmarks
        newAzimuth = atan2( vPos[1], vPos[0] ) - distance*dir[0];
        newPolar = acos( vPos[2] ) - distance*dir[1];
        sinPolar = sin( newPolar );
        (*outputMap)[i][0] = cos( newAzimuth ) * sinPolar;
        (*outputMap)[i][1] = sin( newAzimuth ) * sinPolar;
        (*outputMap)[i][2] = cos( newPolar );
        output->SetParameterizationModified( i, true );
      }
      else 
      {
        (*outputMap)[i] = vPos;
        output->SetParameterizationModified( i, input->GetParameterizationModified( i ) );
      }
    }
    this->GetOutput()->SetBufferedRegion( 0 );
  }


  template <class TParameterization>
  void
  GaussianWarpSphericalParameterizationFilter<TParameterization>
  ::InitializeOutput( ParameterizationType *input )
  {
    ParameterizationPointer output = this->GetOutput();
    (*output) = input;
    // initialize parameterization
    output->InitializeSphericalMap();
  }

}

#endif
