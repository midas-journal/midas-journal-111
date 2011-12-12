  #ifndef _itkStatisticalShapeModel3DCalculatorWithFeatures_txx
#define _itkStatisticalShapeModel3DCalculatorWithFeatures_txx

#include "itkStatisticalShapeModel3DCalculatorWithFeatures.h"
#include "itkStatisticalShapeModel3DCalculatorWithFeatures.h"
#include <algorithm>


namespace itk
{
  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  StatisticalShapeModel3DCalculatorWithFeatures<TLandmarks, TParameterization, TOutputMeshes>
  ::StatisticalShapeModel3DCalculatorWithFeatures()  
  {
  }

  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  StatisticalShapeModel3DCalculatorWithFeatures<TLandmarks, TParameterization, TOutputMeshes>
  ::~StatisticalShapeModel3DCalculatorWithFeatures()
  {
  }

  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculatorWithFeatures<TLandmarks, TParameterization, TOutputMeshes>
  ::PrintSelf( std::ostream& os, Indent indent ) const
  {
    Superclass::PrintSelf(os, indent);
  }

  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculatorWithFeatures<TLandmarks, TParameterization, TOutputMeshes>
  ::SetLandmarks( LandmarkPointer landmarks )
  {
    typename LandmarkType::PointDataContainer::Pointer temp = LandmarkType::PointDataContainer::New(); 
    landmarks->SetPointData (temp) ;
    this->m_LandmarkRotator->SetInput( landmarks );
  }

  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void 
  StatisticalShapeModel3DCalculatorWithFeatures<TLandmarks, TParameterization, TOutputMeshes> 
  ::ComputeFeatureMeans() const
  {
    // Compute the mean for each feature 
    unsigned int nDims = this->GetNumberOfComponents() ;
    unsigned int nSubjects = this->GetNumberOfInputs() ;
    unsigned int nPoints = this->GetNumberOfLandmarks() ;

    unsigned int sub, pts ;
    FeatureVectorType sum, data ;
    
    for (unsigned int i = 0 ; i < nDims ; i++ )
    {
      sum[i] = 0 ;
    }

    for ( sub = 0 ; sub < nSubjects ; sub++ )
    {
      OutputMeshPointer mesh = this->m_Transforms[sub]->GetOutput() ;
      for ( pts = 0 ; pts < nPoints ; pts++ )
      {
        mesh->GetPointData(pts, &data) ;
        for ( unsigned int dim = 0 ; dim < nDims ; dim ++ ) 
        {
          sum[dim] += data[dim] ;
        }
      }
    }

    //cout << sum[0] << " " << sum[1] << endl ;
  
    for ( unsigned int d = 0 ; d < nDims ; d++ ) 
    {
      this->m_FeatureMeans[d] = sum[d] / ( nSubjects * nPoints ) ;
      //cout << m_FeatureMeans[d] << " " ;
    }
  }

  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculatorWithFeatures<TLandmarks, TParameterization, TOutputMeshes>
  ::InitializeMatrix( MatrixType *matrix, unsigned int sampleIdx, OutputMeshPointer mesh ) const
  {
    unsigned int nDim = this->GetNumberOfComponents() ;
    unsigned int row, i ;
        
    this->DebugHelper () ;
    mesh->Update();
    this->DebugHelper () ;
    this->ComputeFeatureMeans() ;

    FeatureVectorType sample ;
    
    for (row=0, i = 0; row<matrix->rows(); row+= nDim, i++)
    {
      mesh->GetPointData ( i, &sample ) ;
       
      for ( unsigned int dim=0 ; dim < nDim ; dim++ )
      {
        (*matrix)[row+dim][sampleIdx] = sample[dim] - this->m_FeatureMeans[dim];
      }
    }
  }
}

#endif
