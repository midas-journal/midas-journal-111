#ifndef _itkStatisticalShapeModel3DCalculator_txx
#define _itkStatisticalShapeModel3DCalculator_txx

#include "itkStatisticalShapeModel3DCalculator.h"
#include "itkShapeModelCalculatorCostFunction.h"
#include <algorithm>


namespace itk
{

  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::StatisticalShapeModel3DCalculator() 
  {
    m_NumberOfInputs = 0;
    m_OptimizeParameterizationStart = true;
    m_ParameterizationWarpStepLength = 0.00002;//0.000015;
    m_ParameterizationStartStepLength = 0.00004;
    m_UseTangentScale = true;
    m_Convergence = 0.01;
    m_LandmarkRotator = TransformMeshType::New();
    m_LandmarkRotator->ReleaseDataBeforeUpdateFlagOff();
    m_Landmarks = m_LandmarkRotator->GetOutput();
    m_OutputLandmarks = 0 ;
    m_ProcrustesFilter = ProcrustesAlignFilterType::New();
    m_MultiThreader = MultiThreader::New();
    m_NumberOfThreads = MultiThreader::GetGlobalDefaultNumberOfThreads();
    m_CostFunction = 0;
    m_ResampleResult = 0 ;
    srand( 0 );
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::~StatisticalShapeModel3DCalculator()
  {
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::PrintSelf( std::ostream& os, Indent indent ) const
  {
    Superclass::PrintSelf(os, indent);
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::SetNumberOfInputs( unsigned int num ) 
  {
    int nDim = this->GetNumberOfComponents() ;
    m_NumberOfInputs = num;
    m_ProcrustesFilter->SetNumberOfInputs( num );
    m_Parameterizations.resize( num );
    m_ParameterizationRotator.resize( num );
    m_ReparameterizationFilters.resize( num );
    m_ShapeGenerators.resize( num );
    m_Transforms.resize( num );
    m_GradientShapeGenerators.resize( num );
    m_GradientTransforms.resize( num );
    m_SingularValues.resize( num ); 
    m_EigenValues.resize( num ); 
    for (unsigned int i=0; i<num; i++)
    {
      m_Parameterizations[i] = ParameterizationType::New();
      m_ParameterizationRotator[i] = RotateParameterizationFilterType::New();
      m_ReparameterizationFilters[i] = ReparameterizationFilterType::New();
      m_ShapeGenerators[i] = RemeshingType::New();
      m_ShapeGenerators[i]->SetLandmarks( m_Landmarks );
      m_ShapeGenerators[i]->SetPointDataDimension (nDim) ;
      m_ProcrustesFilter->SetInput( i, m_ShapeGenerators[i]->GetOutput() );
      m_Transforms[i] = TransformMeshType::New();
      m_Transforms[i]->SetInput( m_ProcrustesFilter->GetOutput( i ) );
      m_GradientShapeGenerators[i] = RemeshingType::New();
      m_GradientShapeGenerators[i]->SetLandmarks( m_Landmarks );
      m_GradientShapeGenerators[i]->SetPointDataDimension (nDim) ;
      m_GradientShapeGenerators[i]->SetInput( m_ReparameterizationFilters[i]->GetOutput() );
      m_GradientTransforms[i] = TransformMeshType::New();
      m_GradientTransforms[i]->SetInput( m_GradientShapeGenerators[i]->GetOutput() );
    }
    this->Modified();
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::SetInput( unsigned int idx, ParameterizationPointer parameterization )
  {
    // copy input parameterization:
    parameterization->Update();
    *(m_Parameterizations[idx]) = parameterization.GetPointer();
    m_Parameterizations[idx]->InitializeSphericalMap();
    m_Parameterizations[idx]->UpdateParameterization( parameterization.GetPointer() );
    // connect to pipeline
    m_ParameterizationRotator[idx]->SetInput( m_Parameterizations[idx] );
    m_ReparameterizationFilters[idx]->SetInput( m_Parameterizations[idx] );
    m_ShapeGenerators[idx]->SetInput( m_Parameterizations[idx] );
    this->Modified();
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::SetLandmarks( LandmarkPointer landmarks )
  {
    m_LandmarkRotator->SetInput( landmarks );
  }

  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::SetOutputLandmarks( LandmarkPointer landmarks )
  {
    m_OutputLandmarks = landmarks ;
  }

  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::GenerateData()
  {
    if (m_CostFunction == 0)
    {
      itkExceptionMacro( "Cost Function not set!" );
    }
    m_Landmarks->Update();
    m_Landmarks->DisconnectPipeline();
    m_LandmarkRotator->SetInput( m_Landmarks );

    unsigned int numRows = this->GetNumberOfLandmarks() * this->GetNumberOfComponents();
    unsigned int numCols = this->GetNumberOfInputs();
    if ( m_AlignedData.rows() != numRows || m_AlignedData.columns() != numCols )
    {
      m_AlignedData.set_size( numRows, numCols );
    }
    
    this->DebugHelper();     
    this->BuildModel();
    this->DebugHelper() ;
    double delta = m_Convergence + 1.0;
    double oldCosts = m_CostFunction->GetValue();
    m_CurrentNumberOfIterations = 0;
    std::cout << "Initial costs = " << m_CostFunction->GetValue() << std::endl;
    if (m_OptimizeParameterizationStart)
    {
      std::cout << "Optimizing initial orientation." << std::endl;
      do 
      {
        m_CurrentNumberOfIterations++;
        this->OptimizeParameterizationStarts();
        this->BuildModel();
        double costs = m_CostFunction->GetValue();
        std::cout << "Iteration " << m_CurrentNumberOfIterations << ": costs = " << costs << std::endl;
        if (m_CurrentNumberOfIterations % 10 == 0)
        {
          delta = oldCosts - costs;
          oldCosts = costs;
        }
      } while (delta > m_Convergence);
      std::cout << "Finished optimizing orientation." << std::endl;
      delta = m_Convergence + 1.0;
      oldCosts += delta;
    }
    this->DebugHelper() ;
    do 
    {
      m_CurrentNumberOfIterations++;
      this->OptimizeParameterizations();
      if (m_OptimizeParameterizationStart)
      {
        this->BuildModel();
        this->OptimizeParameterizationStarts();
      }
      this->BuildModel();
      this->RotateParameterizationsAndLandmarks();
      double costs = m_CostFunction->GetValue();
      std::cout << "Iteration " << m_CurrentNumberOfIterations << ": costs = " << costs << std::endl;
      if (m_CurrentNumberOfIterations % 50 == 0)
      {
        delta = oldCosts - costs;
        oldCosts = costs;
        if (delta <= m_Convergence)
        {
          // model didn't improve enough, so we check if we can increase the level of detail
          unsigned int lod = m_ReparameterizationFilters[0]->GetLevelOfDetail();
          if (lod < m_ReparameterizationFilters[0]->GetMaximumLevelOfDetail())
          {
            lod++;
            std::cout << std::endl << "Switching to level of detail " << lod << std::endl;
            for (unsigned int i=0; i<m_NumberOfInputs; i++)
            {
              m_ReparameterizationFilters[i]->SetLevelOfDetail( lod );
            }
            // keep on going with new level of detail
            delta = m_Convergence + 1.0;
          }
        }
      }
    } while (delta > m_Convergence);
    this->DebugHelper() ;
    this->BuildModel();
    this->DebugHelper() ;
    std::cout << "Finished optimization at iteration " << m_CurrentNumberOfIterations << std::endl;
    std::cout << "Final costs = " << m_CostFunction->GetValue() << std::endl;
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::OptimizeParameterizations()
  {
    MatrixArray dataGradients( m_NumberOfInputs );
    this->ComputeAllReparameterizationGradients( dataGradients );
    m_CostFunction->PrepareGradients( dataGradients );
    // update landmarks:
    for (unsigned int sampleIdx=0; sampleIdx<m_NumberOfInputs; sampleIdx++) 
    {
      unsigned int numControlPoints = m_ReparameterizationFilters[sampleIdx]->GetNumberOfControlPoints();
      for (unsigned int cp=0; cp<numControlPoints; cp++) 
      {
        VectorType moveDir = m_CostFunction->GetGradientVector( sampleIdx, cp );
        m_ReparameterizationFilters[sampleIdx]->SetActiveControlPoint( cp );
        m_ReparameterizationFilters[sampleIdx]->SetDirection( moveDir * m_ParameterizationWarpStepLength );
        m_Parameterizations[sampleIdx]->UpdateParameterization( 
          m_ReparameterizationFilters[sampleIdx]->GetOutput() );
      }
    }
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::OptimizeParameterizationStarts()
  {
    MatrixArray dataGradients( m_NumberOfInputs );
    this->ComputeAllParameterizationStartGradients( dataGradients );
    m_CostFunction->PrepareGradients( dataGradients );
    // update landmarks:
    typename TransformType::OutputVectorType axisVec;
    typename TransformType::Pointer transform;
    for (unsigned int sampleIdx=0; sampleIdx<m_NumberOfInputs; sampleIdx++) 
    {
      for (unsigned int axis=0; axis<3; axis++) 
      {
        double move = m_CostFunction->GetGradient( sampleIdx, axis );
        axisVec.Fill( 0 );
        axisVec[axis]=1;
        transform = TransformType::New();
        transform->Rotate3D( axisVec, m_ParameterizationStartStepLength*move );
        m_ParameterizationRotator[sampleIdx]->SetTransform( transform );
        m_ParameterizationRotator[sampleIdx]->Update();
        m_Parameterizations[sampleIdx]->UpdateParameterization( m_ParameterizationRotator[sampleIdx]->GetOutput() );
      }
    }
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::BuildModel()
  {
    m_ProcrustesFilter->Update();
    this->DebugHelper() ;
    double scaleFactor = 1.0;
    vnl_vector<double> mean( m_AlignedData.rows() );
    mean.fill( 0 );
    for (unsigned int i=0; i<this->GetNumberOfInputs(); i++)
    {
      if (m_UseTangentScale) { scaleFactor = 1.0 / this->GetTangentSize( i ); }
      m_Transforms[i]->GetTransform()->SetIdentity();
      m_GradientTransforms[i]->SetTransform( m_ProcrustesFilter->GetTransform( i ) );
      if ((ITK_VERSION_MAJOR==2 && ITK_VERSION_MINOR>=2) || ITK_VERSION_MAJOR>2)
      {
        // Starting with ITK 2.2 transform has to prescale
        m_Transforms[i]->GetTransform()->Scale( scaleFactor, true );
        m_GradientTransforms[i]->GetTransform()->Scale( scaleFactor, true );
      }
      else
      {
        m_Transforms[i]->GetTransform()->Scale( scaleFactor );
        m_GradientTransforms[i]->GetTransform()->Scale( scaleFactor );
      }
      this->DebugHelper() ;
      this->InitializeMatrix( &m_AlignedData, i, m_Transforms[i]->GetOutput() );
      mean += m_AlignedData.get_column( i );
      // Cache in the ParameterizedRemeshing filters is up-to-date now
      m_Parameterizations[i]->SetParameterizationModified( false );
    }
    mean /= (double)this->GetNumberOfInputs(); 
    
    MatrixType pcaMatrix( m_AlignedData.rows(), m_AlignedData.cols() );
    double c = 1.0 / sqrt( (double)(this->GetNumberOfInputs()-1) );
    for (unsigned int i=0; i<this->GetNumberOfInputs(); i++)
    {
      pcaMatrix.set_column( i, c * (m_AlignedData.get_column( i ) - mean) );
    }
    vnl_svd<double> svd( pcaMatrix );
    for (unsigned int i=0; i<this->GetNumberOfInputs()-1; i++) 
    {
      double sv = svd.W( i );
      m_SingularValues[i] = sv;
      m_EigenValues[i] = sv*sv;
    }
    m_SVDMatrixU = svd.U();
    m_SVDMatrixV = svd.V();
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  double
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::GetTangentSize( unsigned int idx ) const
  {
    const typename OutputMeshType::PointsContainer *meanPoints = m_ProcrustesFilter->GetMean()->GetPoints();
    const typename OutputMeshType::PointsContainer *samplePoints = m_ProcrustesFilter->GetOutput( idx )->GetPoints();
    typename OutputMeshType::PointsContainer::ConstIterator sample, mean;
    double tangentSize = 0;
    sample = samplePoints->Begin();
    for (mean=meanPoints->Begin(); mean!=meanPoints->End(); ++mean)
    { 
      tangentSize += mean.Value()[0] * sample.Value()[0]
                  +  mean.Value()[1] * sample.Value()[1]
                  +  mean.Value()[2] * sample.Value()[2];
      ++sample;
    }
    return tangentSize;
  }

  
  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::ComputeReparameterizationGradients( unsigned int sampleIdx, MatrixType *gradients ) const
  {
    unsigned int numControlPoints = m_ReparameterizationFilters[sampleIdx]->GetNumberOfControlPoints();
    double delta = 0.005;
    VectorType dir;
    MatrixType modifiedShape( m_AlignedData.rows(), 4 );

    for (unsigned int cp=0; cp<numControlPoints; cp++) 
    {
      m_ReparameterizationFilters[sampleIdx]->SetActiveControlPoint( cp );
      m_GradientShapeGenerators[sampleIdx]->CloneCache( m_ShapeGenerators[sampleIdx] );
      // modify control point to estimate derivative:
      dir[0] = delta;  dir[1] = 0;
      m_ReparameterizationFilters[sampleIdx]->SetDirection( dir );
      this->InitializeMatrix( &modifiedShape, 0, m_GradientTransforms[sampleIdx]->GetOutput() );
      dir[0] = -delta;  dir[1] = 0;
      m_ReparameterizationFilters[sampleIdx]->SetDirection( dir );
      this->InitializeMatrix( &modifiedShape, 1, m_GradientTransforms[sampleIdx]->GetOutput() );
      dir[0] = 0;  dir[1] = delta;
      m_ReparameterizationFilters[sampleIdx]->SetDirection( dir );
      this->InitializeMatrix( &modifiedShape, 2, m_GradientTransforms[sampleIdx]->GetOutput() );
      dir[0] = 0;  dir[1] = -delta;
      m_ReparameterizationFilters[sampleIdx]->SetDirection( dir );
      this->InitializeMatrix( &modifiedShape, 3, m_GradientTransforms[sampleIdx]->GetOutput() );
      // store gradients:
      for (unsigned int g=0; g<4; g++) 
      {
        for (unsigned int i=0; i<m_AlignedData.rows(); i++) 
        {
          (*gradients)[i][4*cp+g] = (modifiedShape[i][g] - m_AlignedData[i][sampleIdx]) / delta;
        }
      }
    }
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  ITK_THREAD_RETURN_TYPE
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::ComputeReparameterizationGradientsThreaderCallback( void *arg )
  {
    ThreadDataStruct *str = (ThreadDataStruct*)(((MultiThreader::ThreadInfoStruct*)(arg))->UserData);
    str->modelCalculator->ComputeReparameterizationGradients( str->index, str->gradients );
    return ITK_THREAD_RETURN_VALUE;
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::ComputeAllReparameterizationGradients( MatrixArray &contourGradient ) const
  {
    std::vector<ThreadDataStruct> ts( m_NumberOfThreads );
    unsigned int numThreads = m_NumberOfThreads;
    m_MultiThreader->SetNumberOfThreads( numThreads );
    unsigned int sampleIdx = 0;
    do 
    {
      // adjust number of threads if we don't need all
      unsigned int remainingSamples = m_NumberOfInputs - sampleIdx;
      if (remainingSamples < numThreads) 
      {
        numThreads = remainingSamples;
        m_MultiThreader->SetNumberOfThreads( numThreads );
      }
      // init thread structure for all threads
      for (unsigned int i=0; i<numThreads; i++) 
      {
        ts[i].modelCalculator = this;
        ts[i].index = sampleIdx+i;
        contourGradient[sampleIdx+i].set_size( 
          this->GetNumberOfLandmarks()*this->GetNumberOfComponents(), 
          4 * (m_ReparameterizationFilters[sampleIdx]->GetNumberOfControlPoints()) );
        ts[i].gradients = &(contourGradient[sampleIdx+i]);
        m_MultiThreader->SetMultipleMethod( 
          i, this->ComputeReparameterizationGradientsThreaderCallback, (void*)&(ts[i]) );
      }
      // execute
      m_MultiThreader->MultipleMethodExecute();
      sampleIdx += numThreads;
    } while (sampleIdx < m_NumberOfInputs);
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::ComputeAllParameterizationStartGradients( MatrixArray &contourGradient ) const
  {
    typename TransformType::OutputVectorType axisVec;
    typename TransformType::Pointer transform;
    double delta = 0.002;   // ~0.5 degrees
    unsigned int numSamples = m_NumberOfInputs;
    unsigned int numRows = this->GetNumberOfLandmarks() * this->GetNumberOfComponents();
    unsigned int numCols = 6; // 3 axes in two directions
    contourGradient.resize( numSamples );
    MatrixType modifiedShape( numRows, 2 );

    for (unsigned int sampleIdx=0; sampleIdx<numSamples; sampleIdx++) 
    {
      m_GradientShapeGenerators[sampleIdx]->SetInput( m_ParameterizationRotator[sampleIdx]->GetOutput() );
      m_GradientShapeGenerators[sampleIdx]->CloneCache( m_ShapeGenerators[sampleIdx] );
      m_ParameterizationRotator[sampleIdx]->GetOutput()->SetParameterizationModified( true );
      (contourGradient[sampleIdx]).set_size( numRows, numCols );
      for (unsigned int axis=0; axis<3; axis++) 
      {
        // build rotation transform around specific axis
        axisVec.Fill( 0 );
        axisVec[axis]=1;
        transform = TransformType::New();
        transform->Rotate3D( axisVec, delta );
        m_ParameterizationRotator[sampleIdx]->SetTransform( transform );
        this->InitializeMatrix( &modifiedShape, 0, m_GradientTransforms[sampleIdx]->GetOutput() );
        // rotate in other direction:
        transform = TransformType::New();
        transform->Rotate3D( axisVec, -delta );
        m_ParameterizationRotator[sampleIdx]->SetTransform( transform );
        this->InitializeMatrix( &modifiedShape, 1, m_GradientTransforms[sampleIdx]->GetOutput() );
        // store gradients:
        for (unsigned int g=0; g<2; g++) 
        {
          for (unsigned int i=0; i<modifiedShape.rows(); i++) 
          {
            double m = modifiedShape[i][g];
            double a = m_AlignedData[i][sampleIdx];
            double cg = m - a;
            (contourGradient[sampleIdx])[i][2*axis+g] = cg / delta;
          }
        }
      }
      m_GradientShapeGenerators[sampleIdx]->SetInput( m_ReparameterizationFilters[sampleIdx]->GetOutput() );
    }
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::RotateParameterizationsAndLandmarks()
  {
    // set up a random rotation matrix
    const double twoPi = 6.2831853;
    double a1 = (double)rand() / (double)RAND_MAX;
    double a2 = (double)rand() / (double)RAND_MAX;
    double a3 = (double)rand() / (double)RAND_MAX;

    Matrix<double,3,1> v;
    v[0][0] = cos( twoPi*a1 ) * sqrt( a2 );
    v[1][0] = sin( twoPi*a1 ) * sqrt( a2 );
    v[2][0] = sqrt( 1 - a2 );

    Matrix<double,3,3> h, x, m, rot;
    x = (v * v.GetTranspose()) * 2.0;
    h.SetIdentity();

    h = h - x;
    m[0][0] = cos( twoPi*a3 );  m[0][1] = sin( twoPi*a3 );  m[0][2] = 0;
    m[1][0] =-sin( twoPi*a3 );  m[1][1] = cos( twoPi*a3 );  m[1][2] = 0;
    m[2][0] = 0;                m[2][1] = 0;                m[2][2] = 1;
    rot = h*m;
    rot *= -1.0;
    
    TransformPointer randomRotation = TransformType::New();
    randomRotation->SetMatrix( rot );

    // rotate landmarks
    m_LandmarkRotator->SetTransform( randomRotation );
    m_LandmarkRotator->Update();
    // copy output to m_landmarks
    typename LandmarkType::PointsContainer *outputPoints = m_LandmarkRotator->GetOutput()->GetPoints();
    typename LandmarkType::PointsContainer *landmarkPoints = m_Landmarks->GetPoints();
    for (unsigned int i=0; i<outputPoints->size(); i++) 
    {
      (*landmarkPoints)[i] = (*outputPoints)[i];
    }
    // rotate parameterizations
    for (unsigned int i=0; i<m_NumberOfInputs; i++)
    {
      m_ParameterizationRotator[i]->SetTransform( randomRotation );
      m_ParameterizationRotator[i]->Update();
      m_Parameterizations[i]->UpdateParameterization( m_ParameterizationRotator[i]->GetOutput() );
      m_Parameterizations[i]->SetParameterizationModified( false );
    }
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::InitializeMatrix( MatrixType *matrix, unsigned int sampleIdx, OutputMeshPointer mesh ) const
  {
    mesh->Update();
    typename OutputMeshType::PointsContainer *meshPoints = mesh->GetPoints();
    typename OutputMeshType::PointsContainer::ConstIterator sample = meshPoints->Begin();
    for (unsigned int row=0; row<matrix->rows(); row+=3)
    {
      (*matrix)[row  ][sampleIdx] = sample.Value()[0];
      (*matrix)[row+1][sampleIdx] = sample.Value()[1];
      (*matrix)[row+2][sampleIdx] = sample.Value()[2];
      ++sample;
    }
  }


  template <class TLandmarks, class TParameterization, class TOutputMeshes>
  void
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::DebugHelper() const
  {
    /*
    int i ;
    itk::Vector < double, 2 > data ;
    
    for ( i = 0 ; i < 0 ; i++ )
    {
      cout << i << endl ;
      m_Parameterizations[i]->GetPointData(0, &data) ;
      cout << "m_Parametrizations " << data[0] << " " << data[1] << endl ;

      m_Transforms[i]->GetOutput()->GetPointData(0, &data) ;
      cout << "m_Transforms " << data[0] << " " << data[1] << endl ;
    
      m_ProcrustesFilter->GetOutput(i)->GetPointData(0, &data) ;
      cout << "m_Procrustes " << data[0] << " " << data[1] << endl ;
    
      m_GradientTransforms[i]->GetOutput()->GetPointData(0, &data) ;
      cout << "m_GradientTransforms " << data[0] << " " << data[1] << endl ;
    
      m_ShapeGenerators[i]->GetOutput()->GetPointData(0, &data) ;
      cout << "m_ShapeGenerators " << data[0] << " " << data[1] << endl ;
    
      cout << endl ;
    }*/
  }

template <class TLandmarks, class TParameterization, class TOutputMeshes>
 typename  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>::OutputMeshPointer
  StatisticalShapeModel3DCalculator<TLandmarks, TParameterization, TOutputMeshes>
  ::GetResampledOutputMesh( unsigned int idx ) 
  {
    if ( m_OutputLandmarks.IsNull() )
    {
      return m_ShapeGenerators[idx]->GetOutput() ;
    }
    
    RemeshingPointer remeshFilter ; 
    OutputMeshPointer mesh ;
    ParameterizationPointer parMesh ;

    remeshFilter = RemeshingType::New() ;
    remeshFilter->SetLandmarks ( m_OutputLandmarks ) ;
    remeshFilter->SetPointDataDimension ( this->GetNumberOfComponents() ) ;
  
    mesh = this->GetOutputMesh ( idx ) ;
    mesh->BuildCellLinks () ;
    parMesh = ParameterizationType::New () ;
    *(parMesh) = mesh.GetPointer();
    parMesh->InitializeSphericalMap () ;
    parMesh->UpdateParameterization((typename ParameterizationType::Superclass::Pointer)m_Parameterizations[idx]);
  
    remeshFilter->SetInput ( parMesh ) ;
    remeshFilter->Update () ;
    

    this->m_ResampleResult = remeshFilter->GetOutput() ;
    
    remeshFilter->Delete () ;
    parMesh->Delete () ;

    return this->m_ResampleResult ;
  }
}

#endif
