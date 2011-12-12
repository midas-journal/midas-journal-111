#ifndef _itkProcrustesAlign3DMeshFilter_txx
#define _itkProcrustesAlign3DMeshFilter_txx

#include "itkProcrustesAlign3DMeshFilter.h"


namespace itk
{

  template <class TInputMesh, class TOutputMesh>
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::ProcrustesAlign3DMeshFilter()
  {
    m_Convergence = 0.0001;
    m_Mean = OutputMeshType::New();
    m_OldMean = OutputMeshType::New();
    m_AlignRotation = true;
    m_AlignScale = true;
    m_AlignTranslation = true;
    m_UseScaling = true;
    m_UseInitialAverage = false;
    m_UseSingleIteration = false;
    m_MeanScale = 0.0;
  }


  template <class TInputMesh, class TOutputMesh>
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::~ProcrustesAlign3DMeshFilter()
  {
  }

  
  template<class TInputMesh, class TOutputMesh>
  typename ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>::DataObjectPointer
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::MakeOutput(unsigned int /* idx */ )
  {
    return static_cast<DataObject*>(TOutputMesh::New().GetPointer());
  }


  template<class TInputMesh, class TOutputMesh>
  typename ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>::OutputMeshType *
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::GetOutput(unsigned int idx)
  {
    return m_MeshTransform[idx]->GetOutput();
  }


  template <class TInputMesh, class TOutputMesh>
  void
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::SetNumberOfInputs( unsigned int num )
  {
    this->ProcessObject::SetNumberOfInputs( num );
    this->ProcessObject::SetNumberOfRequiredInputs( num );
    this->ProcessObject::SetNumberOfOutputs( num );
    this->ProcessObject::SetNumberOfRequiredOutputs( num );
    m_MeshTransform.resize( num );
    m_Center.resize( num );
    for (unsigned int i=0; i<num; i++)
    {
      OutputMeshPointer output = static_cast<TOutputMesh*>(this->MakeOutput(i).GetPointer()); 
      this->ProcessObject::SetNthOutput( i, output.GetPointer() );
      m_MeshTransform[i] = TransformMeshType::New();
      m_MeshTransform[i]->GraftOutput( this->GetOutput( i ) );
    }
  }


  template <class TInputMesh, class TOutputMesh>
  void
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::SetInput( unsigned int idx, InputMeshPointer mesh )
  {
    this->ProcessObject::SetNthInput( idx, mesh );
    m_MeshTransform[idx]->SetInput( mesh );
  }


  template <class TInputMesh, class TOutputMesh>
  void
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::GenerateData()
  {
    // find mesh centers and store them
    m_MeanCenter.Fill( 0 );
    for (unsigned int i=0; i<this->GetNumberOfInputs(); i++)
    {
      m_Center[i] = GetMeshCenter( this->GetInput( i ) );
      m_MeanCenter += m_Center[i];
    }
    m_MeanCenter /= (double)this->GetNumberOfInputs();
            
    // check if mean and oldMean have to be initialized
    InputMeshPointer mesh = this->GetInput( 0 );
    unsigned int numPoints = mesh->GetNumberOfPoints();
    if (m_Mean->GetNumberOfPoints() != numPoints)
    {
      m_Mean->GetPoints()->resize( numPoints );
      m_OldMean->GetPoints()->resize( numPoints );
    }

    OutputPointsContainer *meanPoints = m_Mean->GetPoints();
    OutputPointsContainer *oldMeanPoints = m_OldMean->GetPoints();
    typename OutputPointsContainer::Iterator meanIt, oldMeanIt;
    if (m_UseInitialAverage) 
    {
      // initialize mean shape to the actual mean shape
      CalculateMean();
    }
    else
    {
      // initialize mean shape to first input mesh
      typename InputMeshType::PointsContainer::ConstIterator meshIt;
      meshIt = mesh->GetPoints()->Begin();
      for (meanIt=meanPoints->Begin(); meanIt!=meanPoints->End(); ++meanIt) 
      { 
        for (int dim=0; dim<3; dim++) { meanIt.Value()[dim] = meshIt.Value()[dim]; }
        ++meshIt;
      }
    }

    // iteratively update mean
    OutputPointType zeroPnt;
    for (int dim=0; dim<3; dim++) { zeroPnt[dim] = 0; }
    CoordRepType diff, squaredDiff;
    do
    {
      // find current transformations to match the mean
      for (unsigned int i=0; i<this->GetNumberOfInputs(); i++)
      {
        TransformPointer transform = GetProcrustesMatch( i, m_Mean );
        m_MeshTransform[i]->SetTransform( transform );
      }
      // copy mean to oldMean and set mean to 0
      oldMeanIt = oldMeanPoints->Begin();
      for (meanIt=meanPoints->Begin(); meanIt!=meanPoints->End(); ++meanIt)
      { 
        oldMeanIt.Value() = meanIt.Value();
        meanIt.Value() = zeroPnt;
        ++oldMeanIt;
      }
      // calculate new mean
      CalculateMean();
      // calculate average point distance between old mean and new mean
      squaredDiff = 0;
      oldMeanIt = oldMeanPoints->Begin();
      for (meanIt=meanPoints->Begin(); meanIt!=meanPoints->End(); ++meanIt)
      { 
        for (int dim=0; dim<3; dim++) 
        {
          diff = meanIt.Value()[dim] - oldMeanIt.Value()[dim];
          squaredDiff += diff*diff;
        }
        ++oldMeanIt;
      }
      diff = sqrt( squaredDiff );
    } while (diff > m_Convergence && !m_UseSingleIteration);
    if (!m_AlignScale && m_UseScaling)
    {
      // scale factors are not set in the transforms, so we have to set them now
      for (unsigned int i=0; i<this->GetNumberOfInputs(); i++)
      {
        if ((ITK_VERSION_MAJOR==2 && ITK_VERSION_MINOR>=2) || ITK_VERSION_MAJOR>2)
        {
          // Starting with ITK 2.2 transform has to prescale
          m_MeshTransform[i]->GetTransform()->Scale( m_MeanScale, true );
        }
        else
        {
          m_MeshTransform[i]->GetTransform()->Scale( m_MeanScale );
        }
        m_MeshTransform[i]->Modified();
        m_MeshTransform[i]->GetOutput()->Update();
      }
    }
    for (unsigned int i=0; i<this->GetNumberOfInputs(); i++)
    {
      this->GetOutput( i )->SetBufferedRegion( this->GetOutput( i )->GetRequestedRegion() );
    }
  }

  
  template <class TInputMesh, class TOutputMesh>
  typename ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>::TranslationType 
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::GetMeshCenter( InputMeshPointer mesh )
  {
    TranslationType center;
    center.Fill( 0 );
    InputPointType pnt;
    typename InputMeshType::PointsContainer *meshPoints = mesh->GetPoints();
    typename InputMeshType::PointsContainer::ConstIterator pntIt;
    for (pntIt=meshPoints->Begin(); pntIt!=meshPoints->End(); ++pntIt)
    {
      for (int dim=0; dim<3; dim++)
      {
        center[dim] += pntIt.Value()[dim];
      }
    }
    center /= (double)meshPoints->Size();
    return center;
  }


  template <class TInputMesh, class TOutputMesh>
  typename ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>::TransformPointer 
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::GetProcrustesMatch( unsigned int idx, OutputMeshPointer targetMesh )
  {
    // copy source mesh coordinates to source matrix, translating to zero
    MatrixType source;
    InputMeshPointer sourceMesh = this->GetInput( idx );
    source.set_size( sourceMesh->GetNumberOfPoints(), 3 );
    typename InputMeshType::PointsContainer::ConstIterator inputIt;
    unsigned int i=0;
    for (inputIt=sourceMesh->GetPoints()->Begin(); inputIt!=sourceMesh->GetPoints()->End(); ++inputIt)
    {
      for (int dim=0; dim<3; dim++) 
      { 
        source[i][dim] = inputIt.Value()[dim] - m_Center[idx][dim]; 
      }
      i++;
    }
    // copy target mesh coordinates to target matrix
    MatrixType target;
    target.set_size( 3, targetMesh->GetNumberOfPoints() );
    typename OutputMeshType::PointsContainer::ConstIterator outputIt;
    i = 0;
    for (outputIt=targetMesh->GetPoints()->Begin(); outputIt!=targetMesh->GetPoints()->End(); ++outputIt)
    {
      for (int dim=0; dim<3; dim++) 
      { 
        target[dim][i] = outputIt.Value()[dim]; 
      }
      i++;
    }
    // do procruste matching
    MatrixType x1 = target * source / (target.fro_norm()*source.fro_norm());
    vnl_svd<CoordRepType> svd( x1 );
    MatrixType postTrans = svd.V() * svd.U().transpose();
    MatrixType x2 = target * source * postTrans;
    CoordRepType x2Trace = 0;
    for (unsigned int i=0; i<x2.rows(); i++) { x2Trace += x2[i][i]; }
    MatrixType x3 = source.transpose() * source;
    CoordRepType x3Trace = 0;
    for (unsigned int i=0; i<x3.rows(); i++) { x3Trace += x3[i][i]; }
    CoordRepType scale = x2Trace / x3Trace;

    // set up transformation
    TransformPointer result = TransformType::New();
    typename TransformType::InputPointType center;
    for (int dim=0; dim<3; dim++) { center[dim] = m_Center[idx][dim]; }
    result->SetCenter( center );

    if (m_AlignTranslation)
    {
      result->SetTranslation( -m_Center[idx] );
    }
    else
    {
      // the results have to lie at zero origin, so use the same translation if
      // shapes should not be aligned individually
      result->SetTranslation( -m_MeanCenter );
    }
    if (m_AlignRotation)
    {
      typename TransformType::MatrixType rotMatrix;
      // change matrix element order because TransformType uses pre-multiply, while we calculated post-multiply 
      for (int r=0; r<3; r++) 
      {
        for (int c=0; c<3; c++) { rotMatrix[r][c] = postTrans[c][r]; }
      }
      result->SetMatrix( rotMatrix );
    }
    if (m_AlignScale && m_UseScaling) 
    { 
      if ((ITK_VERSION_MAJOR==2 && ITK_VERSION_MINOR>=2) || ITK_VERSION_MAJOR>2)
      {
        // Starting with ITK 2.2 transform has to prescale
        result->Scale( scale, true ); 
      }
      else
      {
        result->Scale( scale ); 
      }
    }
    return result;
  }


  template <class TInputMesh, class TOutputMesh>
  void
  ProcrustesAlign3DMeshFilter<TInputMesh, TOutputMesh>
  ::CalculateMean()
  {
    typename OutputMeshType::PointsContainer *meanPoints = m_Mean->GetPoints();
    typename OutputMeshType::PointsContainer::Iterator meanIt;
    double meanWeight = 1.0 / (double)(this->GetNumberOfInputs());
    for (unsigned int i=0; i < this->GetNumberOfInputs(); i++)
    {
      // For some weird reason, we have to call Update() explicitly here...
      m_MeshTransform[i]->Update();
      OutputMeshPointer transformedMesh = m_MeshTransform[i]->GetOutput();
      typename OutputMeshType::PointsContainer *transformedPoints = transformedMesh->GetPoints();
      typename OutputMeshType::PointsContainer::ConstIterator pntIt;
      pntIt=transformedPoints->Begin();
      for (meanIt = meanPoints->Begin(); meanIt!=meanPoints->End(); ++meanIt)
      { 
        for (int dim=0; dim<3; dim++) 
        {
          meanIt.Value()[dim] += meanWeight * pntIt.Value()[dim];
        }
        ++pntIt;
      }
    }
    // scale mean to get a norm of 1
    double squaredNorm = 0;
    for (meanIt = meanPoints->Begin(); meanIt!=meanPoints->End(); ++meanIt)
    {
      squaredNorm += meanIt.Value().GetVectorFromOrigin().GetSquaredNorm();
    }
    m_MeanScale = 1.0 / sqrt( squaredNorm );
    for (meanIt = meanPoints->Begin(); meanIt!=meanPoints->End(); ++meanIt)
    { 
      for (int dim=0; dim<3; dim++) 
      {
        meanIt.Value()[dim] *= m_MeanScale; 
      }
    }
    if (m_AlignScale)
    {
      // m_MeanScale has to be adjusted by the scalings used in the transformations
      double meanTransformScale = 0.0;
      for (unsigned int i=0; i<this->GetNumberOfInputs(); i++)
      {
        meanTransformScale += m_MeshTransform[i]->GetTransform()->GetScale()[0];
      }
      m_MeanScale *= meanTransformScale * meanWeight;
    }
  }

}

#endif
