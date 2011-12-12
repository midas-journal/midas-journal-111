#ifndef _itkVarianceBasedCostFunction_txx
#define _itkVarianceBasedCostFunction_txx

#include "itkVarianceBasedCostFunction.h"


namespace itk
{

  template <class TShapeModelCalculator>
  VarianceBasedCostFunction<TShapeModelCalculator>
  ::VarianceBasedCostFunction()
  {
    m_MultiThreader = MultiThreader::New();
    m_NumberOfThreads = MultiThreader::GetGlobalDefaultNumberOfThreads();
  }


  template <class TShapeModelCalculator>
  VarianceBasedCostFunction<TShapeModelCalculator>
  ::~VarianceBasedCostFunction()
  {
  }


  template <class TShapeModelCalculator>
  void
  VarianceBasedCostFunction<TShapeModelCalculator>
  ::PrepareGradients( MatrixArray &dataGradients )
  {
    m_EVGradients.resize( this->m_Model->GetNumberOfInputs() );
    std::vector<ThreadDataStruct> ts( m_NumberOfThreads );
    unsigned int numThreads = m_NumberOfThreads;
    m_MultiThreader->SetNumberOfThreads( numThreads );
    unsigned int evIdx = 0;
    do {
      // adjust number of threads if we don't need all
      unsigned int remainingEVs = this->m_Model->GetNumberOfInputs()-1 - evIdx;
      if (remainingEVs < numThreads) {
        numThreads = remainingEVs;
        m_MultiThreader->SetNumberOfThreads( numThreads );
      }
      // init thread structure for all threads
      for (unsigned int i=0; i<numThreads; i++) {
        ts[i].costFunction = this;
        ts[i].index = evIdx+i;
        m_EVGradients[evIdx+i].set_size( 
          this->m_Model->GetNumberOfInputs(), 
          dataGradients[0].columns() );
        ts[i].gradients = &(m_EVGradients[evIdx+i]);
        ts[i].cgArray = &dataGradients;
        m_MultiThreader->SetMultipleMethod( 
          i, this->ComputeEigenValueGradientsThreaderCallback, (void*)&(ts[i]) );
      }
      // execute
      m_MultiThreader->MultipleMethodExecute();
      evIdx += numThreads;
    } while (evIdx < this->m_Model->GetNumberOfInputs()-1);
  }


  template <class TShapeModelCalculator>
  double
  VarianceBasedCostFunction<TShapeModelCalculator>
  ::GetMatrixColumnDotProduct( const MatrixType *m1, unsigned int col1, const MatrixType *m2, unsigned int col2 ) const
  {
    double result = 0.0;
    unsigned int numRows = m1->rows();
    for (unsigned int i=0; i<numRows; i++) { result += (*m1)[i][col1] * (*m2)[i][col2]; }
    return result;
  }

  
  template <class TShapeModelCalculator>
  void
  VarianceBasedCostFunction<TShapeModelCalculator>
  ::ComputeEigenValueGradients( unsigned int ev, MatrixType *evGradients, MatrixArray *dataGradients ) const 
  {
    for (unsigned int sampleIdx=0; sampleIdx<this->m_Model->GetNumberOfInputs(); sampleIdx++) 
    {
      double f = 2.0 * this->GetSingularValue( ev ) * this->GetSVDMatrixV()->get( sampleIdx, ev ); 
      for (unsigned int i=0; i<(*dataGradients)[sampleIdx].columns(); i++) 
      {
        (*evGradients)[sampleIdx][i] = f * this->GetMatrixColumnDotProduct( 
                                  this->GetSVDMatrixU(), ev, &((*dataGradients)[sampleIdx]), i );
      }
    }
  }


  template <class TShapeModelCalculator>
  ITK_THREAD_RETURN_TYPE 
  VarianceBasedCostFunction<TShapeModelCalculator>
  ::ComputeEigenValueGradientsThreaderCallback( void *arg )
  {
    ThreadDataStruct *str = (ThreadDataStruct*)(((MultiThreader::ThreadInfoStruct*)(arg))->UserData);
    dynamic_cast<const Self*>(str->costFunction)->ComputeEigenValueGradients( str->index, str->gradients, str->cgArray );
    return ITK_THREAD_RETURN_VALUE;
  }

}

#endif
