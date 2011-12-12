#ifndef _itkParameterizedTriangleMesh_txx
#define _itkParameterizedTriangleMesh_txx

#include "itkParameterizedTriangleMesh.h"


namespace itk
{

  template <typename TPixelType, typename TIndex, typename TCoordRepType, typename TMeshTraits>
  typename ParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType, TMeshTraits>::Pointer
  ParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType, TMeshTraits>
  ::operator=( IndexedTriangleMeshConstPointer mesh )
  {
    // copy Mesh data

    // copy IndexedTriangleMesh data

    return this;
  }

}

#endif