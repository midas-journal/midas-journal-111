
#ifndef _itkMeshDataSpecialization_h_
#define _itkMeshDataSpecialization_h_

#include <metaMesh.h>
#include <itkVector.h>

typedef itk::Vector<double, 2> VectorType ;


template<>
class MeshData <VectorType>: public MeshDataBase
{
public:

  MeshData() {m_Id=-1;}
  ~MeshData() {};

  virtual MET_ValueEnumType GetMetaType()
    {
    return MET_DOUBLE_ARRAY;
    }

  virtual void Write( std::ofstream* /* stream */ )
    {
    }

  virtual unsigned int GetSize(void)
    {
     unsigned int size = sizeof(int);
     size += sizeof(m_Data);
     return size;
    }

  VectorType m_Data;
};

#endif
