#ifndef _itkParameterizedMeshASCIIReader_h
#define _itkParameterizedMeshASCIIReader_h

#include "vnl/vnl_matrix_fixed.h"
#include "itkMeshSource.h"
#include "itkVector.h"
#include "itkTriangleCell.h"


namespace itk
{

/** \class ParameterizedMeshASCIIReader
 * \brief Reads a spherical parameterized mesh from disk.
 *
 * The mesh data is stored separately in three ASCII files. In addition to the 
 * pts- and fce-file that itk::MeshASCIIReader is using as well, there is a 
 * third file with a par-extension which stores the position of all points in 
 * parameter space. For spherical parameterized meshes, this position is saved 
 * as x, y and z coordinates of the point on the unit sphere, separated by 
 * spaces.
 *
 * \author Tobias Heimann. Division Medical and Biological Informatics, 
 *         German Cancer Research Center, Heidelberg, Germany.
 */
template <class TOutputMesh>
class ParameterizedMeshASCIIReader : public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef ParameterizedMeshASCIIReader           Self;
  typedef MeshSource<TOutputMesh>   Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(ParameterizedMeshASCIIReader, MeshSource);

  /** Some convenient typedefs. */
  typedef TOutputMesh OutputMeshType;
  typedef typename OutputMeshType::MeshTraits   OMeshTraits;
  typedef typename OutputMeshType::PointType    OPointType;
  typedef typename OMeshTraits::PixelType       OPixelType;  
  typedef typename OutputMeshType::Pointer OutputMeshPointer;
  typedef typename OutputMeshType::CellTraits CellTraits;
  typedef typename OutputMeshType::PointsContainerPointer PointsContainerPointer;
  typedef typename OutputMeshType::PointsContainer   PointsContainer;
  typedef CellInterface<OPixelType, CellTraits>   CellInterfaceType;
  typedef TriangleCell<CellInterfaceType>         TriCellType;
  typedef typename TriCellType::SelfAutoPointer       TriCellAutoPointer;
  typedef typename TriCellType::CellAutoPointer       CellAutoPointer;

  /** Set/Get the file prefix (name without extension). */
  itkSetStringMacro(FilePrefix);
  itkGetStringMacro(FilePrefix);

  /** Can be called after Update() to check for errors. */
  bool GetReadError()
  {
    return m_Error;
  }
  
protected:
  ParameterizedMeshASCIIReader();
  ~ParameterizedMeshASCIIReader() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

    
  /** file prefix (name without extension) */
  std::string m_FilePrefix;
  
private:

  ParameterizedMeshASCIIReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool m_Error;

};

} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParameterizedMeshASCIIReader.txx"
#endif
#endif
