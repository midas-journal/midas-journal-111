#ifndef ITK_MESH_ASCII_READER_INCLUDED
#define ITK_MESH_ASCII_READER_INCLUDED

#include "vnl/vnl_matrix_fixed.h"
#include "itkMeshSource.h"
#include "itkVector.h"
#include "itkTriangleCell.h"


namespace itk
{

/** \class MeshASCIIReader
 *  \brief Reads a mesh from disk.
 *
 * The mesh data is stored separately in two ASCII files. The first file has a 
 * pts-extension and holds the vertices of the mesh: Each line represents one 
 * point by the x, y and z coordinates, separated by spaces. The second file 
 * has a fce-extension and stores the faces of the mesh: Each line represents 
 * one face by the indices of the involved vertices, separated by spaces. 
 * The first point in the pts-file has the index 0.
 * Currently, only triangular faces are supported by the reader. After 
 * specifying the file prefix (the name without the extension) with 
 * SetFilePrefix() and calling Update(), GetReadError() can be used to query 
 * the success of the operation. If no error occurred, GetOutput() returns the 
 * imported mesh.
 *
 * \author Tobias Heimann. Division Medical and Biological Informatics, 
 *         German Cancer Research Center, Heidelberg, Germany.
 */
template <class TOutputMesh>
class MeshASCIIReader : public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef MeshASCIIReader           Self;
  typedef MeshSource<TOutputMesh>   Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshASCIIReader, MeshSource);

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
  MeshASCIIReader();
  ~MeshASCIIReader() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

    
  /** file prefix (name without extension) */
  std::string m_FilePrefix;
  
  
private:
  MeshASCIIReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool m_Error;

};

} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshASCIIReader.txx"
#endif
#endif
