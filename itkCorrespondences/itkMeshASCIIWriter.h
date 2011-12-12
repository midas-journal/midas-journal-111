#ifndef MESHASCIIWRITER_H_HEADER_INCLUDED
#define MESHASCIIWRITER_H_HEADER_INCLUDED

#include "itkMeshFileWriter.h"

namespace itk {


  /** \class MeshASCIIWriter
  *   \brief Stores a mesh in ASCII format.
  *
  * The file format is described in the documentation for the 
  * itk::MeshASCIIReader class. Currently, only triangular faces are stored. 
  * Use SetFilePrefix() to specify the filename, which will be extended with 
  * the pts- and fce-extensions.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
  *         German Cancer Research Center, Heidelberg, Germany.
  */
template <class TMeshType>
class MeshASCIIWriter : public MeshFileWriter<TMeshType>
{
  public:

    typedef MeshASCIIWriter<TMeshType> Self;
    typedef MeshFileWriter<TMeshType> Superclass;
    typedef SmartPointer < Self > Pointer;
    typedef SmartPointer < const Self > ConstPointer;

    itkTypeMacro( MeshASCIIWriter<TMeshType>, MeshFileWriter<TMeshType> );
    itkNewMacro( Self );  


    virtual void GenerateData()
    {
      if (strlen( this->m_filePrefix ) == 0) {
        itkExceptionMacro(<<"GenerateData:Please specify a file prefix!");
        return;
      }
      FILE *file;
      char fileName[1024];
      typedef typename TMeshType::CellsContainer::ConstIterator CellIterator;
      typedef typename TMeshType::CellType CellType;
      typename TMeshType::Pointer mesh = dynamic_cast<TMeshType*> (this->GetInput());
      // write triangle indices:
      sprintf( fileName, "%s.fce", this->m_filePrefix );
      file = fopen( fileName, "w" );
      if (file) {
        CellIterator cellIterator = mesh->GetCells()->Begin();
        CellIterator cellEnd      = mesh->GetCells()->End();
        typename TMeshType::CellIdentifier cId = 0;
        while( cellIterator != cellEnd ) {
          CellType *cell = cellIterator.Value();
          if (cell->GetType() == CellType::TRIANGLE_CELL) {
            typename CellType::PointIdConstIterator pntId = cell->GetPointIds();
            fprintf( file, "%u %u %u\n", (unsigned int)pntId[0], (unsigned int)pntId[1], (unsigned int)pntId[2] );
          }
          cId++;
          cellIterator++;
        }
        fclose( file );
      }
      // write point coordinates:
      sprintf( fileName, "%s.pts", this->m_filePrefix );
      file = fopen( fileName, "w" );
      if (file) {
        typename TMeshType::PointType meshPoint;
	      for (int v=0; v < mesh->GetNumberOfPoints(); v++) {
		      mesh->GetPoint( v, &meshPoint );
		      fprintf( file, "%f %f %f\n", meshPoint[0], meshPoint[1], meshPoint[2] );
        }
        fclose( file );
      }
    }

   
protected:
  MeshASCIIWriter() {};
  virtual ~MeshASCIIWriter() {};

};

} // namespace itk
#endif /* MESHASCIIWRITER_H_HEADER_INCLUDED */
