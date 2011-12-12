#ifndef PARAMETERIZEDMESHASCIIWRITER_H_HEADER_INCLUDED
#define PARAMETERIZEDMESHASCIIWRITER_H_HEADER_INCLUDED

#include "itkMeshFileWriter.h"

namespace itk {


  /** \class ParameterizedMeshASCIIWriter
  *   \brief Stores a spherical parameterized mesh in ASCII format.
  *
  * The file format is described in the documentation for the 
  * itk::ParameterizedMeshASCIIReader class.
  * Use SetFilePrefix() to specify the filename, which will be extended with 
  * the pts-, fce- and par-extensions.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
  *         German Cancer Research Center, Heidelberg, Germany.
  */
template <class TMeshType>
class ParameterizedMeshASCIIWriter : public MeshFileWriter<TMeshType>
{
  public:

    typedef ParameterizedMeshASCIIWriter<TMeshType> Self;
    typedef MeshFileWriter<TMeshType> Superclass;
    typedef SmartPointer < Self > Pointer;
    typedef SmartPointer < const Self > ConstPointer;

    itkTypeMacro( ParameterizedMeshASCIIWriter<TMeshType>, MeshFileWriter<TMeshType> );
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
	      for (unsigned int v=0; v < mesh->GetNumberOfPoints(); v++) {
		      mesh->GetPoint( v, &meshPoint );
		      fprintf( file, "%f %f %f\n", meshPoint[0], meshPoint[1], meshPoint[2] );
        }
        fclose( file );
      }
      // write parameterization coordinates:
      typename TMeshType::SphericalMapType *sphericalMap = mesh->GetSphericalMap();
      sprintf( fileName, "%s.par", this->m_filePrefix );
      file = fopen( fileName, "w" );
      if (file) {
        typename TMeshType::PointType paramPoint;
	      for (unsigned int v=0; v < sphericalMap->size(); v++) {
		      paramPoint = (*sphericalMap)[v];
		      fprintf( file, "%f %f %f\n", paramPoint[0], paramPoint[1], paramPoint[2] );
        }
        fclose( file );
      }
    }

   
protected:
  ParameterizedMeshASCIIWriter() {};
  virtual ~ParameterizedMeshASCIIWriter() {};

};

} // namespace itk
#endif /* MESHASCIIWRITER_H_HEADER_INCLUDED */
