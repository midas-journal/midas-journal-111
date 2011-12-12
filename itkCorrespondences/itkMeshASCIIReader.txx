#ifndef _itkMeshASCIIReader_txx
#define _itkMeshASCIIReader_txx

#include "itkMeshASCIIReader.h"

namespace itk
{

/*
 *
 */
template<class TOutputMesh>
MeshASCIIReader<TOutputMesh>
::MeshASCIIReader()
{
  /*
   * Create the output
   */
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer());
  m_FilePrefix = "";
  m_Error = false;
}

/*
 *
 */
template<class TOutputMesh>
void
MeshASCIIReader<TOutputMesh>
::GenerateData()
{
  m_Error = true;

  // open pts file:
  unsigned long numPts = -1;  // gets() alway reads one last line before eof
  std::string ptsName = m_FilePrefix + ".pts";
  FILE *ptsFile = fopen( ptsName.c_str(), "r" );
  if (ptsFile==NULL) return;
  
  // count number of points:
  char buffer[4096];
  while( feof( ptsFile ) == 0 ) {
    fgets( buffer, 4095, ptsFile );
    numPts++;
  }
      
  // memory allocation for nodes
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();  
  outputMesh->GetPoints()->Reserve( numPts );
  PointsContainerPointer  myPoints = outputMesh->GetPoints();
  typename PointsContainer::Iterator   point = myPoints->Begin();
  OPointType pnt;

  // read points:
  float pos[3];
  rewind( ptsFile );
  while( point != myPoints->End() ) {
    fscanf( ptsFile, "%f %f %f\n", &pos[0], &pos[1], &pos[2] );
    for (int d=0; d<3; d++) pnt[d] = pos[d];
    point.Value() = pnt;
    ++point;
  }
  fclose( ptsFile );

  // open fce file:
  unsigned long numTris = -1; // gets() alway reads one last line before eof
  std::string triName = m_FilePrefix + ".fce";
  FILE *triFile = fopen( triName.c_str(), "r" );
  if (triFile==NULL) return;
  
  // count number of triangles:
  while( feof( triFile ) == 0 ) {
    fgets( buffer, 4095, triFile );
    numTris++;
  }

  // the temporary container of nodes' connectness
  unsigned int tripoints[3] = {0,1,2};
  outputMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );
  CellAutoPointer testCell;

  // read triangles:
  rewind( triFile );
  for (unsigned long i=0; i<numTris; i++) {
    fscanf( triFile, "%u %u %u\n", &tripoints[0], &tripoints[1], &tripoints[2] );
    unsigned long longIds[3];
    for (unsigned int p=0; p<3; p++) { longIds[p] = (unsigned long)tripoints[p]; } 
    testCell.TakeOwnership( new TriCellType );
    testCell->SetPointIds( longIds );
    outputMesh->SetCell( i, testCell );
//    outputMesh->SetCellData(i, (OPixelType)3.0);
  }

  fclose( triFile );
  outputMesh->BuildCellLinks();
  outputMesh->SetBufferedRegion( 0 );
  m_Error = false;
}


template<class TOutputMesh>
void
MeshASCIIReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FilePrefix: " << m_FilePrefix << std::endl;
}

} // end namespace itk

#endif
