#ifndef _itkParameterizedMeshASCIIReader_txx
#define _itkParameterizedMeshASCIIReader_txx

#include "itkParameterizedMeshASCIIReader.h"
#include <fstream>
using namespace std ;

namespace itk
{

/*
 *
 */
template<class TOutputMesh>
ParameterizedMeshASCIIReader<TOutputMesh>
::ParameterizedMeshASCIIReader()
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
ParameterizedMeshASCIIReader<TOutputMesh>
::GenerateData()
{
  m_Error = true;

  // open pts file:
  unsigned int numPts = -1;  // gets() alway reads one last line before eof
  std::string ptsName = m_FilePrefix + ".pts";
  FILE *ptsFile = fopen( ptsName.c_str(), "r" );
  if (ptsFile==NULL) return;
  
  // count number of points:
  char buffer[4096];
  while( feof ( ptsFile )  == 0 ) {
     fgets(buffer, 4095, ptsFile) ;
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
  unsigned int numTris = -1; // gets() alway reads one last line before eof
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
  for (unsigned int i=0; i<numTris; i++) {
    fscanf( triFile, "%u %u %u\n", &tripoints[0], &tripoints[1], &tripoints[2] );
    unsigned long longIds[3];
    for (unsigned int p=0; p<3; p++) { longIds[p] = (unsigned long)tripoints[p]; } 
    testCell.TakeOwnership( new TriCellType );
    testCell->SetPointIds( longIds );
    outputMesh->SetCell( i, testCell );
//    outputMesh->SetCellData(i, (OPixelType)3.0);
  }
  fclose( triFile );

  // open par file:
  std::string parName = m_FilePrefix + ".par";
  FILE *parFile = fopen( parName.c_str(), "r" );
  if (parFile==NULL) return;
  
  // memory allocation 
  outputMesh->InitializeSphericalMap();
  typename OutputMeshType::SphericalMapType *map = outputMesh->GetSphericalMap();
  
  // read points:
  rewind( parFile );
  for (unsigned int i=0; i<numPts; i++) {
    fscanf( parFile, "%f %f %f\n", &pos[0], &pos[1], &pos[2] );
    for (int d=0; d<3; d++) (*map)[i][d] = pos[d];
  }
  fclose( parFile );

  outputMesh->BuildCellLinks();


  // <ipek>
  // read point data:
  std::string dataName = m_FilePrefix + ".txt";
  ifstream dataFile ;
  dataFile.open (dataName.c_str());
  bool features = dataFile.is_open () ;

  if ( ! features )
  {
    outputMesh->SetBufferedRegion( 0 );
    m_Error = false;
    return ;
  }
  // This is a quite sloppy implementation of reading keyword/value pairs!
  bool found;
  char * valuePtr;
  char typeString[1001], line[1001];
  unsigned int nPtsFile, nDim ;

  dataFile.seekg(0,std::ios::beg);
  found = false ;
  while ( !found && !dataFile.eof())
  { dataFile.getline ( line, 1000 ) ;
    if (line[0] != '#' && strstr ( line, "NUMBER_OF_POINTS" )) found = true;
  }
  valuePtr=strchr(line, '=');
  if (!valuePtr) return;
  valuePtr++;
  sscanf(valuePtr, " %d ", &nPtsFile);

  dataFile.seekg(0,std::ios::beg);
  found = false ;
  while ( !found && !dataFile.eof())
  { dataFile.getline ( line, 1000 ) ;
    if (line[0] != '#' && strstr ( line, "DIMENSION" )) found = true;
  }
  valuePtr=strchr(line, '=');
  if (!valuePtr) return;
  valuePtr++;
  sscanf(valuePtr, " %d ", &nDim);

  dataFile.seekg(0,std::ios::beg);
  found = false ;
  while ( !found && !dataFile.eof())
  { dataFile.getline ( line, 1000 ) ;
    if (line[0] != '#' && strstr ( line, "FEATURES" )) found = true;
  }
  valuePtr=strchr(line, '=');
  if (!valuePtr) return;
  valuePtr++;
  sscanf(valuePtr, " %s ", typeString);

  assert ( nPtsFile == numPts ) ;
  assert ( nDim > 0 ) ;

  dataFile.seekg(0,std::ios::beg);
  const int numEntries = 3;
  int counter = 0;
  while ( counter < numEntries && !dataFile.eof())
  { dataFile.getline ( line, 1000 ) ;
    if ((line[0] != '#')) counter++;
  }

  double *data;
  data = new double [nDim] ;

  for (unsigned int i = 0 ; i < numPts ; i++ )
  {
    for (unsigned int j = 0 ; j < nDim ; j++ )
    {
      dataFile >> data[j] ;
    }
    outputMesh->SetPointData(i, data);
  }

  dataFile.close();
  
  // </ipek>

  outputMesh->SetBufferedRegion( 0 );
  m_Error = false;
}


template<class TOutputMesh>
void
ParameterizedMeshASCIIReader<TOutputMesh>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "FilePrefix: " << m_FilePrefix << std::endl;
}

} // end namespace itk

#endif
