#include "itkMeshDataSpecialization.h"

#include "itkMeshASCIIReader.h"
#include "itkParameterizedMeshASCIIReader.h"
#include "itkMeshSTLWriter.h"
#include "itkParameterizedMeshASCIIWriter.h"
#include "itkIndexedTriangleMesh.h"
#include "itkConformalSphericalParameterizationFilter.h"
#include "itkRemeshParameterizedMeshFilter.h"
#include "itkStatisticalShapeModel3DCalculatorWithFeatures.h"
#include "itkSimplifiedMDLCostFunction.h"
#include "itkExceptionObject.h"

#include <itkMetaMeshConverter.h>
#include <iostream>
#include <fstream>
#include <string>

#include "itkMeshDataSpecialization.h"

typedef itk::Vector<double, 2> PixelType;
typedef unsigned short IndexType;
typedef itk::IndexedTriangleMesh<PixelType, IndexType>      LandmarkType;
typedef itk::MeshASCIIReader<LandmarkType>                  LandmarkReader;
typedef itk::IndexedTriangleMesh<PixelType, IndexType>      InputMeshType;
typedef itk::MeshASCIIReader<InputMeshType>                 InputMeshReader;
typedef itk::SphericalParameterizedTriangleMesh
<PixelType, IndexType>                         ParameterizationType;
typedef itk::ParameterizedMeshASCIIReader
<ParameterizationType>                         ParameterizationReader;
typedef itk::ParameterizedMeshASCIIWriter
<ParameterizationType>                         ParameterizationWriter;
typedef itk::ConformalSphericalParameterizationFilter
<InputMeshType, ParameterizationType>          ParameterizationFilterType;
typedef LandmarkType                                        OutputMeshType;
typedef itk::MeshSTLWriter<OutputMeshType>                  OutputMeshWriter;
typedef itk::RemeshParameterizedMeshFilter
<LandmarkType, ParameterizationType, OutputMeshType> LandmarkFilterType;
typedef std::vector< std::string >                StringVector;


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


StringVector fileNames;
StringVector fileExtensions;
double modelSize = 0;
unsigned int numSamples = 0;
std::string landmarkName = "";
std::string outputLandmarkName = "";

int main( int argc, char *argv[] )
{
  // parse command line arguments
  if (argc < 4)
  {
    std::cout << "Usage: " << argv[0] << " MeshListFile LandmarkFile ModelRadius [OutputLandmarkFile]" << std::endl;
    return -1;
  }
  std::string buffer, extension, outputExtension;
  std::string shapeListName = argv[1];
  landmarkName = argv[2];
  std::ifstream shapeListFile;
  shapeListFile.open( shapeListName.c_str(), std::ios_base::in );
  if (!shapeListFile)
  {
    std::cerr << "Unable to open shape list \"" << shapeListName << "\"!" << std::endl;
    return 1;
  }

  int lastPoint = landmarkName.rfind( '.' );
  extension.erase();
  if (lastPoint>0) { extension = landmarkName.substr( lastPoint ); }
  if ( extension.compare( ".pts" )!=0 )
  {
    std::cerr << "Landmark file \"" << landmarkName << "\" has no valid extension!" << std::endl;
    return 2;
  }
  landmarkName = landmarkName.erase( lastPoint, 4 );
  
  modelSize = atof( argv[3] );
  if (modelSize <= 0)
  {
    std::cerr << "Invalid ModelRadius!" << std::endl;
    return 3;
  }
  
  if (argc > 4) 
  {
    outputLandmarkName = argv[4];
    int outputLastPoint = outputLandmarkName.rfind( '.' );
    outputExtension.erase();
    if (outputLastPoint>0) { outputExtension = outputLandmarkName.substr( outputLastPoint ); }
    if ( outputExtension.compare( ".pts" )!=0 )
    {
      std::cerr << "Landmark file \"" << outputLandmarkName << "\" has no valid extension!" << std::endl;
      return 4;
    }
    outputLandmarkName = outputLandmarkName.erase( outputLastPoint, 4 );
  }

  // parse the shape list file
  numSamples = 0;
  while (!shapeListFile.eof())
  {
    std::getline( shapeListFile, buffer );
    #ifndef WIN32
      if (buffer.empty()) continue;
      // files created on Windows and read on Linux show additional '\r's -> delete them evil bastards!
      int lastCharPos = buffer.length()-1;
      if (buffer[lastCharPos]=='\r') buffer.erase( lastCharPos );
      // delete leading g: if detected
      if (buffer.compare( 0, 2, "g:")==0) buffer.erase( 0, 2 );
    #endif
    if (buffer.empty() || buffer[0]=='#') { continue; }
    fileNames.resize( numSamples+1 );
    fileExtensions.resize( numSamples+1 );
    // check file extension
    int lastPoint = buffer.rfind( '.' );
    extension.erase();
    if (lastPoint>0) { extension = buffer.substr( lastPoint ); }
    if ( extension.compare( ".pts" )!=0 )
    {
      std::cerr << "Shape file \"" << buffer << "\" has no valid extension!" << std::endl;
      return 5;
    }
    // delete extension (for reader)
    fileNames[numSamples] = buffer.erase( lastPoint, 4 );
    fileExtensions[numSamples] = extension;
    numSamples++;
  }
  shapeListFile.close();
  

  // determine whether the data includes local features
  bool hasFeatures = true ;
  ifstream tempFile ;
  for ( unsigned int i = 0 ; hasFeatures && ( i < numSamples ) ; i++ )
  {
    std::string currentFileName ;
    currentFileName = fileNames[i] + ".txt" ;
    tempFile.open ( currentFileName.c_str () ) ;
    if ( ! tempFile.is_open () )
    {
      hasFeatures = false ;
    }
  }
  if ( hasFeatures ) 
  {
    std::cout << "Has Features" << std::endl ;
  }
  else
  {
    std::cout << "Has NO Features" << std::endl ;
  }
  
  

  // declare both a standard and a feature-based correspondence optimizer
  // the appropriate optimizer will be used in the rest of the program
  typedef itk::StatisticalShapeModel3DCalculator<LandmarkType, ParameterizationType, OutputMeshType> SSMCalculatorType;
  typedef itk::SimplifiedMDLCostFunction<SSMCalculatorType>     CostFunctionType;
  
  typedef itk::StatisticalShapeModel3DCalculatorWithFeatures<LandmarkType, ParameterizationType, OutputMeshType> SSMCalculatorTypeWF;
  typedef itk::SimplifiedMDLCostFunction<SSMCalculatorTypeWF>     CostFunctionTypeWF;
  
  CostFunctionType::Pointer mdlCosts ;
  CostFunctionTypeWF::Pointer mdlCostsWF ;
  SSMCalculatorType::Pointer ssmCalculator ;
  SSMCalculatorTypeWF::Pointer ssmCalculatorWF ; 

  // set up the pipeline for the appropriate correspondence optimizer
  if ( hasFeatures ) 
  {
    mdlCostsWF = CostFunctionTypeWF::New();
    mdlCostsWF->SetVarianceCutForModelRadius( modelSize );
    ssmCalculatorWF = SSMCalculatorTypeWF::New();
    ssmCalculatorWF->SetCostFunction( mdlCostsWF );
    ssmCalculatorWF->SetNumberOfInputs( numSamples );
  }
  else
  {
    mdlCosts = CostFunctionType::New() ;
    mdlCosts->SetVarianceCutForModelRadius( modelSize );
    ssmCalculator = SSMCalculatorType::New();
    ssmCalculator->SetCostFunction( mdlCosts );
    ssmCalculator->SetNumberOfInputs( numSamples );
  }  

  // read in the landmarks
  LandmarkReader::Pointer landmarkReader = LandmarkReader::New();
  landmarkReader->SetFilePrefix( landmarkName.c_str() );
  std::cout << "Reading \"" << landmarkName << "\"" << std::endl;
  landmarkReader->Update();
  if (landmarkReader->GetReadError())
  {
    std::cerr << "Unable to open \"" << landmarkName << "\"!" << std::endl;
    return 6;
  }
  LandmarkType::Pointer landmarks = landmarkReader->GetOutput();

  // set the landmarks for the appropriate correspondence optimizer
  if ( hasFeatures ) 
  {
    ssmCalculatorWF->SetLandmarks( landmarks );
  }
  else
  {
    ssmCalculator->SetLandmarks( landmarks );
  }
  
  // read in the shapes
  for (unsigned int i=0; i<numSamples; i++) 
  {
    ParameterizationReader::Pointer parameterizationReader = ParameterizationReader::New();
    parameterizationReader->SetFilePrefix( fileNames[i].c_str() );
    std::cout << "Reading \"" << fileNames[i] << "\"" << std::endl;
    parameterizationReader->Update();
    ParameterizationType::Pointer parameterization = parameterizationReader->GetOutput();
    
    if (parameterizationReader->GetReadError())
    {
      // <ipek>
      std::cerr << "Attempted generating conformal parametrization." << std::endl 
        << "Check parametrization for: " << fileNames[i] << std::endl ;
      // </ipek> 

      InputMeshReader::Pointer sampleReader = InputMeshReader::New();
      sampleReader->SetFilePrefix( fileNames[i].c_str() );
      sampleReader->Update();
      if (sampleReader->GetReadError())
      {
        std::cerr << "Unable to open \"" << fileNames[i] << "\"!" << std::endl;
        return 7;
      }
      InputMeshType::Pointer inputMesh = sampleReader->GetOutput();
      std::cout << "Generating parameterization for \"" << fileNames[i] << "\"." << std::endl;
      ParameterizationFilterType::Pointer parameterizationFilter = ParameterizationFilterType::New();
      parameterizationFilter->SetInput( inputMesh );
      parameterizationFilter->Update();
      parameterization = parameterizationFilter->GetOutput();
      // save parameterization
      std::cout << "Writing \"" << fileNames[i] << ".par\"" << std::endl;
      ParameterizationWriter::Pointer writer = ParameterizationWriter::New();
      writer->SetFilePrefix( fileNames[i].c_str() );
      writer->SetInput( parameterization );
      writer->Write();
    }
    // connect the current shape to the appropriate correspondence optimizer
    if ( hasFeatures ) 
    {
      ssmCalculatorWF->SetInput( i, parameterization );
    }
    else
    {
      ssmCalculator->SetInput( i, parameterization );
    }
  }

  // update the appropriate correspondence optimizer
  std::cout << std::endl << "Starting correspondence optimization with " << numSamples << " shapes." << std::endl; 
  try
  {
    if ( hasFeatures ) 
    {
      ssmCalculatorWF->Update();
    }
    else
    {
      ssmCalculator->Update();
    }
  }
  catch( itk::ExceptionObject& excp ) {
    excp.Print( std::cout );
    char x;
    std::cin >> x;
    return 8 ;
  }
  
  if (outputLandmarkName != "" )
  {
    // read the landmark configuration to be used for the output meshes
    LandmarkReader::Pointer outputLandmarkReader = LandmarkReader::New();
    outputLandmarkReader->SetFilePrefix( outputLandmarkName.c_str() );
    std::cout << "Reading \"" << outputLandmarkName << "\"" << std::endl;
    outputLandmarkReader->Update();
    if (outputLandmarkReader->GetReadError())
    {
      std::cerr << "Unable to open \"" << outputLandmarkName << "\"!" << std::endl;
      return 9 ;
    }
    LandmarkType::Pointer outputLandmarks = outputLandmarkReader->GetOutput();
    if ( hasFeatures ) 
    {
      ssmCalculatorWF->SetOutputLandmarks( outputLandmarks );
    }
    else
    {
      ssmCalculator->SetOutputLandmarks( outputLandmarks );
    }
  }
  // write out the resulting shapes
  typedef itk::MeshSpatialObject < InputMeshType::Superclass > itkMeshSOType ;
  typedef itk::MetaMeshConverter < 3, PixelType, InputMeshType::MeshTraits > MeshConverterType ;
  MeshConverterType * itkConverter = new MeshConverterType() ;
  itkMeshSOType::Pointer meshSO = itkMeshSOType::New() ;
    
  char outFileName[256];
  for (unsigned int i=0; i<numSamples; i++) 
  {  
    sprintf( outFileName, "result%03i.meta", i );
    std::cout << "Writing " << outFileName << std::endl; 
    OutputMeshType::Pointer resultMesh ;
    if ( hasFeatures ) 
    {
      resultMesh = ssmCalculatorWF->GetResampledOutputMesh( i ) ;
    }
    else
    {
      resultMesh = ssmCalculator->GetResampledOutputMesh( i ) ;
    }
    meshSO->SetMesh((InputMeshType::Superclass::Pointer) resultMesh ) ;
    itkConverter->WriteMeta (meshSO, outFileName) ;
  }

  std::cout << std::endl << "Finished." << std::endl; 

  return 0 ;
}

