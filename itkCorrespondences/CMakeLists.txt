INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

#Change PROJECT_NAME to the name of your project
PROJECT(CreateShapeModel)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
SET(Required_Packages
ITK
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
ITKCommon
ITKIO
ITKNumerics
ITKBasicFilters 
ITKStatistics 
)

#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project

SET( SOURCES ${SOURCES}
itkIndexedTriangleMesh.h itkIndexedTriangleMesh.txx
itkParameterizedTriangleMesh.h
itkSphericalParameterizedTriangleMesh.h itkSphericalParameterizedTriangleMesh.txx
itkConformalSphericalParameterizationFilter.h itkConformalSphericalParameterizationFilter.txx
itkRemeshParameterizedMeshFilter.h itkRemeshParameterizedMeshFilter.txx
itkStatisticalShapeModel3DCalculator.h itkStatisticalShapeModel3DCalculator.txx
itkStatisticalShapeModel3DCalculatorWithFeatures.h itkStatisticalShapeModel3DCalculatorWithFeatures.txx
itkProcrustesAlign3DMeshFilter.h itkProcrustesAlign3DMeshFilter.txx
itkGaussianWarpSphericalParameterizationFilter.h itkGaussianWarpSphericalParameterizationFilter.txx
itkRotateSphericalParameterizationFilter.h itkRotateSphericalParameterizationFilter.txx
itkShapeModelCalculatorCostFunction.h itkShapeModelCalculatorCostFunction.txx
itkVarianceBasedCostFunction.h itkVarianceBasedCostFunction.txx
itkSimplifiedMDLCostFunction.h itkSimplifiedMDLCostFunction.txx
itkMeshASCIIReader.h itkMeshASCIIReader.txx
itkParameterizedMeshASCIIReader.h itkParameterizedMeshASCIIReader.txx
itkMeshFileWriter.h itkMeshSTLWriter.h itkParameterizedMeshASCIIWriter.h
itkMeshDataSpecialization.h
)

SET(CurrentExe "CreateShapeModel")
ADD_EXECUTABLE(${CurrentExe} ${SOURCES} main.cpp) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
#ADD_TEST(CreateImage ImageCopy ${CMAKE_SOURCE_DIR}/img1.png)
#ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/img1.png ${CMAKE_SOURCE_DIR}/img1.png.base.png)

#Once you're done modifying this template, you should rename it to "CMakeLists.txt"
