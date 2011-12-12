#ifndef _itkSphericalParameterizedTriangleMesh_h
#define _itkSphericalParameterizedTriangleMesh_h

#include "itkParameterizedTriangleMesh.h"

namespace itk
{

  /** \class SphericalParameterizedTriangleMesh
  *  \brief Implements a surface parameterization of a spherical (genus zero) 
  *         2-manifold mesh consisting exclusively of triangle cells.
  *
  * Internally, the parameterization is represented as a second, spherical 
  * mesh with the same topology as the original mesh. Coordinates of the 
  * spherical mesh are stored in a STL-array. Parameter coordinates are 
  * represented as 3D coordinates lying on the unit sphere. 
  * The mapping from parameter to object space is implemented by calculating 
  * the intersection of a ray (from the center of the sphere to a point 
  * indicated by the parameter coordinates) and the triangles of the sphere. 
  * If the intersecting triangle cannot be guessed, MapCoordinates() sorts
  * all triangles according to their average distance to the parameter 
  * coordinates and starts intersection testing with the closest. 
  * To map the parameterization to 2D-images, two patches are used (each 
  * mapping one half-sphere by means of the stereographic projection).
  * New methods of the class are GetSphericalMap(), which returns a reference 
  * to the coordinate array of the spherical mesh and InitializeSphericalMap(), 
  * which reserves the necessary memory for the coordinate array and is 
  * usually called once by filters.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
  *         German Cancer Research Center, Heidelberg, Germany.
  */
  template<
    typename TPixelType, 
    typename TIndex = unsigned int,
    typename TCoordRepType = double
  >
  class SphericalParameterizedTriangleMesh : public ParameterizedTriangleMesh
    <TPixelType, TIndex, TCoordRepType>
  {

  public:

    /** Standard class typedefs. */
    typedef SphericalParameterizedTriangleMesh      Self;
    typedef ParameterizedTriangleMesh
      <TPixelType, TIndex, TCoordRepType>           Superclass;
    typedef SmartPointer<Self>                      Pointer;
    typedef SmartPointer<const Self>                ConstPointer;

    /** Convenient typedefs. */
    typedef Superclass                              ParameterizedTriangleMeshType;
    typedef typename Superclass::Pointer            ParameterizedTriangleMeshPointer;
    typedef typename Superclass
                ::IndexedTriangleMeshConstPointer   IndexedTriangleMeshConstPointer;
    typedef typename Superclass
      ::IndexedTriangleMeshPointer                  IndexedTriangleMeshPointer;
    typedef typename Superclass::PointType          PointType;
    typedef typename Superclass::VectorType         VectorType;
    typedef typename Superclass::CoordRepType       CoordRepType;
    typedef typename Superclass::IndexType          IndexType;
    typedef typename Superclass::PatchPointType     PatchPointType;
    typedef std::vector<PointType>                  SphericalMapType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SphericalParameterizedTriangleMesh, ParameterizedTriangleMesh);

    /** Finds the the barycentric coordinates and corresponding triangle for 
    * coordinates in parameter space. 
    * \returns true if parameter coordinates were valid and mapping 
    *          successful, else false
    */
    virtual bool MapCoordinates( const PointType parameterCoordinates, 
      IndexType &faceIndex, double &barycentricP, double &barycentricQ ) const;

    /** Checks if the supplied parameterCoordinates are lying within the 
    * specified face.If yes, the method returns true and writes the exact 
    * barycentric coordinates to p1 and p2.
    */
    virtual bool CoordinatesInFace( int faceId, 
      const PointType parameterCoordinates, double &p1, double &p2 ) const;

    /** Copies the parameterization of source. */
    virtual void UpdateParameterization( ParameterizedTriangleMeshPointer source );
    
    /** This class uses sterographic projection of two half-spheres to 
    * calculate the patch coordinates, i.e. there are two patches.
    */
    virtual IndexType GetNumberOfPatches() const
    {
      return 2;
    }

    /** Returns the index of the patch which represents the specified
    * triangle with the lowest possible distortion.
    */
    virtual IndexType GetPatchIndexForFace( IndexType faceId ) const;

    /** Maps the specified coordinates from parameter space to the given patch.
    * The returned coordinates lie between [0..1] for both dimensions. 
    */
    virtual PatchPointType MapParameterizationToPatch( const PointType parameterCoordinates, IndexType patchIdx ) const;

    /** Maps the specified patch coordinates from the given patch to
    * parameter space.
    * Patch coordinates have to lie between [0..1] for both dimensions. 
    */
    virtual PointType MapPatchToParameterization( const PatchPointType patchCoordinates, IndexType patchIdx ) const;

    /** Copies all pointers of itk::Mesh (shallow copy) and all indexing 
    * information of itk::IndexedTriangleMesh (deep copy). 
    * Does not copy the parameterization data (the spherical map)!
    */
    Pointer operator=( IndexedTriangleMeshPointer mesh )
    {
      Superclass::operator=( mesh );
      return this;
    }

    /** Returns the spherical map for the parameterized mesh. */
    SphericalMapType* GetSphericalMap() 
    { 
      return &m_SphericalMap; 
    }

    /** Reserves the necessary memory to hold all map points. */
    void InitializeSphericalMap();

   
  private:

    /** Holds the spherical map. */
    SphericalMapType    m_SphericalMap;
    

  };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSphericalParameterizedTriangleMesh.txx"
#endif

#endif
