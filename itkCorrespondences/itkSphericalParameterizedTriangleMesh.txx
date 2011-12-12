#ifndef _itkSphericalParameterizedTriangleMesh_txx
#define _itkSphericalParameterizedTriangleMesh_txx

#include "itkSphericalParameterizedTriangleMesh.h"
#include <vnl/vnl_cross.h>
#include <map>


// macros for fast triangle intersection tests:
#define EPSILON 0.0000000001
#define CROSS(dest,v1,v2) \
  dest[0]=(v1[1])*v2[2]-(v1[2])*v2[1]; \
  dest[1]=(v1[2])*v2[0]-(v1[0])*v2[2]; \
  dest[2]=(v1[0])*v2[1]-(v1[1])*v2[0];
#define DOT(v1,v2) \
  (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) \
  dest[0]=v1[0]-v2[0]; \
  dest[1]=v1[1]-v2[1]; \
  dest[2]=v1[2]-v2[2]; 


namespace itk
{
  
  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  bool
  SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::MapCoordinates( const PointType parameterCoordinates, IndexType &faceIndex, double &barycentricP, double &barycentricQ ) const
  {
    // sort faces by distance to parameterCoordinates:
    typedef std::map<double, IndexType> MapType;
    typedef std::pair <double, IndexType> MapPair;
    MapType faceDistances;
    double distThresh = 0.5;
    CoordRepType faceXDiff, faceYDiff, faceZDiff;
    int v0, v1, v2;
    for (IndexType f=0; f<this->GetNumberOfFaces(); f++) {
      v0 = this->GetPointIndexForFace( f, 0 );
      v1 = this->GetPointIndexForFace( f, 1 );
      v2 = this->GetPointIndexForFace( f, 2 );
      faceXDiff = ( m_SphericalMap[v0][0] + m_SphericalMap[v1][0] + m_SphericalMap[v2][0] ) / 3.0 - parameterCoordinates[0];
      faceYDiff = ( m_SphericalMap[v0][1] + m_SphericalMap[v1][1] + m_SphericalMap[v2][1] ) / 3.0 - parameterCoordinates[1];
      faceZDiff = ( m_SphericalMap[v0][2] + m_SphericalMap[v1][2] + m_SphericalMap[v2][2] ) / 3.0 - parameterCoordinates[2];
      double dist = faceXDiff*faceXDiff + faceYDiff*faceYDiff + faceZDiff*faceZDiff;
      if (dist < distThresh) 
      {
        faceDistances.insert( MapPair( dist, f ) );
      }
    }

    // find intersecting triangle and barycentric coordinates:
    for (typename MapType::iterator it=faceDistances.begin(); it!=faceDistances.end(); it++) {
      //double dist = it->first;
      faceIndex = it->second;
      if (this->CoordinatesInFace( faceIndex, parameterCoordinates, barycentricP, barycentricQ )) {
        return true;
      }
    }

    // last resort: brute force search across all faces (probably distThresh was too low)
    for (faceIndex=0; faceIndex<this->GetNumberOfFaces(); faceIndex++) {
      if (this->CoordinatesInFace( faceIndex, parameterCoordinates, barycentricP, barycentricQ )) {
        return true;
      }
    }

    // this is very bad :-(
    return false;
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  bool
  SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::CoordinatesInFace( int faceId, const PointType parameterCoordinates, double &p1, double &p2 ) const
  {
    /* This code is a customized version of the algorithm by Tomas Moller and
    * Ben Trumbore, presented in the Journal of Graphics Tools (JGT).
    * http://www.acm.org/jgt/papers/MollerTrumbore97/
    */
    double edge1[3], edge2[3], pvec[3], qvec[3];
    double det, inv_det;

    int v0 = this->GetPointIndexForFace( faceId, 0 );
    int v1 = this->GetPointIndexForFace( faceId, 1 );
    int v2 = this->GetPointIndexForFace( faceId, 2 );

    // find vectors for two edges sharing vert0 
    SUB(edge1, m_SphericalMap[v1], m_SphericalMap[v0]);
    SUB(edge2, m_SphericalMap[v2], m_SphericalMap[v0]);

    // begin calculating determinant - also used to calculate U parameter 
    CROSS(pvec, parameterCoordinates, edge2);
    // if determinant is near zero, ray lies in plane of triangle
    det = DOT(edge1, pvec);
    if (det > EPSILON)
    {
      // calculate U parameter and test bounds
      p1 = -DOT(m_SphericalMap[v0], pvec);
      if (p1 < -EPSILON || p1 > det+EPSILON) return false;
      // prepare to test V parameter
      CROSS(qvec, -(m_SphericalMap[v0]), edge1);
      // calculate V parameter and test bounds
      p2 = DOT(parameterCoordinates, qvec);
      if (p2 < -EPSILON || p1+p2 > det+EPSILON) return false;
    }
    else if(det < -EPSILON)
    {
      // calculate U parameter and test bounds
      p1 = -DOT(m_SphericalMap[v0], pvec);
      if (p1 > EPSILON || p1 < det-EPSILON) return false;
      // prepare to test V parameter
      CROSS(qvec, -m_SphericalMap[v0], edge1);
      // calculate V parameter and test bounds
      p2 = DOT(parameterCoordinates, qvec) ;
      if (p2 > EPSILON || p1+p2 < det-EPSILON) return false;
    }
    else return false;  // ray is parallel to the plane of the triangle

    inv_det = 1.0 / det;
    if (DOT(edge2, qvec) * inv_det <= 0) return false;  // triangle is on the other side of ray

    // calculate u+v, ray intersects triangle
    p1 *= inv_det;
    p2 *= inv_det;
    return true;
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  void
  SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::UpdateParameterization( ParameterizedTriangleMeshPointer source )
  {
    source->Update();
    bool modified = false;
    SphericalMapType *newMap = (dynamic_cast<Self*>(source.GetPointer()))->GetSphericalMap();
    for (unsigned int i=0; i<m_SphericalMap.size(); i++)
    {
      if (m_SphericalMap[i] != (*newMap)[i])
      {
        m_SphericalMap[i] = (*newMap)[i];
        this->SetParameterizationModified( i, true );
        modified = true;
      }
    }
    if (modified) { this->Modified(); }
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  void
  SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::InitializeSphericalMap()
  {
    m_SphericalMap.resize( this->GetNumberOfPoints() );
    this->m_PointModified.resize( this->GetNumberOfPoints(), true );
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  typename SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>::IndexType
  SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::GetPatchIndexForFace( IndexType faceId ) const
  {
    IndexType vertexId[3];
    for (IndexType i=0; i<3; i++) 
    {
      vertexId[i] = this->GetPointIndexForFace( faceId, i );
    }
    PointType p0 = m_SphericalMap[vertexId[0]];
    PointType p1 = m_SphericalMap[vertexId[1]];
    PointType p2 = m_SphericalMap[vertexId[2]];
    IndexType patchId;
    CoordRepType meanX = p0[0] + p1[0] + p2[0];
    if (meanX >= 0) patchId = 0;
    else            patchId = 1;
    return patchId;
  }

  
  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  typename SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>::PatchPointType 
  SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::MapParameterizationToPatch( const PointType parameterCoordinates, IndexType patchIdx ) const
  {
    const double halfPi = 1.5707963267948966;
    PatchPointType patchCoords;
    CoordRepType azimuth, polar, k;
    if (patchIdx == 0) 
    {
      azimuth = atan2( parameterCoordinates[1], parameterCoordinates[0] );
    }
    else
    {
      azimuth = atan2( parameterCoordinates[1], -parameterCoordinates[0] );
    }
    // make sure we don't have floating point accuracy problems here:
    if (parameterCoordinates[2] < -1.0) polar = acos( -1.0 ) - halfPi;
    else if (parameterCoordinates[2] > 1.0) polar = acos( 1.0 ) - halfPi;
    else polar = acos( parameterCoordinates[2] ) - halfPi;
    k = 0.5 / ( 1.0 + cos( polar )*cos( azimuth ) ) ;
    patchCoords[0] = k * cos( polar ) * sin( azimuth );
    patchCoords[1] = k * sin( polar );
    return patchCoords;
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  typename SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>::PointType 
  SphericalParameterizedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::MapPatchToParameterization( const PatchPointType patchCoordinates, IndexType patchIdx ) const
  {
    const double halfPi = 1.5707963267948966;
    CoordRepType azimuth, polar;
    CoordRepType x = 4.0*patchCoordinates[0];
    CoordRepType y = 4.0*patchCoordinates[1];
    CoordRepType p = sqrt( x*x + y*y );
    if (p==0) {
      azimuth = 0;
      polar = halfPi;
    }
    else {
      CoordRepType c = 2.0 * atan( 0.5*p );
      azimuth = atan2( x*sin(c), (p*cos(c)) );
      polar = asin( y*sin(c) / p ) + halfPi;  
    }
    PointType cartCoords;
    cartCoords[0] = cos( azimuth ) * sin( polar );
    cartCoords[1] = sin( azimuth ) * sin( polar );
    cartCoords[2] = cos( polar );
    if (patchIdx == 1) { cartCoords[0] *= -1.0; }
    return cartCoords;
  }

}

#endif
