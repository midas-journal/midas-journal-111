#ifndef _itkIndexedTriangleMesh_h
#define _itkIndexedTriangleMesh_h

#include "vnl/vnl_vector_fixed.h"
#include "itkMesh.h"
#include <vector>


namespace itk
{
  /** \class IndexedTriangleMesh
  *   \brief The class represents a 3D, 2-manifold mesh consisting exclusively 
  *          of triangle cells and allows advanced queries. 
  *
  * The IndexedTriangleMesh stores adjacency information for points, edges and 
  * faces. The information is calculated in BuildCellLinks(). 
  * The template parameter TIndex specifies the type of the link indices. 
  * To save memory, this can be set to unsigned short for meshes with less than 
  * 65.536 faces, edges and vertices.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
  *         German Cancer Research Center, Heidelberg, Germany.
  */
  template<
    typename TPixelType, 
    typename TIndex = unsigned int,
    typename TCoordRepType = double
  >
  class IndexedTriangleMesh : public Mesh
    <
    TPixelType, 
    3, 
    DefaultStaticMeshTraits<TPixelType, 3, 2, TCoordRepType, 
                            TCoordRepType, TPixelType>
    >
  {

  public:

    /** Standard typedefs. */
    typedef IndexedTriangleMesh                           Self;
    typedef Mesh<TPixelType, 3, DefaultStaticMeshTraits
      <TPixelType, 3, 2, TCoordRepType, 
       TCoordRepType, TPixelType> >                       Superclass;
    typedef SmartPointer<Self>                            Pointer;
    typedef SmartPointer<const Self>                      ConstPointer;

    /** Convenient typedefs. */
    typedef TIndex                                        IndexType;
    typedef TCoordRepType                                 CoordRepType;
    typedef Superclass                                    MeshType;
    typedef Vector<CoordRepType, 3>                       VectorType;
    typedef typename MeshType::PointType                  PointType;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Standard part of every itk Object. */
    itkTypeMacro(IndexedTriangleMesh, Mesh);

    /** Returns a reference to a point of the mesh. */
    PointType& GetPoint( IndexType pointId ) 
    { 
      return (*(this->m_PointsContainer))[pointId]; 
    }

    /** Inherited from itk::Mesh class. */
    bool GetPoint( typename Superclass::PointIdentifier pointId, PointType *point ) const
    {
      return Superclass::GetPoint( pointId, point );
    }

    /** Returns the number of edges in the mesh. 
     * Adjacent faces share an edge.
     */
    unsigned int GetNumberOfEdges() const 
    { 
      return m_NumberOfEdges; 
    }

    /** Returns a vector along the specified edge from 
     *  local point 0 to local point 1. 
     */
    VectorType GetEdge( IndexType edgeId ) const 
    { 
      return const_cast<Self*>(this)->GetPoint( GetPointIndexForEdge( edgeId, 1 ) ) - 
             const_cast<Self*>(this)->GetPoint( GetPointIndexForEdge( edgeId, 0 ) ); 
    }

    /** Returns the length of the specified edge . */
    CoordRepType GetEdgeLength( IndexType edgeId ) const 
    { 
      return GetEdge( edgeId ).GetNorm(); 
    }
    
    /** Returns the number of edges connected to the specified point. */
    unsigned int GetNumberOfEdgesForPoint( IndexType pointId ) const 
    { 
      return m_EdgesPerVertex[pointId].size(); 
    }

    /** Returns the index of an edge connected to the specified point.
    * A value of 0 for localEdgeId will return the first connected edge,
    * increasing values the other ones.
    */
    IndexType GetEdgeIndexForPoint( IndexType pointId, IndexType localEdgeId ) 
      const 
    { 
      return m_EdgesPerVertex[pointId][localEdgeId]; 
    }

    /** Returns the point indices that form the specified edge.
    * Valid values for localVertexId are 0 and 1.
    */
    IndexType GetPointIndexForEdge( IndexType edgeId, IndexType localPointId )
      const 
    { 
      return m_EdgePoints[edgeId].edgeVertex[localPointId]; 
    }

    /** Returns the number of faces/triangles in the mesh. */
    unsigned int GetNumberOfFaces() const 
    { 
      return m_NumberOfFaces; 
    }

    /** Returns the number of faces adjacent to the specified point. */
    unsigned int GetNumberOfFacesForPoint( IndexType pointId ) const
    { 
      return m_FaceLinks[pointId].size(); 
    }

    /** Returns the index of a face adjacent to the specified point. 
    * A value of 0 for localFaceId will return the first adjacent face,
    * increasing values the other ones. 
    */
    IndexType GetFaceIndexForPoint( IndexType pointId, IndexType localFaceId ) 
      const
    {
      return m_FaceLinks[pointId][localFaceId];
    }
    
    /** Returns the point indices that form the specified face.
    * Since all faces are triangles, valid values for localPointId are 0..2.
    */
    IndexType GetPointIndexForFace( IndexType faceId, IndexType localPointId ) 
      const 
    { 
      return m_FaceVertices[faceId][localPointId]; 
    }

    /** Returns the edge indices that form the specified face.
    * Since all faces are triangles, valid values for localEdgeId are 0..2.
    */
    IndexType GetEdgeIndexForFace( IndexType faceId, IndexType localEdgeId ) 
      const;

    /** Returns the number of faces adjacent to the specified face.
    * Adjacent means that the two faces share a common vertex,
    * but not necessarily a common edge.
    */
    IndexType GetNumberOfAdjacentFaces( IndexType faceId ) const
    {
      return m_AdjacentFaces[faceId].size();
    }

    /** Returns the index of a face adjacent to the specified face. 
    * A value of 0 for localFaceId will return the first adjacent face,
    * increasing values the other ones. 
    */
    IndexType GetAdjacentFaceIndex( IndexType faceId, IndexType localFaceId ) 
      const
    {
      return m_AdjacentFaces[faceId][localFaceId];
    }

    /** Returns the index of the edge connecting the two specified points.
    * If there is no edge connecting the points, the method returns -1.
    */
    long int GetConnectingEdgeIndex( IndexType pointId0, IndexType pointId1 ) 
      const;
    
    /** Returns the index of a point that is forming a triangular face with the specified edge. 
    * Valid values for localFaceId are 0 and 1.
    * If the edge is adjacent to a hole in the mesh, i.e. there is only one adjacent face,
    * a localFaceId of 1 will return -1.
    */
    long int GetMissingPointIndex( IndexType edgeId, IndexType localFaceId ) 
      const 
    { 
      return m_EdgePoints[edgeId].missingVertex[localFaceId]; 
    }

    /** Returns the index of the point on the other side of the specified edge.
    * If the edge is not connected to the specified point, 
    * the method returns -1.
    */
    long int GetConnectedPointIndex( IndexType pointId, IndexType edgeId ) 
      const;
    
    /** Calculates and stores the indexing information for the mesh.
    * This function has to be called after all topological changes.
    */
    void BuildCellLinks();

    /** Copies all pointers of itk::Mesh (shallow copy) and all indexing 
    * information (deep copy). 
    */
    Pointer operator=( Pointer mesh );


  protected:

    /** Standard constructor. */
    IndexedTriangleMesh();

    /** Standard destructor. */
    ~IndexedTriangleMesh();

    /** Prints information about the object. */
    void PrintSelf( std::ostream& os, Indent indent ) const;

    /** Inserts a new edge between vertices v1 and v2 and returns the 
    * new edgeId.
    */
    IndexType InsertEdge( IndexType v1, IndexType v2 );

    /** Builds the indices for neighbouring faces. */
    void InitializeFaceLinks();


    /** Helper class which holds edge information of the mesh. */
    class E2VStruct
    {
    public:
      /** Equality operator */
      bool operator==( const E2VStruct& o ) const 
      { 
        return ( ( o.edgeVertex[0]==edgeVertex[0] && 
                   o.edgeVertex[1]==edgeVertex[1]    ) || 
                 ( o.edgeVertex[0]==edgeVertex[1] && 
                   o.edgeVertex[1]==edgeVertex[0]    )    ); 
      }

      /** Indices of the two points that form the edge */
      IndexType	edgeVertex[2];	    
      /** Indices of the two points that form triangles which the edge is 
      * part of (or -1 if there is no triangle)
      */
      long int  missingVertex[2];   
    };


  private:

    typedef std::vector<IndexType> V2EStruct;
    typedef std::vector<IndexType> IndexVector;
    typedef vnl_vector_fixed<IndexType,3> IndexVector3;


    /** Stores the number of edges in the mesh. */
    int						            m_NumberOfEdges;
    /** Stores the number of faces/triangles in the mesh. */
    int                       m_NumberOfFaces;
    /** Contains the point indices at and nearby each edge. */
    std::vector<E2VStruct>    m_EdgePoints;		          
    /** Contains a vector with edge indices for each point. */
    std::vector<V2EStruct>    m_EdgesPerVertex;	    
    /** Stores the three point indices for each face. */
    std::vector<IndexVector3> m_FaceVertices;       
    /** Stores all faces a specific point is part of. */
    std::vector<IndexVector>  m_FaceLinks;      
    /** Contains for each face a vector with indices of adjacent faces. */
    std::vector<IndexVector>  m_AdjacentFaces;  

  };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIndexedTriangleMesh.txx"
#endif

#endif
