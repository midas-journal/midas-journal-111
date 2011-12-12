#ifndef _itkIndexedTriangleMesh_txx
#define _itkIndexedTriangleMesh_txx

#include "itkIndexedTriangleMesh.h"
#include <itkVertexCell.h>


namespace itk
{

  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::IndexedTriangleMesh()
  {
    m_NumberOfEdges = 0;
    m_NumberOfFaces = 0;
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::~IndexedTriangleMesh()
  {
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  void
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::PrintSelf( std::ostream& os, Indent indent ) const
  {
    Superclass::PrintSelf(os, indent);
    os << indent << "Number Of Edges: " << m_NumberOfEdges << std::endl;
    os << indent << "Number Of Faces: " << m_NumberOfFaces << std::endl;
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  typename IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>::IndexType
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::GetEdgeIndexForFace( IndexType faceId, IndexType localEdgeId ) const 
  {
    switch( localEdgeId ) 
    {
      case 0: return (IndexType)(this->GetConnectingEdgeIndex( 
                      m_FaceVertices[faceId][0], m_FaceVertices[faceId][1] ));
      case 1: return (IndexType)(this->GetConnectingEdgeIndex( 
                      m_FaceVertices[faceId][1], m_FaceVertices[faceId][2] ));
      case 2: return (IndexType)(this->GetConnectingEdgeIndex( 
                      m_FaceVertices[faceId][2], m_FaceVertices[faceId][0] ));
    }
    return 0;
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  long int
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::GetConnectingEdgeIndex( IndexType pointId0, IndexType pointId1 ) const 
  {
    for (unsigned int i=0; i<this->GetNumberOfEdgesForPoint( pointId0 ); i++) 
    {
      IndexType edgeId = this->GetEdgeIndexForPoint( pointId0, i );
      if (this->GetConnectedPointIndex( pointId0, edgeId ) == pointId1) 
      { 
        return edgeId;
      }
    }
    return -1;
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  long int
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::GetConnectedPointIndex( IndexType pointId, IndexType edgeId ) const
  { 
    if (m_EdgePoints[edgeId].edgeVertex[0]==pointId) 
    {
      return m_EdgePoints[edgeId].edgeVertex[1];
    }
    else if (m_EdgePoints[edgeId].edgeVertex[1]==pointId) 
    {
      return m_EdgePoints[edgeId].edgeVertex[0];
    }
    return -1;
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  void
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::BuildCellLinks()
  {
    Superclass::BuildCellLinks();
    typename MeshType::CellsContainerPointer cells = this->GetCells();
    if (cells.IsNull()) return;

    // clear stored data:
    m_NumberOfEdges = 0;
    m_NumberOfFaces = 0;
    m_EdgePoints.clear();
    m_EdgesPerVertex.clear();
    m_FaceVertices.clear();
    m_FaceLinks.clear();
    m_AdjacentFaces.clear();

    itkDebugMacro( "Indexing mesh data..." );
    typedef typename MeshType::CellsContainer::ConstIterator CellIterator;
    typedef typename MeshType::CellType CellType;
    typedef typename CellType::CellAutoPointer CellAutoPointer;
    typedef VertexCell<CellType> VertexType;

    m_EdgesPerVertex.resize( this->GetNumberOfPoints() );
    // extract basic mesh edge data:
    typename MeshType::CellIdentifier cId = 0;
    CellIterator cellIterator = this->GetCells()->Begin();
    CellIterator cellEnd      = this->GetCells()->End();
    while( cellIterator != cellEnd ) 
    {
      CellType *cell = cellIterator.Value();
      if (!cell) 
      {
        cId++;
        cellIterator++;
        continue;
      }
      if (cell->GetType() == CellType::LINE_CELL) 
      {
        int vertexId[2];
        for (int v=0; v<2;v++) 
        {
          typename MeshType::CellIdentifier neighbourId;
          bool found = GetBoundaryAssignment( 0, cId, v, &neighbourId );
          if (!found) { itkExceptionMacro( "Missing boundary information!" ); }
          CellAutoPointer neighbourCell;
          GetCell( neighbourId, neighbourCell );
          typename CellType::PointIdConstIterator pntId = 
                                                  neighbourCell->GetPointIds();
          vertexId[v] = (*pntId);
        }
        // init edge struct
        E2VStruct edgeData;
        edgeData.edgeVertex[0] = vertexId[0];	
        edgeData.edgeVertex[1] = vertexId[1];
        edgeData.missingVertex[0] = -1;			
        edgeData.missingVertex[1] = -1;
        // check if edge already in mesh
        bool edgeNew = true;
        for (int i=0; i<m_NumberOfEdges; i++) 
        {
          if (edgeData == m_EdgePoints[i]) 
          {
            edgeNew = false;
            break;
          }
        }
        if (edgeNew) 
        {
          m_EdgePoints.push_back( edgeData );
          m_EdgesPerVertex[vertexId[0]].push_back( m_NumberOfEdges );
          m_EdgesPerVertex[vertexId[1]].push_back( m_NumberOfEdges );
          m_NumberOfEdges++;
        }
      }
      cId++;
      cellIterator++;
    }

    // missing vertex data:
    cellIterator = this->GetCells()->Begin();
    cellEnd      = this->GetCells()->End();
    while( cellIterator != cellEnd ) 
    {
      CellType *cell = cellIterator.Value();
      if (cell->GetType() == CellType::TRIANGLE_CELL) 
      {
        int vertexId[3];
        typename CellType::PointIdConstIterator pntId = cell->GetPointIds();
        for (int v=0; v<3; v++) 
        {
          vertexId[v] = *pntId;  
          pntId++;
        }
        for (int ec=0; ec<3; ec++) 
        {
          int v1 = vertexId[ec];
          int v2 = vertexId[(ec+1)%3];
          int e = GetConnectingEdgeIndex( v1, v2 );
          if (e<0) e = InsertEdge( v1, v2 );
          int missingId = vertexId[(ec+2)%3];
          int missing0 = this->GetMissingPointIndex( e, 0 );
          if ( missing0 < 0) m_EdgePoints[e].missingVertex[0] = missingId;
          else if (missingId != missing0) 
          {
            int missing1 = m_EdgePoints[e].missingVertex[1];
            if (missing1 >= 0 && missing1 != missingId)             
            {
              itkExceptionMacro( "Mesh geometry is non-manifold!" );
            }
            m_EdgePoints[e].missingVertex[1] = missingId;
          }
        }
        m_NumberOfFaces++;
      }
      cellIterator++;
    }
    itkDebugMacro( "Mesh contains " << this->GetNumberOfPoints() << 
                   " vertices, " << m_NumberOfEdges << " edges, " << 
                   m_NumberOfFaces << " faces." << std::endl );

    this->InitializeFaceLinks();
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  typename IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>::IndexType
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::InsertEdge( IndexType v1, IndexType v2 )
  {
    E2VStruct edgeData;
    edgeData.edgeVertex[0] = v1;	
    edgeData.edgeVertex[1] = v2;
    edgeData.missingVertex[0] = -1;			
    edgeData.missingVertex[1] = -1;
    m_EdgePoints.push_back( edgeData );
    m_EdgesPerVertex[v1].push_back( m_NumberOfEdges );
    m_EdgesPerVertex[v2].push_back( m_NumberOfEdges );
    m_NumberOfEdges++;
    return (IndexType)(m_NumberOfEdges-1);
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  void
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::InitializeFaceLinks()
  {
    itkDebugMacro( "Computing normal vectors..." );

    typedef typename MeshType::CellsContainer::ConstIterator CellIterator;
    typedef typename MeshType::CellType CellType;

    m_FaceLinks.resize( this->GetNumberOfPoints() );
    m_FaceVertices.resize( m_NumberOfFaces );
    
    // init normals for all triangles:
    CellIterator cellIterator = this->GetCells()->Begin();
    CellIterator cellEnd      = this->GetCells()->End();
    int faceId = 0;
    while( cellIterator != cellEnd ) 
    {
      CellType *cell = cellIterator.Value();
      if (cell->GetType() == CellType::TRIANGLE_CELL) 
      {
        int vertexId[3];
        typename CellType::PointIdConstIterator pntId = cell->GetPointIds();
        for (int v=0; v<3; v++) 
        {
          vertexId[v] = *pntId;  
          m_FaceVertices[faceId][v] = *pntId;
          m_FaceLinks[*pntId].push_back( faceId );
          pntId++;
        }
        faceId++;
      }
      cellIterator++;
    }

    m_AdjacentFaces.resize( m_NumberOfFaces );
    for (int f=0; f<m_NumberOfFaces; f++) 
    {
      for (int v=0; v<3; v++) 
      {
        unsigned int vertexId = m_FaceVertices[f][v];
        for (unsigned int i=0; i<m_FaceLinks[vertexId].size(); i++) 
        {
          int adjFaceId = m_FaceLinks[vertexId][i];
          if (adjFaceId != f) m_AdjacentFaces[f].push_back( adjFaceId );
        }
      }
    }
    itkDebugMacro( " Done!" << std::endl );
  }


  template <typename TPixelType, typename TIndex, typename TCoordRepType>
  typename IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>::Pointer
  IndexedTriangleMesh<TPixelType, TIndex, TCoordRepType>
  ::operator=( Pointer mesh )
  {
    // copy PointSet data
    this->SetPoints( mesh->GetPoints() );
    this->SetPointData( mesh->GetPointData() );
    // copy Mesh data
    this->SetCells( mesh->GetCells() );
    this->SetCellLinks( mesh->GetCellLinks() );
    this->SetCellData( mesh->GetCellData() );
    this->SetCellsAllocationMethod( mesh->GetCellsAllocationMethod() );
    for( unsigned int dim=0; dim<MeshType::MaxTopologicalDimension; dim++ )
    {
      this->SetBoundaryAssignments( dim, mesh->GetBoundaryAssignments( dim ) );
    }
    // copy IndexedTriangleMesh data
    m_NumberOfEdges = mesh->GetNumberOfEdges();
    m_NumberOfFaces = mesh->GetNumberOfFaces();
    // copy EdgesPerVertex
    m_EdgesPerVertex.resize( mesh->GetNumberOfPoints() );
    for (unsigned int i=0; i<m_EdgesPerVertex.size(); i++)
    {
      m_EdgesPerVertex[i].resize( mesh->GetNumberOfEdgesForPoint( i ) );
      for (unsigned int j=0; j<m_EdgesPerVertex[i].size(); j++)
      {
        m_EdgesPerVertex[i][j] = mesh->GetEdgeIndexForPoint( i, j );
      }
    }
    // copy EdgePoints
    m_EdgePoints.resize( mesh->GetNumberOfEdges() );
    for (unsigned int i=0; i<m_EdgePoints.size(); i++)
    {
      m_EdgePoints[i].edgeVertex[0] = mesh->GetPointIndexForEdge( i, 0 );
      m_EdgePoints[i].edgeVertex[1] = mesh->GetPointIndexForEdge( i, 1 );
      m_EdgePoints[i].missingVertex[0] = mesh->GetMissingPointIndex( i, 0 );
      m_EdgePoints[i].missingVertex[1] = mesh->GetMissingPointIndex( i, 1 );
    }
    // copy FaceVertices
    m_FaceVertices.resize( mesh->GetNumberOfFaces() );
    for (unsigned int i=0; i<m_FaceVertices.size(); i++)
    {
      for (unsigned int j=0; j<3; j++) 
      {
        m_FaceVertices[i][j] = mesh->GetPointIndexForFace( i, j );
      }
    }
    // copy FaceLinks
    m_FaceLinks.resize( mesh->GetNumberOfPoints() );
    for (unsigned int i=0; i<m_FaceLinks.size(); i++)
    {
      m_FaceLinks[i].resize( mesh->GetNumberOfFacesForPoint( i ) );
      for (unsigned int j=0; j<m_FaceLinks[i].size(); j++)
      {
        m_FaceLinks[i][j] = mesh->GetFaceIndexForPoint( i, j );
      }
    }
    // copy AdjacentFaces
    m_AdjacentFaces.resize( mesh->GetNumberOfFaces() );
    for (unsigned int i=0; i<m_AdjacentFaces.size(); i++)
    {
      m_AdjacentFaces[i].resize( mesh->GetNumberOfAdjacentFaces( i ) );
      for (unsigned int j=0; j<m_AdjacentFaces[i].size(); j++)
      {
        m_AdjacentFaces[i][j] = mesh->GetAdjacentFaceIndex( i, j );
      }
    }
    return this;
  }

}

#endif
