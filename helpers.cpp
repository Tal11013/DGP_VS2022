#include "stdafx.h"
#include "helpers.h"
#include <queue>

double getAngleBetweenVertices(const MPoint& left, const MPoint& center, const MPoint& right) {
	const MVector u = center - left;
	const MVector v = center - right;

	return u.angle(v);
}

void getGaussianCurvature(const MFnMesh& meshFn, std::map<int, double>& curvature)
{
	MItMeshVertex vertex_it{meshFn.object()};
	MItMeshPolygon face_it = meshFn.object();
	while (!vertex_it.isDone())
	{
		if (vertex_it.onBoundary())
		{
			curvature[vertex_it.index()] = M_PI;
		}
		else
		{
			curvature[vertex_it.index()] = 2 * M_PI;
		}
		vertex_it.next();
	}

	MIntArray vertices;
	MPointArray points;
	while (!face_it.isDone())
	{
		face_it.getVertices(vertices);
		const int numVertices = vertices.length();
		points.setLength(numVertices);
		for (int i = 0; i < numVertices; ++i)
		{
			meshFn.getPoint(vertices[i], points[i]);
		}
		for (int i = 0; i < numVertices; ++i)
		{
			curvature[vertices[i]] -= getAngleBetweenVertices(
				points[(i + numVertices - 1) % numVertices],
				points[i], points[(i + 1) % numVertices]);
		}
		face_it.next();
	}
}

int connectedComponents(const MFnMesh &meshFn,
                        bool onlyBoundaries) {
	// Note that an edge is a boundary edge iff all vertices it is connected to
	// are boundary vertices, but a boundary vertex may be connected to non-boundary
	// edges.

	std::queue<int> edgeQueue;
	std::map<int, bool> visitedEdges;

	int components = 0, currIndex, prevIndex, unvisited = 0;
	MIntArray connectedEdges;

	MItMeshEdge edge_it = meshFn.object();


	while (!edge_it.isDone())
	{
		currIndex = edge_it.index();
		if (!onlyBoundaries || edge_it.onBoundary()) {
			++unvisited;
			visitedEdges[currIndex] = false;
		}
		edge_it.next();
	}


	edge_it.reset();

	while (unvisited > 0) {
		++components;

		int firstEdge;
		for (const std::pair<const int, bool>& visited_edge : visitedEdges) {
			if (visited_edge.second == false){
				firstEdge = visited_edge.first;
			 	break;
			}
		}

		edgeQueue.push(firstEdge);
		visitedEdges[firstEdge] = true;
		--unvisited;

		while(!edgeQueue.empty()) {
			currIndex = edgeQueue.front();
			edgeQueue.pop();

			edge_it.setIndex(currIndex, prevIndex);
			edge_it.getConnectedEdges(connectedEdges);

			for (int edge : connectedEdges) {
				edge_it.setIndex(edge, prevIndex);
				if ((!onlyBoundaries || edge_it.onBoundary()) && visitedEdges[edge] == false) {
					visitedEdges[edge] = true;
					edgeQueue.push(edge);
					--unvisited;
				}
			}
		}
	}

	return components;
}


// Gets the boundary of the mesh as a linked list of vertex indices.
// Assumes that the boundary of a mesh is isomorphic to a plane polygon, otherwise
// it wouldn't work.
void getBoundary(const MFnMesh& meshFn,	std::list<int> &vertices) {
	vertices.clear();

	int startIndex, prevIndex = 0, currIndex, _;
	MIntArray connectedVertices;
	MItMeshVertex vertex_it = meshFn.object();
	while (!vertex_it.isDone()) {
		if (vertex_it.onBoundary()) {
			currIndex = prevIndex = startIndex = vertex_it.index();
			vertices.push_back(startIndex);
			break;
		}
		vertex_it.next();
	}

	do {
		currIndex = vertex_it.index();
		vertex_it.getConnectedVertices(connectedVertices);
		for (int vertex : connectedVertices) {
			vertex_it.setIndex(vertex, _);
			if (vertex != prevIndex && vertex_it.onBoundary()) {
				prevIndex = currIndex;
				currIndex = vertex;
				vertices.push_back(vertex);
				break;
			}
		}
	} while (currIndex != startIndex);
}

void createUV(MFnMesh& meshFn, const MFloatArray& u, const MFloatArray& v, const char* name) {
	const MString uvName = name;
	meshFn.createUVSetWithName(uvName);
	meshFn.setUVs(u, v, &uvName);
	MIntArray uvCounts, uvIds;
	meshFn.getVertices(uvCounts, uvIds);
	meshFn.assignUVs(uvCounts, uvIds, &uvName);
}
