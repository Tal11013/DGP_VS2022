#pragma once
#include "stdafx.h"
#include <list>

// Get the angle between 2 points using trig functions
double getAngleBetweenVertices(const MPoint& left, const MPoint& center, const MPoint& right);

// Get curvature of each vertex in the mesh
void getGaussianCurvature(const MFnMesh& meshFn, std::map<int, double>& curvature);

// Get connected components in the mesh. Specify onlyBoundaries=true to get
// only the connected components of the boundary.
int connectedComponents(const MFnMesh &meshFn, bool onlyBoundaries);

// Get a linked list representing all boundary vertices in the mesh.
void getBoundary(const MFnMesh& meshFn, std::list<int> &vertices);

// Create a UV set with given name and set it properly
void createUV(MFnMesh& meshFn, const MFloatArray& u, const MFloatArray& v, const char* name);
