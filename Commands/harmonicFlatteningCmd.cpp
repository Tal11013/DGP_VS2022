#include "harmonicFlatteningCmd.h"
#include "stdafx.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/GMM_Macros.h"
#include "helpers.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/MatlabInterface.h"

// Calculates the (u,v) coordinates of the boundary vertices by chord length.
unsigned int calculateBoundaryCoords(const MFnMesh& meshFn, MFloatArray& u, MFloatArray& v) {
	std::list<int> boundary;
	getBoundary(meshFn, boundary);

	// Get the differences between the coords of pairs of vertices on the boundary
	std::list<MPoint> boundaryCoords = {};
	MPoint coords;
	for (int vertex : boundary) {
		meshFn.getPoint(vertex, coords);
		boundaryCoords.emplace_back(coords);
	}
	std::list<MVector> differences = {};
	std::adjacent_difference(boundaryCoords.begin(), boundaryCoords.end(), 
							std::back_inserter(differences));


	// Sum the lengths of the differences to get the total boundary length
	double totalBoundaryLength = 0;
	auto i = differences.begin();
	for (++i; i != differences.end(); ++i) {
		totalBoundaryLength += (*i).length();
	}

	// Initialize iterators to go over the boundaries except the first and last
	auto j = boundary.begin(), j_end = boundary.end();
	auto i_end = differences.end();
	--j_end;
	--i_end;
	i = differences.begin();

	// Set the boundary coords, at last
	u[*j] = 1;
	v[*j] = 0;
	double cumulativeBoundaryLength = 0; 
	for (++i, ++j ; i != i_end && j != j_end; ++i, ++j) {
		cumulativeBoundaryLength += (*i).length() / totalBoundaryLength;
		u[*j] = cos(2 * M_PI * cumulativeBoundaryLength);
		v[*j] = sin(2 * M_PI * cumulativeBoundaryLength);
	}

	return meshFn.numVertices() - boundary.size() + 1;
}

// Initialize a map between rows in coordinate matrix and indexes in the mesh,
// and the other way around
void createMaps(const MFnMesh& meshFn, MIntArray& rowMap,
                std::map<int, int>& indexMap) {
	MItMeshVertex vertex_it(meshFn.object());
	int row = 0;
	while (!vertex_it.isDone()) {
		if (!vertex_it.onBoundary()) {
			rowMap[row] = vertex_it.index();
			indexMap[vertex_it.index()] = row;
			++row;
		}
		vertex_it.next();
	}
}


// Creates (u,v) coordinates for all vertices using a uniform-weight matrix,
// Assuming the coordinates of the boundary are fixed and a row map is known
void uniformWeights(const MFnMesh& meshFn, MFloatArray& u, MFloatArray& v, 
					const MIntArray& rowMap, std::map<int,int>& indexMap, 
					unsigned int rowCount) {
	GMMSparseRowMatrix weight_matrix(rowCount, rowCount);
	GMMDenseColMatrix rhs(rowCount, 2);

	MItMeshVertex vertex_it(meshFn.object());

	
	// Fill in the weight matrix and the rhs vector with uniform weights
	MIntArray connectedVertices;
	for (unsigned int currRow = 0; currRow < rowCount; ++currRow) {
		int _;
		int rowSum = 0;
		vertex_it.setIndex(rowMap[currRow], _);
		vertex_it.getConnectedVertices(connectedVertices);
		for (const int vertex : connectedVertices) {
			vertex_it.setIndex(vertex, _);
			++rowSum;

			if (vertex_it.onBoundary()) {
				rhs(currRow, 0) -= u[vertex];
				rhs(currRow, 1) -= v[vertex];
			} else {
				const unsigned int col = indexMap[vertex];
				weight_matrix(currRow, col) = 1;
			}
		}
		weight_matrix(currRow, currRow) = -rowSum;
	}

	// Transfer matrices to MatLab
	int result = MatlabGMMDataExchange::SetEngineDenseMatrix("rhs", rhs);
	result = MatlabGMMDataExchange::SetEngineSparseMatrix("weights", weight_matrix);



	// Solve for the coordinates
	MatlabInterface::GetEngine().Eval("weights = weights * -1");
	MatlabInterface::GetEngine().Eval("rhs = rhs * -1");
	MatlabInterface::GetEngine().Eval("coords = solve_linear_system_with_cholesky(weights, rhs)");

	// Get the coords back from MatLab
	GMMDenseColMatrix coord_matrix(rowCount, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("coords", coord_matrix);

	// Set the coords to u, v
	for (unsigned int currRow = 0; currRow < rowCount; ++currRow) {
		const int index = rowMap[currRow];
		const double uValue = coord_matrix(currRow, 0);
		const double vValue = coord_matrix(currRow, 1);
		u[index] = uValue;
		v[index] = vValue;
	}
}

// Creates (u,v) coordinates for all vertices using a cotangent-weight matrix,
// Assuming the coordinates of the boundary are fixed and a row map is known
void cotangentWeights(const MFnMesh& meshFn, MFloatArray& u, MFloatArray& v, 
					const MIntArray& rowMap, std::map<int,int>& indexMap, 
					unsigned int rowCount) {
	GMMSparseRowMatrix weight_matrix(rowCount, rowCount);
	GMMDenseColMatrix rhs(rowCount, 2);


	MItMeshVertex vertex_it(meshFn.object());
	MItMeshEdge edge_it(meshFn.object());
	MItMeshPolygon polygon_it(meshFn.object());

	MIntArray connectedEdges, connectedFaces, polygonVertices;

	// Initialize the matrix
	for (unsigned int currRow = 0; currRow < rowCount; ++currRow) {
		const int currVertex = rowMap[currRow];

		int _;
		double rowSum = 0;
		vertex_it.setIndex(currVertex, _);
		vertex_it.getConnectedEdges(connectedEdges);
		for (const int edge : connectedEdges) {
			double value = 0;

			// Get the vertex on the other side of the edge
			int adjVertex, thirdVertex;
			int2 edgeVertices;
			meshFn.getEdgeVertices(edge, edgeVertices);
			if (edgeVertices[0] == currVertex) {
				adjVertex = edgeVertices[1];
			} else {
				adjVertex = edgeVertices[0];
			}


			edge_it.setIndex(edge, _);
			edge_it.getConnectedFaces(connectedFaces);
			for (const int face : connectedFaces) {
				polygon_it.setIndex(face, _);
				polygon_it.getVertices(polygonVertices);

				// Get the third polygon vertex
				for (const int vertex : polygonVertices) {
					if (vertex != currVertex && vertex != adjVertex) {
						thirdVertex = vertex;
					}
				}
				// Get the physical points on the triangle
				MPoint left, center, right;
				meshFn.getPoint(adjVertex, left);
				meshFn.getPoint(thirdVertex, center);
				meshFn.getPoint(currVertex, right);

				// Calculate the value to add to the sum
				value += (center - left) * (center - right) / (center - left ^ center - right).length();
			}

			
			// If the 2nd vertex is on the boundary, subtract the sum from the rhs,
			// and if not set it as the matrix element
			vertex_it.setIndex(adjVertex, _);
			if (vertex_it.onBoundary()) {
				rhs(indexMap[currVertex], 0) -= u[adjVertex] * value;
				rhs(indexMap[currVertex], 1) -= v[adjVertex] * value;
			} else {
				weight_matrix(indexMap[currVertex], indexMap[adjVertex]) = value;
			}

			rowSum += value;
		}
		weight_matrix(currRow, currRow) = -rowSum;
	}

	// Transfer matrices to MatLab
	int result = MatlabGMMDataExchange::SetEngineDenseMatrix("rhs", rhs);
	result = MatlabGMMDataExchange::SetEngineSparseMatrix("weights", weight_matrix);

	// Solve for the coordinates
	MatlabInterface::GetEngine().Eval("weights = weights * -1");
	MatlabInterface::GetEngine().Eval("rhs = rhs * -1");
	MatlabInterface::GetEngine().Eval("coords = solve_linear_system_with_cholesky(weights, rhs)");

	// Get the coords back from MatLab
	GMMDenseColMatrix coord_matrix(rowCount, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("coords", coord_matrix);

	// Set the coords to u, v
	for (unsigned int currRow = 0; currRow < rowCount; ++currRow) {
		const int index = rowMap[currRow];
		const double uValue = coord_matrix(currRow, 0);
		const double vValue = coord_matrix(currRow, 1);
		u[index] = uValue;
		v[index] = vValue;
	}
}

// Creates (u,v) coordinates for all vertices using a mean-value-weight matrix,
// Assuming the coordinates of the boundary are fixed and a row map is known
void meanValueWeights(const MFnMesh& meshFn, MFloatArray& u, MFloatArray& v, 
					const MIntArray& rowMap, std::map<int,int>& indexMap, 
					unsigned int rowCount) {
	GMMSparseRowMatrix weight_matrix(rowCount, rowCount);
	GMMDenseColMatrix rhs(rowCount, 2);

	MItMeshVertex vertex_it(meshFn.object());
	MItMeshEdge edge_it(meshFn.object());
	MItMeshPolygon polygon_it(meshFn.object());

	MIntArray connectedEdges, connectedFaces, polygonVertices;

	// Initialize the matrix
	for (unsigned int currRow = 0; currRow < rowCount; ++currRow) {
		const int currVertex = rowMap[currRow];

		int _;
		double rowSum = 0;
		vertex_it.setIndex(currVertex, _);
		vertex_it.getConnectedEdges(connectedEdges);
		for (const int edge : connectedEdges) {
			double value = 0;

			// Get the vertex on the other side of the edge
			int adjVertex, thirdVertex;
			int2 edgeVertices;
			meshFn.getEdgeVertices(edge, edgeVertices);
			if (edgeVertices[0] == currVertex) {
				adjVertex = edgeVertices[1];
			} else {
				adjVertex = edgeVertices[0];
			}


			edge_it.setIndex(edge, _);
			edge_it.getConnectedFaces(connectedFaces);
			for (const int face : connectedFaces) {
				polygon_it.setIndex(face, _);
				polygon_it.getVertices(polygonVertices);

				// Get the third polygon vertex
				for (const int vertex : polygonVertices) {
					if (vertex != currVertex && vertex != adjVertex) {
						thirdVertex = vertex;
					}
				}
				// Get the physical points on the triangle
				MPoint left, center, right;
				meshFn.getPoint(adjVertex, right);
				meshFn.getPoint(thirdVertex, left);
				meshFn.getPoint(currVertex, center);


				// tan = (|u||v| - u*v ) / |u x v|

				// Calculate the value to add to the sum
				value += ((center - left).length() * (center - right).length()
						- (center - left) * (center - right)) / 
					((center - left ^ center - right).length() *
							(center - right).length());
			}

			
			// If the 2nd vertex is on the boundary, subtract the sum from the rhs,
			// and if not set it as the matrix element
			vertex_it.setIndex(adjVertex, _);
			if (vertex_it.onBoundary()) {
				rhs(indexMap[currVertex], 0) -= u[adjVertex] * value;
				rhs(indexMap[currVertex], 1) -= v[adjVertex] * value;
			} else {
				weight_matrix(indexMap[currVertex], indexMap[adjVertex]) = value;
			}

			rowSum += value;
		}
		weight_matrix(currRow, currRow) = -rowSum;
	}

	// Transfer matrices to MatLab
	int result = MatlabGMMDataExchange::SetEngineDenseMatrix("rhs", rhs);
	result = MatlabGMMDataExchange::SetEngineSparseMatrix("weights", weight_matrix);


	// Solve for the coordinates, using LU as the matrix isn't symmetric positive semi-definite
	MatlabInterface::GetEngine().Eval("coords = solve_linear_system_with_LU(weights, rhs)");

	// Get the coords back from MatLab
	GMMDenseColMatrix coord_matrix(rowCount, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("coords", coord_matrix);

	// Set the coords to u, v
	for (unsigned int currRow = 0; currRow < rowCount; ++currRow) {
		const int index = rowMap[currRow];
		const double uValue = coord_matrix(currRow, 0);
		const double vValue = coord_matrix(currRow, 1);
		u[index] = uValue;
		v[index] = vValue;
	}
}

harmonicFlatteningCmd::harmonicFlatteningCmd() = default;

void* harmonicFlatteningCmd::creator() {
	return new harmonicFlatteningCmd();
}

MString harmonicFlatteningCmd::commandName() {
	return "harmonicFlatteningCmd";
}

bool harmonicFlatteningCmd::isUndoable() const {
	return false;
}

MStatus harmonicFlatteningCmd::doIt(const MArgList& argList) {
	// Checks to see if everything is alright

	MStatus stat = MS::kSuccess;

	const MSyntax commandSyntax = syntax();

	const MArgDatabase argData(commandSyntax, argList, &stat);
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());

	MSelectionList objectsList;
	stat = argData.getObjects(objectsList);
	MCHECKERROR(stat, "Can't access object list");

	MObject object;
	stat = objectsList.getDependNode(0, object);
	MCHECKERROR(stat, "Can't access object");

	MObject meshObject;
	stat = Maya_Utils::getMe_a_Mesh(object, meshObject);
	MCHECKERROR(stat, "Object is not a mesh");

	MFnMesh meshFn(meshObject, &stat);
	MCHECKERROR(stat, "Can't access mesh");

	// Check that the mesh is indeed a triangle mesh
	const int numPolygons = meshFn.numPolygons(&stat);

	MItMeshPolygon poly(meshObject);
	if(!poly.isPlanar(&stat) || poly.isLamina(&stat) || poly.isHoled(&stat))
	{
		
		MCHECKERROR(MS::kFailure, "The given polygon shape is either self intersecting, holed or non-planar which are not supported");
	}

	for (int i=0; i<numPolygons; i++)
	{
		const unsigned int temp = poly.polygonVertexCount();
		if ( 3 != temp )
			MCHECKERROR(MS::kFailure, "this is not a triangle mesh!");
		poly.next();
	}

	// Check that the mesh is a topological disk
	const int components = connectedComponents(meshFn, false);
	const int boundaries = connectedComponents(meshFn, true);
	const int eulerCharacteristic = meshFn.numVertices() + meshFn.numPolygons()
		- meshFn.numEdges();
	const int genus = components - (eulerCharacteristic + boundaries) / 2;

	if (components != 1 || boundaries != 1 || genus != 0) {
		MGlobal::displayError("mesh is not a topological disk");
		return MS::kFailure;
	}

	MFloatArray u, v;
	u.setLength(meshFn.numVertices());
	v.setLength(meshFn.numVertices());

	const unsigned int rowCount = calculateBoundaryCoords(meshFn, u, v);

	// Create 3 sets of u,v arrays, copy to them the already-initialized
	// boundary values. This should take quite some time
	MFloatArray unif_u, unif_v, cot_u, cot_v, mean_u, mean_v;
	unif_u.copy(u);
	unif_v.copy(v);
	cot_u.copy(u);
	cot_v.copy(v);
	mean_u.copy(u);
	mean_v.copy(v);

	MIntArray rowMap(rowCount);
	std::map<int, int> indexMap;
	createMaps(meshFn, rowMap, indexMap);

	// Create the 3 (u,v) sets
	uniformWeights(meshFn, unif_u, unif_v, rowMap, indexMap, rowCount);
	cotangentWeights(meshFn, cot_u, cot_v, rowMap, indexMap, rowCount);
	meanValueWeights(meshFn, mean_u, mean_v, rowMap, indexMap, rowCount);

	// Assign the 3 created (u,v) mappings to the mesh
	createUV(meshFn, unif_u, unif_v, "Uniform Weight Harmonic");
	createUV(meshFn, cot_u, cot_v, "Cotangent Weight Harmonic");
	createUV(meshFn, mean_u, mean_v, "Mean Value Weight Harmonic");

	return MS::kSuccess;
}

MSyntax harmonicFlatteningCmd::syntax() {
	MSyntax commandSyntax;

	const MStatus stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");
	commandSyntax.useSelectionAsDefault(true);
	return commandSyntax;
}