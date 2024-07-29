#include "LSCMCmd.h"
#include "stdafx.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/GMM_Macros.h"
#include "helpers.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/MatlabInterface.h"


int isNotFixed(const int& index, const int2& fixed) {
	return index != fixed[0] && index != fixed[1];
}


// Get the indices of the 2 farthest points on the boundary, which will be
// fixed in the parametrization
void getFixedPoints(const MFnMesh& meshFn, int2& indices) {
	// Get boundary vertices
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


	// Compute the cumulative boundary length, to get the 2 vertices that
	// are the most apart from each other (meaning dist = totalBoundary / 2)
	indices[0] = *j;
	double cumulativeBoundaryLength = 0; 
	for (++i, ++j ; i != i_end && j != j_end; ++i, ++j) {
		cumulativeBoundaryLength += (*i).length() / totalBoundaryLength;
		if (cumulativeBoundaryLength > 0.499) {
			indices[1] = *j;
			break;
		}
	}
}


// Initialize a map between rows in coordinate matrix and indexes in the mesh,
// and the other way around
void createMaps(const MFnMesh& meshFn, MIntArray& rowMap, 
				std::map<int,int>& indexMap, const int2& fixed) {
	MItMeshVertex vertex_it(meshFn.object());
	int row = 0;
	while (!vertex_it.isDone()) {
		if (isNotFixed(vertex_it.index(), fixed)) {
			rowMap[row] = vertex_it.index();
			indexMap[vertex_it.index()] = row;
			++row;
		}
		vertex_it.next();
	}
}

// Calculate the weight matrix then solve it
void LSCMWeights(const MFnMesh& meshFn, MFloatArray& u, MFloatArray& v, 
					const MIntArray& rowMap, std::map<int,int>& indexMap, 
					const int2& fixed) {
	const unsigned int n = meshFn.numVertices(), m = 2;
	GMMSparseComplexRowMatrix weights(n - m ,n - m);
	GMMSparseComplexRowMatrix rhs(n - m, 1);



	MIntArray polygonVertices;
	MItMeshPolygon poly_it(meshFn.object());
	MItMeshVertex vertex_it(meshFn.object());
	while(!poly_it.isDone()) {
		// Get points and actual coordinates of the vertices of the triangle
		poly_it.getVertices(polygonVertices);
		const int first_v = polygonVertices[0], second_v = polygonVertices[1],
				  third_v = polygonVertices[2];

		MPoint first, second, third;
		meshFn.getPoint(first_v, first);
		meshFn.getPoint(second_v, second);
		meshFn.getPoint(third_v, third);

		// Get the cotangents of the angles of the triangle
		double area;
		poly_it.getArea(area);
		const double alpha = (second - first)  * (third  - first) / (2 * area),
					 beta  = (first  - third)  * (second - third) / (2 * area),
					 gamma = (first  - second) * (third  - second) / (2 * area);

		// Set the values in the matrix to the correct values
		// some ugly code but this is better than looping
		if (isNotFixed(first_v, fixed)) {
			weights(indexMap[first_v], indexMap[first_v]) += - beta - gamma;
			if (isNotFixed(second_v, fixed)) {
				weights(indexMap[first_v], indexMap[second_v]) += {beta, -1};
			} else {
				rhs(indexMap[first_v], 0) -= std::complex<double>(u[second_v], v[second_v])
												* std::complex<double>(beta, -1);
			}
			if (isNotFixed(third_v, fixed)) {
				weights(indexMap[first_v], indexMap[third_v]) += {gamma, 1};
			} else {
				rhs(indexMap[first_v], 0) -= std::complex<double>(u[third_v], v[third_v]) *
												  std::complex<double>(gamma, 1);
			}
		}

		if (isNotFixed(second_v, fixed)) {
			weights(indexMap[second_v], indexMap[second_v]) += - beta - alpha;
			if (isNotFixed(third_v, fixed)) {
				weights(indexMap[second_v], indexMap[third_v]) += {alpha, -1};
			} else {
				rhs(indexMap[second_v], 0) -= std::complex<double>(u[third_v], v[third_v])
												 * std::complex<double>(alpha, -1);
			}
			if (isNotFixed(first_v, fixed)) {
				weights(indexMap[second_v], indexMap[first_v]) += {beta, 1};
			} else {
				rhs(indexMap[second_v], 0) -= std::complex<double>(u[first_v], v[first_v])
												 * std::complex<double>(beta, 1);
			}
		}

		if (isNotFixed(third_v, fixed)) {
			weights(indexMap[third_v], indexMap[third_v]) += - gamma - alpha;
			if (isNotFixed(second_v, fixed)) {
				weights(indexMap[third_v], indexMap[second_v]) += {alpha, 1};
			} else {
				rhs(indexMap[third_v], 0) -= std::complex<double>(u[second_v], v[second_v])
												* std::complex<double>(alpha, 1);
			}
			if (isNotFixed(first_v, fixed)) {
				weights(indexMap[third_v], indexMap[first_v]) += {gamma, -1};
			} else {
				rhs(indexMap[third_v], 0) -= std::complex<double>(u[first_v], v[first_v])
												* std::complex<double>(gamma, -1);;
			}
		}
		poly_it.next();
	}

	int result = MatlabGMMDataExchange::SetEngineSparseMatrix("rhs", rhs);
	result = MatlabGMMDataExchange::SetEngineSparseMatrix("weights", weights);



	// Solve for the coordinates
	MatlabInterface::GetEngine().Eval("weights = weights * -1");
	MatlabInterface::GetEngine().Eval("rhs = rhs * -1");
	MatlabInterface::GetEngine().Eval("coords = solve_linear_system_with_cholesky(weights, rhs)");

		// Get the coords back from MatLab
	GMMDenseComplexColMatrix coords(n - m, 1);
	MatlabGMMDataExchange::GetEngineDenseMatrix("coords", coords);

	// Set the coords to u, v
	for (unsigned int currRow = 0; currRow < n - m; ++currRow) {
		const int index = rowMap[currRow];
		u[index] = coords(currRow, 0).real();
		v[index] = coords(currRow, 0).imag();
	}
}


LSCMCmd::LSCMCmd() = default;

void* LSCMCmd::creator() {
	return new LSCMCmd();
}

MString LSCMCmd::commandName() {
	return "LSCMCmd";
}

bool LSCMCmd::isUndoable() const {
	return false;
}

MStatus LSCMCmd::doIt(const MArgList& argList) {
	// Checks to see if everything is alright

	MStatus stat = MS::kSuccess;

	const MSyntax commandSyntax = syntax();

	MArgDatabase argData(commandSyntax, argList, &stat);
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
		int numPolygons = meshFn.numPolygons(&stat);

	MItMeshPolygon poly(meshObject);
	if(!poly.isPlanar(&stat) || poly.isLamina(&stat) || poly.isHoled(&stat)) {
		
		MCHECKERROR(MS::kFailure, "The given polygon shape is either self intersecting, holed or non-planar which are not supported");
	}

	for (int i=0; i<numPolygons; i++) {
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

	int2 fixed;
	getFixedPoints(meshFn, fixed);

	// Set u,v coords for the fixed vertices
	u[fixed[0]] = 0;
	v[fixed[0]] = 0;
	u[fixed[1]] = 1;
	v[fixed[1]] = 0;

	// Initialize the complex matrices for the weights and the RHS
	const int n = meshFn.numVertices(), m = 2;

	MIntArray rowMap(meshFn.numVertices() - 2);

	std::map<int, int> indexMap;
	createMaps(meshFn, rowMap, indexMap, fixed);


	LSCMWeights(meshFn, u, v, rowMap, indexMap, fixed);

	createUV(meshFn, u, v, "LSCM Conformal");

	return MS::kSuccess;

}


MSyntax LSCMCmd::syntax() {
	MSyntax commandSyntax;

	const MStatus stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");
	commandSyntax.useSelectionAsDefault(true);
	return commandSyntax;
}