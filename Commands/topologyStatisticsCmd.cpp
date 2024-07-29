#include "topologyStatisticsCmd.h"
#include "stdafx.h"
#include "Utils/Maya_Macros.h"
#include <Utils/Maya_Utils.h>
#include "..\helpers.h"

topologyStatisticsCmd::topologyStatisticsCmd() {
	
}

void* topologyStatisticsCmd::creator() {
    return new topologyStatisticsCmd();
}

MString topologyStatisticsCmd::commandName() {
    return "topologyStatisticsCmd";
}

bool topologyStatisticsCmd::isUndoable() const {
	return false;
}


MStatus topologyStatisticsCmd::doIt(const MArgList& argList) {
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


	MString message = "";
	message += "Mesh name: " + meshFn.name() + "\n";

	MItMeshPolygon poly(meshObject);
	bool isTriangular = true;
	for (int i=0; i<meshFn.numPolygons(); i++) {
		if (poly.polygonVertexCount() != 3) {
			isTriangular = false;
			break;
		}
		poly.next();
	}
	message += "Is triangle mesh: ";
	message += isTriangular ? "yes\n" : "no\n";

	message += "Number of vertices: " + 
		MString(std::to_string(meshFn.numVertices()).c_str()) + "\n";
	message += "Number of faces: " +
		MString(std::to_string(meshFn.numPolygons()).c_str()) + "\n";
	message += "Number of edges: " +
		MString(std::to_string(meshFn.numEdges()).c_str()) + "\n";

	const int components = connectedComponents(meshFn, false);
	const int boundaries = connectedComponents(meshFn, true);
	// x=V+F-E=2(C-g)-b
	const int eulerCharacteristic = meshFn.numVertices() + meshFn.numPolygons()
		- meshFn.numEdges();
	//g=C-((x+b)/2)
	const int genus = components - (eulerCharacteristic + boundaries) / 2;
	message += "Genus: " + MString(std::to_string(genus).c_str()) + "\n";
	message += "Number of connected components: " + 
		MString(std::to_string(components).c_str()) + "\n";
	message += "Number of boundaries: " + 
		MString(std::to_string(boundaries).c_str()) + "\n";
	message += "Euler characteristic: " + 
				MString(std::to_string(eulerCharacteristic).c_str()) + "\n";

	std::map <int, double> curvatures;
	getGaussianCurvature(meshFn, curvatures);
	double totalCurvature = 0;
	for (const auto& curvature : curvatures)
	{
		totalCurvature += curvature.second;
	}
	
	message += "Euler characteristic based on discrete Gauss-Bonnet: ";
	MString curvatureString;
	curvatureString.set((totalCurvature) / (2 * M_PI), 10);
	message += curvatureString + "\n";

	if (!isTriangular)
	{
		message += "Note that the mesh isn't triangular, so due to technical"
			 " constraints the value will be off \n";
	}

	MGlobal::displayInfo(message);

	return MS::kSuccess;
}

MSyntax topologyStatisticsCmd::syntax() {
	MStatus stat = MS::kSuccess;
	MSyntax commandSyntax;

	stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");
	commandSyntax.useSelectionAsDefault(true);

	return commandSyntax;
}