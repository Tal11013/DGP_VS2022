#include "stdafx.h"

#include "Utils/STL_Macros.h"
#include "TriangulatePolygonCmd.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/triangulation.h"


//flags
#define kNumVerticesShort       "-v"
#define kNumVerticesLong        "-numWantedVertices"
#define kSubsampleBoundaryShort "-s"
#define kSubsampleBoundaryLong  "-subsampleBoundary"



TriangulatePolygonCmd::TriangulatePolygonCmd()
{

}


void* TriangulatePolygonCmd::creator()
{
	return new TriangulatePolygonCmd;
}



MString TriangulatePolygonCmd::commandName()
{
	return "triangulatePolygon"; 
}



//usage example:  triangulatePolygon -numWantedVertices 1000 -subsampleBoundary false MeshShape1;
//-numWantedVertices is optional
//-subsampleBoundary is optional - default is true
MStatus TriangulatePolygonCmd::doIt(const MArgList& argList)
{
	MStatus stat = MS::kSuccess; 

	MSyntax commandSyntax = syntax();
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

	int numVerticesInMeshCreated = 0;
	int numPolygons = meshFn.numPolygons(&stat);

	if(numPolygons != 1)
	{
		MCHECKERROR(MS::kFailure, "Please select a mesh with a single polygon");
	}

	MItMeshPolygon poly(meshObject);
	if(!poly.isPlanar(&stat) || poly.isLamina(&stat) || poly.isHoled(&stat))
	{
		MCHECKERROR(MS::kFailure, "The given polygon shape is either self intersecting, holed or non-planar which are not supported");
	}

	int numWantedVertices = 0;

	if(argData.isFlagSet(kNumVerticesLong))
	{
		numWantedVertices = argData.flagArgumentInt(kNumVerticesLong, 0, &stat);
		MCHECKERROR(stat, "Can't access flag " + MString(kNumVerticesLong) + " argument");
		if(numWantedVertices < 3)
		{
			MCHECKERROR(MS::kFailure, "Illegal number of wanted vertices has been set");
		}
	}
	else //calculating a default value
	{
		numWantedVertices = 2*meshFn.polygonVertexCount(0, &stat);
		MCHECKERROR(stat, "Failed to get the number of vertices in the given polygon");
	}

	bool subsampleBoundary = true;
	if(argData.isFlagSet(kSubsampleBoundaryLong))
	{
		subsampleBoundary = argData.flagArgumentBool(kSubsampleBoundaryLong, 0, &stat);
		MCHECKERROR(stat, "Can't access flag " + MString(kSubsampleBoundaryLong) + " argument");
	}
	
	MString meshName;
	stat = triangulateMe(meshFn, numWantedVertices, numVerticesInMeshCreated, subsampleBoundary, meshName);
	MCHECKERROR(stat, "Failed to triangulate the given polygon");
	
	//for some reason Maya doesn't automatically flush the output streams after a command is issued, so I do it manually
	cout.flush();
	cerr.flush();

	setResult(meshName);

	return MS::kSuccess;
}


//this function assumes that the mesh has only one polygon
MStatus TriangulatePolygonCmd::triangulateMe(MFnMesh& meshFn, int numWantedVertices, int& numVerticesInMeshCreated, bool subsampleBoundary, MString& meshName)
{
	MStatus stat = MS::kSuccess;


	MItMeshPolygon poly(meshFn.object(&stat));

	MVector normal;
	poly.getNormal(normal);

	const MVector zAxis(0.0, 0.0, 1.0);

	MQuaternion quaternion(zAxis, normal);

	MMatrix rotMat = quaternion.asMatrix();

	double polygonArea = 0.0;
	stat = poly.getArea(polygonArea);

	MIntArray vertexList;
	meshFn.getPolygonVertices(0, vertexList);

	int numVertices = vertexList.length();

	if(numVertices < 3)
	{
		return MS::kFailure;
	}

	std::vector<Complex> polygon(numVertices);
	MPointArray vertexArray;
	stat = meshFn.getPoints(vertexArray);

	for(int i = 0; i < numVertices; i++)
	{
		MPoint p = vertexArray[vertexList[i]];
		p = rotMat*p;
		polygon[i] = Complex(p.x, p.y);
	}


	std::vector<Complex> meshVertices;
	std::vector<std::pair<int, Complex> > boundaryVertices;
	std::vector<unsigned int> triangleIndices;

	bool result = triangulatePolygon(polygon, meshVertices, triangleIndices, boundaryVertices, polygonArea / (double)numWantedVertices, subsampleBoundary);
	
	if(!result)
	{
		return MS::kFailure;
	}
	int numIndices = triangleIndices.size();
	int numTriangles = numIndices/3;
	numVerticesInMeshCreated = meshVertices.size();
	MPointArray points(numVerticesInMeshCreated);

	MMatrix invRotMat = rotMat.transpose();

	for(int i = 0; i < numVerticesInMeshCreated; i++)
	{
		MPoint p(meshVertices[i].real(), meshVertices[i].imag(), 0.0);
		points[i] = invRotMat*p;
	}

	MIntArray polygonCounts(numTriangles, 3); //[3, 3, 3, 3, 3, ... ,3]
	MIntArray polygonConnects(numIndices);
	for(int i = 0; i < numIndices; i++)
	{
		polygonConnects[i] = triangleIndices[i];
	}

	MObject meshObj = meshFn.object(&stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	MFnDagNode dagNode(meshObj, &stat);

	if(dagNode.parentCount(&stat) != 1 || stat != MS::kSuccess)
	{
		return MS::kFailure;
	}
	
	meshFn.create(numVerticesInMeshCreated, numTriangles, points, polygonCounts, polygonConnects, MObject::kNullObj, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

 	int numV = meshFn.numVertices();

	MFloatArray uArray, vArray;

	uArray.setLength(numV);
	vArray.setLength(numV);

	MItMeshVertex vertexIt(meshFn.object(&stat), &stat);
	for(; !vertexIt.isDone(); vertexIt.next())
	{
		MPoint p = vertexIt.position();
		int i = vertexIt.index();
		uArray[i] = (float)p[0];
		vArray[i] = (float)p[1];
	}

	stat = meshFn.setUVs(uArray, vArray, NULL); //NULL means use the current uvset
	MCHECKERROR(stat, "Can't set new values for uv's");

	MIntArray vertexCount;
	MIntArray uvIndices;
	stat = meshFn.getVertices(vertexCount, uvIndices);
	MCHECKERROR(stat, "Can't access Maya's mesh data");
	stat = meshFn.assignUVs(vertexCount, uvIndices, NULL);
	MCHECKERROR(stat, "Can't assign new values for uv's");

	meshName = meshFn.name(&stat);

	return MS::kSuccess;
}




bool TriangulatePolygonCmd::isUndoable() const
{
	return false;
}



MSyntax TriangulatePolygonCmd::syntax()
{
	MStatus stat = MS::kSuccess;
	MSyntax commandSyntax;

	stat = commandSyntax.addFlag(kSubsampleBoundaryShort, kSubsampleBoundaryLong, MSyntax::kBoolean);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	stat = commandSyntax.addFlag(kNumVerticesShort, kNumVerticesLong, MSyntax::kLong);
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1); //expect exactly one object
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");

	commandSyntax.useSelectionAsDefault(true);

	return commandSyntax;
}

