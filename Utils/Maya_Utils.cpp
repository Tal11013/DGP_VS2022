#include <maya/MStatus.h>
#include <maya/MFnDagNode.h>
#include <maya/MDagPath.h>
#include "triangulation.h"
#include <complex>
#include <vector>

#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MFnMesh.h>

#include "Utils/Maya_Utils.h"



//given an object, the function checks if it actually represents a mesh. If it is not a mesh but it is a dag node,
//the function will try to retrieve a shape that is attached to this dag node and then checks if this shape is a mesh
MStatus Maya_Utils::getMe_a_Mesh(const MObject& inObject, MObject& meshObject)
{
	MStatus stat = MS::kSuccess;

	if(inObject.hasFn(MFn::kMesh))
	{
		meshObject = inObject;
		return MS::kSuccess;
	}


	if(inObject.hasFn(MFn::kDagNode))
	{
		MFnDagNode dagNode(inObject, &stat);
		if(stat != MS::kSuccess) return stat;

		MDagPath dagPath;
		stat = dagNode.getPath(dagPath);
		if(stat != MS::kSuccess) return stat;

		stat = dagPath.extendToShape();
		if(stat != MS::kSuccess) return stat;

		bool isMesh = dagPath.hasFn(MFn::kMesh, &stat);

		if(isMesh)
		{
			meshObject = dagPath.node(&stat);
			return stat;
		}
	}

	return MS::kFailure; //not a mesh
}



//checks if a given mesh has the specified uv set name
bool Maya_Utils::hasThatUVSetName(const MFnMesh& meshFn, const MString& uvSetName)
{
	MStringArray allUVSetNames;
	meshFn.getUVSetNames(allUVSetNames);

	bool hasThatName = false;

	for(int i = 0; i < allUVSetNames.length(); i++)
	{
		if(allUVSetNames[i] == uvSetName)
		{
			hasThatName = true;
			break;
		}
	}
	return hasThatName;
}