#pragma once



class TriangulatePolygonCmd : public MPxCommand
{

public:

	TriangulatePolygonCmd();
	virtual MStatus	doIt(const MArgList& argList);
	static void* creator();
	static MSyntax syntax();
	static MString commandName();
	MStatus triangulateMe(MFnMesh& meshFn, int numWantedVertices, int& numVerticesInMeshCreated, bool subsampleBoundary, MString& meshName);
	virtual bool isUndoable() const;

};

