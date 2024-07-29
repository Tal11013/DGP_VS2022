#pragma once

class colorMeshVerticesCmd : public MPxCommand
{

public:
	colorMeshVerticesCmd();							
	virtual MStatus	doIt(const MArgList& argList);
	static void* creator();							
	static MSyntax syntax();
	static MString commandName();					
	virtual bool isUndoable() const;

};
