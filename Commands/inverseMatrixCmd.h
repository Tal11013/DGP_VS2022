#pragma once
#include "stdafx.h"

class inverseMatrixCmd : public MPxCommand {
public:
    inverseMatrixCmd();
    virtual MStatus doIt(const MArgList& argList);
    static void* creator();
    static MSyntax syntax();
    static MString commandName();
    virtual bool isUndoable() const;
};
