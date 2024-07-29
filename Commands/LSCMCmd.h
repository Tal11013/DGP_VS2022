#pragma once
#include "stdafx.h"

class LSCMCmd : public MPxCommand {
public:
    LSCMCmd();
    MStatus doIt(const MArgList& argList) override;
    static void* creator();
    static MSyntax syntax();
    static MString commandName();
    bool isUndoable() const override;
};