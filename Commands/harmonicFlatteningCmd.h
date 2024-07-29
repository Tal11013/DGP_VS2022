#pragma once
#include "stdafx.h"

class harmonicFlatteningCmd : public MPxCommand {
public:
    harmonicFlatteningCmd();
    MStatus doIt(const MArgList& argList) override;
    static void* creator();
    static MSyntax syntax();
    static MString commandName();
    bool isUndoable() const override;
};