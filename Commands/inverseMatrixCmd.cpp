#include "inverseMatrixCmd.h"

#include "stdafx.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"


inverseMatrixCmd::inverseMatrixCmd() {

}

void* inverseMatrixCmd::creator()
{
    return new inverseMatrixCmd();
}

MSyntax inverseMatrixCmd::syntax() {
    MSyntax commandSyntax;

	for (int i = 0; i < 9; ++i) {
	    commandSyntax.addArg(MSyntax::kDouble);
    }

    return commandSyntax;
}


MString inverseMatrixCmd::commandName() {
    return "inverseMatrixCmd";
}

bool inverseMatrixCmd::isUndoable() const
{
    return false;
}

MStatus inverseMatrixCmd::doIt(const MArgList& argList)
{
	MSyntax commandSyntax = syntax();

	GMMDenseColMatrix M(3, 3);
    for (int i = 0; i < 9; ++i)
    {
	    M(i / 3, i % 3) = argList.asDouble(i);
    }

    int result = MatlabGMMDataExchange::SetEngineDenseMatrix("M", M);
    MatlabInterface::GetEngine().Eval("Q = inv(M)");

    // First check if the matrix is invertible
    MatlabInterface::GetEngine().Eval("d = [det(M)]");
    GMMDenseColMatrix d;
	result = MatlabGMMDataExchange::GetEngineDenseMatrix("d", d);

    if (d(0, 0) == 0.0) {
	    cerr << "Error: Matrix isn't invertible" << endl;
        return MS::kSuccess;
    }

    GMMDenseColMatrix Q;
    result = MatlabGMMDataExchange::GetEngineDenseMatrix("Q", Q);
    cout << "The inverse matrix is: " <<  Q << endl;

    return MS::kSuccess;
}
