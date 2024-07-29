#pragma once



namespace Maya_Utils
{
	typedef std::complex<double> Complex;


	MStatus getMe_a_Mesh(const MObject& inObject, MObject& meshObject);
	bool hasThatUVSetName(const MFnMesh& meshFn, const MString& uvSetName);
}


