#pragma once
// #define GPU_PROC


#include "Utils/STL_Macros.h"
#include "Utils/GMM_Macros.h"

class SpaceDeformer2D : public MPxDeformerNode
{
public:
	SpaceDeformer2D();
	virtual ~SpaceDeformer2D();

	static void* creator();
	static MStatus initialize();

	virtual MStatus deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex);
	MStatus updateCage(MFnMesh& cageMeshFn);
	MStatus doSetup(MItGeometry& iter);


public:
	const static MTypeId mTypeId;
	const static MString mTypeName;

protected:
	static MObject mCageAttr;
	static MObject mCoordinateTypeAttr;

protected:

	bool mIsFirstTime;

	GMMDenseComplexColMatrix mCageVertices; //this matrix is actually a column vector. dimensions are: n x 1
	GMMDenseComplexColMatrix mShiftedCageVertices; //this matrix is actually a column vector. dimensions are: n x 1
	GMMDenseComplexColMatrix mCoordinates; //this matrix has dimensions m x n
	GMMDenseComplexColMatrix mInterpolCoordinates; // this matrix has dimensions m x n
	GMMDenseComplexColMatrix mInternalPoints; //this matrix is actually a column vector. dimensions are: m x 1

};
