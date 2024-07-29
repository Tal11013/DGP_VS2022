#include "stdafx.h"

#include "SpaceDeformer2D.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/MatlabInterface.h"

const MTypeId SpaceDeformer2D::mTypeId(0x6723c);
const MString SpaceDeformer2D::mTypeName("SpaceDeformer2D");

MObject SpaceDeformer2D::mCageAttr;
MObject SpaceDeformer2D::mCoordinateTypeAttr;



SpaceDeformer2D::SpaceDeformer2D() : mIsFirstTime(true) {}

SpaceDeformer2D::~SpaceDeformer2D() = default;

void* SpaceDeformer2D::creator()
{
	return new SpaceDeformer2D();
}

MStatus SpaceDeformer2D::initialize()
{
	MStatus stat;

	MFnTypedAttribute cageAttr;
	mCageAttr = cageAttr.create("cage" ,"cage", MFnData::kMesh, MObject::kNullObj, &stat);
	CHECK_MSTATUS(addAttribute(mCageAttr));
	CHECK_MSTATUS(attributeAffects(mCageAttr, outputGeom));

	MFnEnumAttribute coordinateTypeAttr;
	mCoordinateTypeAttr = coordinateTypeAttr.create("coordinateType" ,"coordinateType", 0, &stat);
	CHECK_MSTATUS(coordinateTypeAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mCoordinateTypeAttr));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Cauchy", 0));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Cauchy Interpolation", 1));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Point to point", 2));
	CHECK_MSTATUS(attributeAffects(mCoordinateTypeAttr, outputGeom));

	return MStatus::kSuccess;
}


MStatus SpaceDeformer2D::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex)
{
	MStatus stat;

	MDataHandle coordHandle = block.inputValue(mCoordinateTypeAttr, &stat);
	short coordinateType = coordHandle.asShort();

	MDataHandle handle = block.inputValue(mCageAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	MObject cageMesh = handle.asMesh();

	MFnMesh cageMeshFn(cageMesh, &stat);
	if(stat != MS::kSuccess) {
		return stat;
	}

	updateCage(cageMeshFn);

	if(mIsFirstTime) {
		stat = doSetup(iter);
		CHECK_MSTATUS_AND_RETURN_IT(stat);
		mIsFirstTime = false;
	}

	///// add your code here //////
	///////////////////////////////
	//compute the deformation of all the internal points.
	//This is done by simply multiplying the coordinate matrix by the cage vertices vector
	#ifdef GPU_PROC
		// Multiply matrices in the GPU
		MatlabGMMDataExchange::SetEngineDenseMatrix("cage", mCageVertices);
		MatlabInterface::GetEngine().Eval("gpu_cage = gpuArray(cage)");
		if (coordinateType == 0 ) {
			MatlabInterface::GetEngine().Eval("gpu_internal = gpu_coords * gpu_cage");
		} else if (coordinateType == 1) {
			MatlabInterface::GetEngine().Eval("gpu_internal = gpu_interpol_coords * gpu_cage");
		}
		MatlabInterface::GetEngine().Eval("internal = gather(gpu_internal)");
		MatlabGMMDataExchange::GetEngineDenseMatrix("internal", mInternalPoints);
	#else
		if (coordinateType == 0) {
			mult(mCoordinates, mCageVertices, mInternalPoints);	
		} else if (coordinateType == 1) {
			mult(mInterpolCoordinates, mCageVertices, mInternalPoints);	
		}
		
	#endif


	///////////////////////////////
	///////////////////////////////


	//update the new deformed position of all the internal vertices
	for(iter.reset(); !iter.isDone(); iter.next())
	{
		int i = iter.index();


		///// add your code here //////
		///////////////////////////////

		//update c to be the deformed position of the i'th vertex

		///////////////////////////////
		///////////////////////////////


		Complex c = mInternalPoints(i, 0);

		iter.setPosition(MPoint(c.real(), c.imag(), 0.0));
	}

	return stat;
}


MStatus SpaceDeformer2D::updateCage(MFnMesh& cageMeshFn)
{
	MStatus stat;

	int numFaces = cageMeshFn.numPolygons(&stat);

	assert(numFaces == 1);

	MIntArray vertexIndices;
	cageMeshFn.getPolygonVertices(0, vertexIndices);
	int numV = vertexIndices.length();
	assert(numV >= 3);

	clear(mCageVertices);
	resize(mCageVertices, numV, 1);

	MPointArray vertexArray;
	stat = cageMeshFn.getPoints(vertexArray);

	assert(numV == vertexArray.length());

	for(int i = 0; i < numV; i++)
	{
		MPoint p = vertexArray[i];
		Complex c(p[0], p[1]);
		mCageVertices(i, 0) = c;
	}
	return MS::kSuccess;
}


MStatus SpaceDeformer2D::doSetup(MItGeometry& iter)
{
	MStatus stat;

	int m = iter.count(&stat); //num internal points
	int n = mCageVertices.nrows(); //num vertices in the cage

	clear(mCoordinates);
	resize(mCoordinates, m, n);

	clear(mInternalPoints);
	resize(mInternalPoints, m, 1);

	clear(mShiftedCageVertices);
	resize(mShiftedCageVertices, n, 1);

	clear(mInterpolCoordinates);
	resize(mInterpolCoordinates, m, n);

	const Complex imag_unit = {0, 1};

	// First, create the shifted mesh
	for (int j = 0; j < n; ++j) {
		constexpr double d = 0.00001;
		Complex z0 = mCageVertices((j - 1 + n) % n, 0);
		Complex z1 = mCageVertices( j, 0);
		Complex z2 = mCageVertices((j + 1) % n, 0);

		// Shift z1 along the angular bisector of the angle (z0-z1-z2)
		Complex u = z0 - z1, v = z2 - z1;
		Complex bisector = abs(u) * v + abs(v) * u;
		Complex normalized_bisector = - bisector / abs(bisector) * d;

		// Compute the direction to shift
		const double alpha = fmod(arg(u) - arg(v) + 2*M_PI, 2*M_PI) ;
		if (alpha >= M_PI) {
			normalized_bisector *= -1;
		}

		mShiftedCageVertices(j, 0) = z1 + normalized_bisector;
	}

	for(iter.reset(); !iter.isDone(); iter.next()) {


		int i = iter.index();
		MPoint pt = iter.position();

		Complex z0(pt[0], pt[1]); //internal point

		mInternalPoints(i, 0) = z0;

		for(int j = 0; j < n; j++) {
			///// add your code here //////
			///////////////////////////////

			//update K to be value of the j'th coordinate at the i'th internal point
			// j'th, (j-1)'th and (j+1)'th cage points
			Complex zj0 = mShiftedCageVertices((n+j-1)%n, 0);
			Complex zj = mShiftedCageVertices(j, 0);
			Complex zj1 = mShiftedCageVertices((j+1) % n, 0);

			// A_j's and B_j's in the formula
			Complex Bj0 = zj0 - z0, Bj = zj - z0, Bj1 = zj1 - z0,
					Aj = zj - zj0, Aj1 = zj1 - zj;

			const Complex K = 1.0 / (2 * M_PI * imag_unit) * (Bj1 / Aj1 * log(Bj1 / Bj) -
				Bj0 / Aj * log(Bj / Bj0));

			///////////////////////////////
			///////////////////////////////
			mCoordinates(i, j) = K;
		}
	}

	// The intermediate matrix. size n x n
	GMMDenseComplexColMatrix interpolMidMatrix(n, n);


	// For part 2, calculate the intermediate matrix
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			Complex z0 = mCageVertices(i, 0);

			Complex zj0 = mShiftedCageVertices((n+j-1)%n, 0);
			Complex zj = mShiftedCageVertices(j, 0);
			Complex zj1 = mShiftedCageVertices((j+1) % n, 0);

			// A_j's and B_j's in the formula
			Complex Bj0 = zj0 - z0, Bj = zj - z0, Bj1 = zj1 - z0,
					Aj = zj - zj0, Aj1 = zj1 - zj;

			const Complex K = 1.0 / (2 * M_PI * imag_unit) * (Bj1 / Aj1 * log(Bj1 / Bj) -
				Bj0 / Aj * log(Bj / Bj0));

			interpolMidMatrix(i, j) = K;
		}
	}

	// Calculate the interpolating coordinates
	MatlabGMMDataExchange::SetEngineDenseMatrix("coords", mCoordinates);
	MatlabGMMDataExchange::SetEngineDenseMatrix("intermediate", interpolMidMatrix);
	MatlabInterface::GetEngine().EvalToCout("interpol_coords = coords / intermediate");
	MatlabGMMDataExchange::GetEngineDenseMatrix("interpol_coords", mInterpolCoordinates);

	#ifdef GPU_PROC
		// Transfer matrices to MatLab and set them as GPUArrays
		MatlabInterface::GetEngine().Eval("gpu_coords = gpuArray(coords)");
		MatlabInterface::GetEngine().Eval("gpu_interpol_coords = gpuArray(interpol_coords)");
	#endif


	return MS::kSuccess;
}
