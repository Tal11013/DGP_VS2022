#pragma once



//////  MAYA API  //////
#include <maya/MStreamUtils.h>
#include <maya/MSelectionList.h>
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MArgList.h>
#include <maya/MPxCommand.h>
#include <maya/MPxHwShaderNode.h>
#include <maya/MMessage.h>
#include <maya/MNodeMessage.h>
#include <maya/MPointArray.h>
#include <maya/MPoint.h>
#include <maya/MBoundingBox.h>
#include <maya/MGlobal.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MImage.h>
#include <maya/MPlugArray.h>
#include <maya/MFnStringData.h>
#include <maya/MEulerRotation.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MMatrix.h>
#include <maya/MPxSurfaceShape.h>
#include <maya/MPxSurfaceShapeUI.h>
#include <maya/MTextureEditorDrawInfo.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItMeshFaceVertex.h>
#include <maya/MItMeshEdge.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFloatPointArray.h>
#include <maya/MPxDeformerNode.h>
#include <maya/MItGeometry.h>
#include <maya/MDagModifier.h>
#include <maya/MQuaternion.h>


//////  STL  //////
#include <complex>
#include <vector>

//////  Runtime  //////
#include <ctime>
#include <cassert>
#include <limits>


//////  GMM  //////
#include "gmm/gmm.h"
