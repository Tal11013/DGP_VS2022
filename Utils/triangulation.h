#pragma once

#include <vector>
#include <complex>


void Exactinit();
void freeTriangleStructure(struct triangulateio& triangleStruct);
double CounterClockWise(const double pa[2], const double pb[2], const double pc[2]);
double CounterClockWise(const std::complex<double>& pa, const std::complex<double>& pb, const std::complex<double>& pc);
bool triangulatePolygon(const std::vector<std::complex<double> >& polygonPoints, std::vector<std::complex<double> >& meshVertices,
						std::vector<unsigned int>& triangleIndices, std::vector<std::pair<int, std::complex<double> > >& boundaryVertices,
						double maxTriangleArea, bool subsampleBoundaryEdges = true);
bool triangulatePolygonWithoutAddingVertices(const std::vector<std::complex<double> >& polygonPoints, std::vector<unsigned int>& triangleIndices);
