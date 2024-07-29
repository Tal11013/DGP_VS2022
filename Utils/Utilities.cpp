#include "Utils\Utilities.h"

#include <ctime>
#include <limits>

#include "jet.h"


void swap(double& a, double& b)
{
	double temp = b;
	b = a;
	a = temp;
}

void swap(void** a, void** b)
{
	void* temp = *b;
	*b = *a;
	*a = temp;
}

int getRandomNumber()
{
	static bool firstTime = true;
	if(firstTime)
	{
		srand((unsigned)time(NULL));
		firstTime = false;
	}
	int randNumber = rand();
	return randNumber;
}


double getRandomNumber(double low, double high)
{
	static bool firstTime = true;
	#define NORMALIZATION_RAND_FACTOR (3.0518509475997192297128208258309e-5) //(1.0 / RAND_MAX)
	if(firstTime)
	{
		srand((unsigned)time(NULL));
		firstTime = false;
	}
	double randNumber = (high - low)*(NORMALIZATION_RAND_FACTOR*rand()) + low;
	return randNumber;
}

double computePolygonArea(const std::vector<std::complex<double> >& polygon)
{
	int numVertices = polygon.size();
	double doubleArea = 0.0;

	std::complex<double> center(0.0, 0.0);
	for(int i = 0; i < numVertices; i++)
	{
		center += polygon[i];
	}
	for(int i = 0; i < numVertices; i++)
	{
		std::complex<double> v1 = polygon[i] - center;
		std::complex<double> v2 = polygon[i < numVertices - 1 ? i+1 : 0] - center;

		doubleArea += (v1.real()*v2.imag() - v1.imag()*v2.real());
	}
	return 0.5*doubleArea;
}

void mapColor(double scalarValue, float& R, float& G, float& B, double minValue, double maxValue)
{
	if(minValue >= maxValue)
	{
		return;
	}
	double t = (scalarValue - minValue) / (maxValue - minValue);
	int t_int = (int)floor(t*JET_MAP_NCOLS);
	if(t_int > JET_MAP_NCOLS - 1)
	{
		t_int = JET_MAP_NCOLS - 1;
	}
	if(t_int < 0)
	{
		t_int = 0;
	}
	double* mapTable = &jet_65536[0][0];
	R = (float)(mapTable[3*t_int]);
	G = (float)(mapTable[3*t_int + 1]);
	B = (float)(mapTable[3*t_int + 2]);
}
