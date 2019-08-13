#include "Header.h"

void VertexInterp(size_t Threshold, Coor p1, Coor p2, long double pv1, long double pv2, size_t resultIndex, vector<Coor> &result);
void Polygonise(vector<long double>& PointValues, double x1, double y1, double z1, double x2, double y2, double z2, size_t Threshold, Mesh& mesh);
void March(const double * Pixel_Array,
	const size_t ZSize,
	const size_t YSize,
	const size_t XSize,
	double ZDist,
	double YDist,
	double XDist,
	long long int Threshold,
	Mesh& mesh
);
size_t MakeOBJ(const char *name, Mesh & mesh);
