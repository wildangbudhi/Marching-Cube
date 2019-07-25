#include "Header.h"

void VertexInterp(size_t Threshold, Coor p1, Coor p2, double pv1, double pv2, size_t resultIndex, vector<Coor> &result);
void CalculateNormal(Mesh &mesh);
void Polygonise(vector<double>& PointValues, double x1, double y1, double z1, double x2, double y2, double z2, size_t Threshold, unordered_map<Coor, size_t, Coor_Hash_Func>& verts, Mesh& mesh);
void March(const double * Pixel_Array,
	const size_t ZSize,
	const size_t YSize,
	const size_t XSize,
	const double ZDist,
	const double YDist,
	const double XDist,
	long long int Threshold,
	Mesh& mesh
);
void MakeOBJ(const char *name, Mesh& mesh);
