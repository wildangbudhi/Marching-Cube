#include "Header.h"

void March(	const double * Pixel_Array,
			const size_t ZSize, const size_t YSize, const size_t XSize,
			double ZDist, double YDist, double XDist,
			long long int Threshold,
			Mesh& mesh
		);
void Smoothing(Mesh & mesh, size_t rounds);
size_t MakeOBJ(const char *name, Mesh & mesh, bool isDoubleSided);