#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "MarchingCube.h"
//#include "Laplacian_Smoothing.h"

namespace py = pybind11;

void ConvertHu(py::array_t<double> pixel_array, const double slope, const double intercept)
{
	auto rA = pixel_array.request();

	double *Voxel = (double *)rA.ptr;

#pragma omp parallel for
	for (int i = 0; i < (rA.shape[0] * rA.shape[1] * rA.shape[2]); i++)
		Voxel[i] = (Voxel[i] * slope) + intercept;
}

void
MarchingCube(py::array_t<double> pixel_array,
	const size_t ZSize, const size_t YSize, const size_t XSize,
	const double ZDist, const double YDist, const double XDist,
	const long long int trs,
	const char * name
)
{
	auto rA = pixel_array.request();
	/*double *Pixel_Array = (double *) rA.ptr;

	int x = 0, y = 0, z = 0;
	size_t xWidth = (XSize);
	size_t xyWidth = xWidth * (YSize); */

	
	Mesh mesh;

	cout << "March Start" << endl;
	March((double *)rA.ptr, ZSize, YSize, XSize, ZDist, YDist, XDist, trs, mesh);
	/*cout << "Smoothing Start" << endl;
	Smoothing(mesh, 2);*/
	//cout << "Make OBJ Start " << name << endl;
	cout << MakeOBJ(name, mesh);
}

void
Coba(const char* name)
{
	printf("hahahh");
}

PYBIND11_MODULE(MarchingCubeCpp, m)
{
	/*py::class_<Coor>(m, "Coor")
		.def_readwrite("x", &Coor::x)
		.def_readwrite("y", &Coor::y)
		.def_readwrite("z", &Coor::z);*/

	m.def("MarchingCube", &MarchingCube, "Making 3D Object from Series of Pixel Array or Voxel");
	m.def("ConvertHu", &ConvertHu, "Convert Pixel Array from Gray Scale to Hu Scale");
	m.def("Coba", &Coba, "Function for Experiment");
}