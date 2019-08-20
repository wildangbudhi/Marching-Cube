#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "MarchingCube.h"

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
	const size_t trs,
	const size_t roundFactor,
	const char * name,
	bool isDoubleSided = false
)
{
	auto rA = pixel_array.request();
	Mesh mesh;

	cout << "March Start" << endl;
	March((double *)rA.ptr, ZSize, YSize, XSize, ZDist, YDist, XDist, trs, mesh);

	Smoothing(mesh, roundFactor);
	MakeOBJ(name, mesh, isDoubleSided);
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

	m.def(	"MarchingCube", &MarchingCube, 
			py::arg("pixel_array"), 
			py::arg("ZSize"), py::arg("YSize"), py::arg("XSize"),
			py::arg("ZDist"), py::arg("YDist"), py::arg("XDist"),
			py::arg("Threshold"),
			py::arg("roundFactor"),
			py::arg("Name"),
			py::arg("isDoubleSided") = false,
			"Making 3D Object from Series of Pixel Array or Voxel");
	m.def("ConvertHu", &ConvertHu, "Convert Pixel Array from Gray Scale to Hu Scale");
	m.def("Coba", &Coba, "Function for Experiment");
}