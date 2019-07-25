#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <iostream>
#include <omp.h>
#include <pybind11/stl.h>
#include "Table.h"

using namespace std;

struct Coor
{
	double x, y, z;
	bool operator==(const Coor& p) const { return (x == p.x) && (y == p.y) && (z == p.z); }
};

struct Coor_Hash_Func
{
	size_t operator() (const Coor& p) const
	{
		return hash<double>()(p.x) ^ hash<double>()(p.y) ^ hash<double>()(p.z);
	}
};

struct Mesh
{
	vector<Coor> Vertex;
	vector<Coor> VertexNormal;
	vector<vector<size_t> > Faces;
};