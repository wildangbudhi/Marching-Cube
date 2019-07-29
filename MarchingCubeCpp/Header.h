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
	long double x, y, z;
	bool operator==(const Coor& p) const { return (x == p.x) && (y == p.y) && (z == p.z); }
	void normalization() 
	{
		long double sum = sqrtl((x*x) + (y*y) + (z*z));
		x /= sum; 
		y /= sum; 
		z /= sum;
	}
};

struct Coor_Hash_Func
{
	size_t operator() (const Coor& p) const
	{
		return hash<long double>()(p.x) ^ hash<long double>()(p.y) ^ hash <long double> ()(p.z);
	}
};

struct Mesh
{
	vector<Coor> Vertex;
	vector<pair<Coor, size_t> > VertexNormal;
	vector<vector<size_t> > Faces;
};