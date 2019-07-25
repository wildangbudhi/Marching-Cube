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
	bool operator==(const Coor& p) const {	return (x == p.x) && (y == p.y) && (z == p.z);	}
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
