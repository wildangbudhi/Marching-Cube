#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <errno.h>
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
		return hash<long double>()(p.x) ^ hash<long double>()(p.y) ^ hash <long double>()(p.z);
	}
};

struct Mesh
{
	unordered_map<Coor, size_t, Coor_Hash_Func> Vertex;
	vector<pair<Coor, size_t> > VertexNormal;
	vector<vector<size_t> > Faces;
	FILE *vCache, *vnCache, *fCache;

	size_t push_Vertex(const Coor key, const size_t val)
	{
		errno = 0;

		if (!Vertex.size())
		{
			vCache = fopen("vCache", "w");
			fprintf(vCache, "\n# List of geometric vertices\n");
		}
		fprintf(vCache, "v %.3Lf %.3Lf %.3Lf\n", key.x, key.y, key.z);
		Vertex[key] = val;

		return (size_t)errno;
	}

	size_t insert_VertexNormal(const size_t idx, const Coor N, const size_t size)
	{
		errno = 0;

		if (VertexNormal.size() <= 3)
		{
			vnCache = fopen("vnCache", "w");
			fprintf(vnCache, "\n# List of vertex normals\n");
		}
		fprintf(vnCache, "vn %.3Lf %.3Lf %.3Lf\n", N.x, N.y, N.z);
		VertexNormal[idx] = { N, size };

		return (size_t)errno;
	}

	size_t push_Faces(vector<size_t> f)
	{
		errno = 0;

		if (!Faces.size())
		{
			fCache = fopen("fCache", "w");
			fprintf(fCache, "\n# Polygonal face element\n");
		}
		fprintf(fCache, "f %lld//%lld %lld//%lld %lld//%lld\n", f[0], f[0], f[1], f[1], f[2], f[2]);
		Faces.push_back(f);

		return (size_t)errno;
	}

	void done()
	{
		fclose(vCache);
		fclose(vnCache);
		fclose(fCache);
		Vertex.clear();
		VertexNormal.clear();
		Faces.clear();
	}
};