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
		if (sum > 0.0) {	x /= sum; y /= sum;	z /= sum;	}
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
	vector<set<size_t>> Adj;
	map<size_t, Coor> newVertex;
	FILE *fCache;

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
		for (auto it = Vertex.begin(); it != Vertex.end();)
		{
			newVertex[it->second - 1] = it->first;
			it = Vertex.erase(it);
		}
		fclose(fCache);
		Faces.clear();
	}

	void remove()
	{
		VertexNormal.clear();
		Faces.clear();
		Adj.clear();
		newVertex.clear();
	}
};

template<typename T>
void Swap(T& a, T&b)
{
	T tmp = a;
	a = b;
	b = tmp;
}