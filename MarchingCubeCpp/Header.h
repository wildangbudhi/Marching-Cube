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

	long double Longitude()
	{
		return sqrtl((x*x) + (y*y) + (z*z));
	}

	void normalization()
	{
		long double sum = this->Longitude();
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

	void done()
	{
		for (auto it = Vertex.begin(); it != Vertex.end();)
		{
			newVertex[it->second - 1] = it->first;
			it = Vertex.erase(it);
		}
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