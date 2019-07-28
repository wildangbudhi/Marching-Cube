#include "MarchingCube.h"

void VertexInterp(size_t Threshold, Coor p1, Coor p2, double pv1, double pv2, size_t resultIndex, vector<Coor>& result)
{
	if (abs(Threshold - pv1) < 0.00001) result[resultIndex] = p1;
	if (abs(Threshold - pv2) < 0.00001) result[resultIndex] = p2;
	if (abs(pv1 - pv2) < 0.00001) result[resultIndex] = p1;

	double mu = (Threshold - pv1) / (pv1 - pv2);

	result[resultIndex] = { (p1.x + mu * (p2.x - p1.x)), 
							(p1.y + mu * (p2.y - p1.y)), 
							(p1.z + mu * (p2.z - p1.z))	};
}

void CalculateNormal(Mesh & mesh)
{
	for (size_t i = 0; i < mesh.Vertex.size(); i++)
		mesh.VertexNormal.push_back({ { 0,0,1 } , 0});

	#pragma omp parallel for
	for (long long int i = 0; i < mesh.Faces.size(); i++)
	{

		// a = p2 - p1; b = p3 - p1
		Coor a = {	mesh.Vertex[mesh.Faces[i][1]].x - mesh.Vertex[mesh.Faces[i][0]].x,
					mesh.Vertex[mesh.Faces[i][1]].y - mesh.Vertex[mesh.Faces[i][0]].y,
					mesh.Vertex[mesh.Faces[i][1]].z - mesh.Vertex[mesh.Faces[i][0]].z };

		Coor b = {	mesh.Vertex[mesh.Faces[i][2]].x - mesh.Vertex[mesh.Faces[i][0]].x,
					mesh.Vertex[mesh.Faces[i][2]].y - mesh.Vertex[mesh.Faces[i][0]].y,
					mesh.Vertex[mesh.Faces[i][2]].z - mesh.Vertex[mesh.Faces[i][0]].z };

		// Cross Product
		Coor N = { (a.y * b.z) - (a.z * b.y),
					(a.z * b.x) - (a.x * b.z),
					(a.x * b.y) - (a.x * b.y) };

		// Normalization
		N.normalization();

		// Combine Normal
		for (int j = 0; j < 3; j++)
		{
			size_t SizeTemp;
			if ((SizeTemp = mesh.VertexNormal[mesh.Faces[i][j]].second) == 0)
				mesh.VertexNormal[mesh.Faces[i][j]] = { N, SizeTemp + 1 };
			else
			{
				// Mean of Vertex Normal
				Coor temp = mesh.VertexNormal[mesh.Faces[i][j]].first;
				N = {	((SizeTemp * temp.x) + (N.x)) / (SizeTemp + 1)	,
						((SizeTemp * temp.y) + (N.y)) / (SizeTemp + 1)	,
						((SizeTemp * temp.z) + (N.z)) / (SizeTemp + 1)	,
					};

				// Normalization
				N.normalization();

				mesh.VertexNormal[mesh.Faces[i][j]] = { N, SizeTemp + 1 };
			}
		}
	}

}

void Polygonise(vector<double>& PointValues, double x1, double y1, double z1, double x2, double y2, double z2, size_t Threshold, unordered_map<Coor, size_t, Coor_Hash_Func>& verts, Mesh& mesh)
{
	size_t cubeIndex = 0;
	int shifter = 1;

	for (int i = 0; i < PointValues.size(); shifter <<= 1, i++)
		if (PointValues[i] < Threshold) cubeIndex |= shifter;

	if (EDGE_TABLE[cubeIndex] == 0) return;

	vector<Coor> vertlist(12);

	if (EDGE_TABLE[cubeIndex] & 1) VertexInterp(Threshold, { x1, y1, z1 }, { x1,y2,z1 }, PointValues[0], PointValues[1], 0, vertlist);
	if (EDGE_TABLE[cubeIndex] & 2) VertexInterp(Threshold, { x1,y2,z1 }, { x2,y2,z1 }, PointValues[1], PointValues[2], 1, vertlist);
	if (EDGE_TABLE[cubeIndex] & 4) VertexInterp(Threshold, { x2,y2,z1 }, { x2,y1,z1 }, PointValues[2], PointValues[3], 2, vertlist);
	if (EDGE_TABLE[cubeIndex] & 8) VertexInterp(Threshold, { x2,y1,z1 }, { x1,y1,z1 }, PointValues[3], PointValues[0], 3, vertlist);
	if (EDGE_TABLE[cubeIndex] & 16) VertexInterp(Threshold, { x1,y1,z2 }, { x1,y2,z2 }, PointValues[4], PointValues[5], 4, vertlist);
	if (EDGE_TABLE[cubeIndex] & 32) VertexInterp(Threshold, { x1,y2,z2 }, { x2,y2,z2 }, PointValues[5], PointValues[6], 5, vertlist);
	if (EDGE_TABLE[cubeIndex] & 64) VertexInterp(Threshold, { x2,y2,z2 }, { x2,y1,z2 }, PointValues[6], PointValues[7], 6, vertlist);
	if (EDGE_TABLE[cubeIndex] & 128) VertexInterp(Threshold, { x2,y1,z2 }, { x1,y1,z2 }, PointValues[7], PointValues[4], 7, vertlist);
	if (EDGE_TABLE[cubeIndex] & 256) VertexInterp(Threshold, { x1,y1,z1 }, { x1,y1,z2 }, PointValues[0], PointValues[4], 8, vertlist);
	if (EDGE_TABLE[cubeIndex] & 512) VertexInterp(Threshold, { x1,y2,z1 }, { x1,y2,z2 }, PointValues[1], PointValues[5], 9, vertlist);
	if (EDGE_TABLE[cubeIndex] & 1024) VertexInterp(Threshold, { x2,y2,z1 }, { x2,y2,z2 }, PointValues[2], PointValues[6], 10, vertlist);
	if (EDGE_TABLE[cubeIndex] & 2048) VertexInterp(Threshold, { x2,y1,z1 }, { x2,y1,z2 }, PointValues[3], PointValues[7], 11, vertlist);

	for (size_t i = 0; TRIANGLE_TABLE[cubeIndex][i] != -1; i += 3)
	{
		size_t v[3];
		
		for (int j = 0; j < 3; j++)
		{
			if (verts.count(vertlist[TRIANGLE_TABLE[cubeIndex][i + j]]))
			{
				v[j] = verts[vertlist[TRIANGLE_TABLE[cubeIndex][i + j]]];
			}
			else 
			{
				v[j] = verts.size() + 1;
				verts[ vertlist[ TRIANGLE_TABLE[cubeIndex][i + j] ] ] = v[j];
				mesh.Vertex.push_back( vertlist[ TRIANGLE_TABLE[cubeIndex][i + j] ] );
			}
		}

		mesh.Faces.push_back(vector<size_t>{ v[0], v[1], v[2] });
	}

}

void March(	const double * Pixel_Array, 
			const size_t ZSize, const size_t YSize, const size_t XSize, 
			const double ZDist, const double YDist, const double XDist,
			long long int Threshold,
			Mesh& mesh
		)
{
	size_t xWidth = (XSize);
	size_t xyWidth = xWidth * (YSize);

	double ZCoor = 0.0, YCoor, XCoor;
	unordered_map<Coor, size_t, Coor_Hash_Func> verts;

	//#pragma omp parallel for
	for (int z = 0; z < (ZSize - 1); z++)
	{
		YCoor = 0.0;
		for (int y = 0; y < (YSize - 1); y++)
		{
			XCoor = 0.0;
			for (int x = 0; x < (XSize - 1); x++)
			{
				vector<double> PointValue	{
												Pixel_Array[z * xyWidth + y * xWidth + x], 
												Pixel_Array[z * xyWidth + (y + 1) * xWidth + x],
												Pixel_Array[z * xyWidth + (y + 1) * xWidth + (x + 1)],
												Pixel_Array[z * xyWidth + y * xWidth + (x + 1)],
												Pixel_Array[(z + 1) * xyWidth + y * xWidth + x],
												Pixel_Array[(z + 1) * xyWidth + (y + 1) * xWidth + x],
												Pixel_Array[(z + 1) * xyWidth + (y + 1) * xWidth + x + 1],
												Pixel_Array[(z + 1) * xyWidth + y * xWidth, (x + 1)]
											};

				Polygonise(PointValue, XCoor, YCoor, ZCoor, XCoor + XDist, YCoor + YDist, ZCoor + ZDist, Threshold, verts, mesh);

				XCoor += XDist;
				PointValue.erase(PointValue.begin(), PointValue.end());
			}
			YCoor += YDist;
		}
		ZCoor += ZDist;
	}

	cout << "Start Calculate Normal " << endl;
	CalculateNormal(mesh);
}

void MakeOBJ(const char * name, Mesh & mesh)
{
	FILE *f1 = fopen(name, "w");

	fprintf(f1, "# OBJ from MarchigCube by 3D BAMAG\n");

	fprintf(f1, "\n# List of geometric vertices\n");
	for (size_t i = 0; i < mesh.Vertex.size(); i++)
		fprintf(f1, "v %lf %lf %lf\n", mesh.Vertex[i].x, mesh.Vertex[i].y, mesh.Vertex[i].z);

	fprintf(f1, "\n# List of vertex normals\n");
	for (size_t i = 0; i < mesh.VertexNormal.size(); i++) 
		fprintf(f1, "vn %Lf %Lf %Lf\n", mesh.VertexNormal[i].first.x, mesh.VertexNormal[i].first.y, mesh.VertexNormal[i].first.z);

	fprintf(f1, "\n# Polygonal face element\n");
	for (size_t i = 0; i < mesh.Faces.size(); i++) 
		fprintf(f1, "f %lld//%lld %lld//%lld %lld//%lld\n", mesh.Faces[i][0], mesh.Faces[i][0], mesh.Faces[i][1], mesh.Faces[i][1], mesh.Faces[i][2], mesh.Faces[i][2]);

	fclose(f1);
}
