#include "MarchingCube.h"

long double GetValueFromTrilinearInterpolation(Coor V, long double V0[])
{
/*	Vxyz =	V000(1 - x) (1 - y) (1 - z) +
			V100 x(1 - y) (1 - z) +
			V010(1 - x) y(1 - z) +
			V001(1 - x) (1 - y) z +
			V101 x(1 - y) z +
			V011(1 - x) y z +
			V110 x y(1 - z) +
			V111 x y z		*/
	
	return 	(V0[0] * (1 - V.x) * (1 - V.y) * (1 - V.z)) +
			(V0[1] * V.x * (1 - V.y) * (1 - V.z)) +
			(V0[2] * (1 - V.x) * V.y *(1 - V.z)) +
			(V0[3] * (1 - V.x) * (1 - V.y) * V.z) +
			(V0[4] * V.x * (1 - V.y) * V.z) +
			(V0[5] * (1 - V.x) * V.y * V.z) +
			(V0[6] * V.x * V.y *(1 - V.z)) +
			(V0[7] * V.x * V.y * V.z);

}

void VertexInterp(size_t Threshold, Coor p1, Coor p2, double pv1, double pv2, size_t resultIndex, vector<Coor>& result)
{
	if (fabs(Threshold - pv1) < 0.00001) result[resultIndex] = p1;
	if (fabs(Threshold - pv2) < 0.00001) result[resultIndex] = p2;
	if (fabs(pv1 - pv2) < 0.00001) result[resultIndex] = p1;

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

void Polygonise(vector<long double>& PointValues, double x1, double y1, double z1, double x2, double y2, double z2, size_t Threshold, unordered_map<Coor, size_t, Coor_Hash_Func>& verts, Mesh& mesh)
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
			double ZDist, double YDist, double XDist,
			long long int Threshold,
			Mesh& mesh
		)
{
	size_t xWidth = (XSize);
	size_t xyWidth = xWidth * (YSize);

	double ZCoor = 0.0, YCoor, XCoor;
	unordered_map<Coor, size_t, Coor_Hash_Func> verts;

	XDist /= 2.0;
	YDist /= 2.0;
	ZDist /= 2.0;

	//#pragma omp parallel for
	for (long double z = 0.0; z < (long double) (ZSize - 1); z += 0.5)
	{
		YCoor = 0.0;
		for (long double y = 0.0; y < (long double) (YSize - 1); y += 0.5)
		{
			XCoor = 0.0;
			for (long double x = 0.0; x < (long double) (XSize - 1); x += 0.5)
			{
				//vector<long double> PointValue	{
													//Pixel_Array[z * xyWidth + y * xWidth + x], 
													//Pixel_Array[z * xyWidth + (y + 1) * xWidth + x],
													//Pixel_Array[z * xyWidth + (y + 1) * xWidth + (x + 1)],
													//Pixel_Array[z * xyWidth + y * xWidth + (x + 1)],
													//Pixel_Array[(z + 1) * xyWidth + y * xWidth + x],
													//Pixel_Array[(z + 1) * xyWidth + (y + 1) * xWidth + x],
													//Pixel_Array[(z + 1) * xyWidth + (y + 1) * xWidth + x + 1],
													//Pixel_Array[(z + 1) * xyWidth + y * xWidth, (x + 1)]
												//};
				long double PointValue[] =	{
												Pixel_Array[(int) z * xyWidth + (int) (y + 0.5) * xWidth + (int) x],
												Pixel_Array[(int) z * xyWidth + (int) (y + 0.5) * xWidth + (int) (x + 0.5)],
												Pixel_Array[(int) (z + 0.5) * xyWidth + (int) (y + 0.5) * xWidth + (int) x],
												Pixel_Array[(int) z * xyWidth + (int) y * xWidth + (int) x],
												Pixel_Array[(int) z * xyWidth + (int) y * xWidth + (int) (x + 0.5)],
												Pixel_Array[(int) z * xyWidth + (int) (y + 0.5) * xWidth + (int) x],
												Pixel_Array[(int) (z + 0.5) * xyWidth + (int) (y + 0.5) * xWidth + (int) (x + 0.5)],
												Pixel_Array[(int) (z + 0.5) * xyWidth + (int) y * xWidth, (int) (x + 0.5)]
											};

				vector<long double> PointValues	{
													GetValueFromTrilinearInterpolation({x, y, z}, PointValue),
													GetValueFromTrilinearInterpolation({z , y + 0.5, x}, PointValue),
													GetValueFromTrilinearInterpolation({x + 0.5, y + 0.5, z}, PointValue),
													GetValueFromTrilinearInterpolation({x + 0.5, y, z}, PointValue),
													GetValueFromTrilinearInterpolation({x, y, z + 0.5}, PointValue),
													GetValueFromTrilinearInterpolation({x, y + 0.5, z + 0.5}, PointValue),
													GetValueFromTrilinearInterpolation({x + 0.5, y + 0.5, z + 0.5}, PointValue),
													GetValueFromTrilinearInterpolation({x + 0.5, y, z + 0.5}, PointValue)
												};

				Polygonise(PointValues, XCoor, YCoor, ZCoor, XCoor + XDist, YCoor + YDist, ZCoor + ZDist, Threshold, verts, mesh);

				XCoor += XDist;
				PointValues.erase(PointValues.begin(), PointValues.end());
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
