#include "MarchingCube.h"

void VertexInterp(size_t Threshold, Coor p1, Coor p2, long double pv1, long double pv2, size_t resultIndex, vector<Coor>& result)
{
	/*if (fabs(Threshold - pv1) < 0.00001) result[resultIndex] = p1;
	if (fabs(Threshold - pv2) < 0.00001) result[resultIndex] = p2;
	if (fabs(pv1 - pv2) < 0.00001) result[resultIndex] = p1;

	double mu = (Threshold - pv1) / (pv1 - pv2);

	result[resultIndex] = { (p1.x + mu * (p2.x - p1.x)),
							(p1.y + mu * (p2.y - p1.y)),
							(p1.z + mu * (p2.z - p1.z)) };*/

	result[resultIndex] = { (p1.x + p2.x) / 2,
							(p1.y + p2.y) / 2,
							(p1.z + p2.z) / 2 };
}

void Polygonise(vector<long double>& PointValues, double x1, double y1, double z1, double x2, double y2, double z2, size_t Threshold, Mesh& mesh)
{
	size_t cubeIndex = 0;
	int shifter = 1;

	for (int i = 0; i < PointValues.size(); shifter <<= 1, i++)
		if (PointValues[i] < Threshold) cubeIndex |= shifter;

	if (EDGE_TABLE[cubeIndex] == 0) return;

	vector<Coor> vertlist(12);

	vector<Coor> Point{
		{x1, y1, z1}, // 0
		{x1, y2, z1}, // 1
		{x2, y2, z1}, // 2
		{x2, y1, z1}, // 3
		{x1, y1, z2}, // 4
		{x1, y2, z2}, // 5
		{x2, y2, z2}, // 6
		{x2, y1, z2}, // 7
	};

	if (EDGE_TABLE[cubeIndex] & 1)		VertexInterp(Threshold, Point[0], Point[1], PointValues[0], PointValues[1], 0, vertlist);
	if (EDGE_TABLE[cubeIndex] & 2)		VertexInterp(Threshold, Point[1], Point[2], PointValues[1], PointValues[2], 1, vertlist);
	if (EDGE_TABLE[cubeIndex] & 4)		VertexInterp(Threshold, Point[2], Point[3], PointValues[2], PointValues[3], 2, vertlist);
	if (EDGE_TABLE[cubeIndex] & 8)		VertexInterp(Threshold, Point[3], Point[0], PointValues[3], PointValues[0], 3, vertlist);
	if (EDGE_TABLE[cubeIndex] & 16)		VertexInterp(Threshold, Point[4], Point[5], PointValues[4], PointValues[5], 4, vertlist);
	if (EDGE_TABLE[cubeIndex] & 32)		VertexInterp(Threshold, Point[5], Point[6], PointValues[5], PointValues[6], 5, vertlist);
	if (EDGE_TABLE[cubeIndex] & 64)		VertexInterp(Threshold, Point[6], Point[7], PointValues[6], PointValues[7], 6, vertlist);
	if (EDGE_TABLE[cubeIndex] & 128)	VertexInterp(Threshold, Point[7], Point[4], PointValues[7], PointValues[4], 7, vertlist);
	if (EDGE_TABLE[cubeIndex] & 256)	VertexInterp(Threshold, Point[0], Point[4], PointValues[0], PointValues[4], 8, vertlist);
	if (EDGE_TABLE[cubeIndex] & 512)	VertexInterp(Threshold, Point[1], Point[5], PointValues[1], PointValues[5], 9, vertlist);
	if (EDGE_TABLE[cubeIndex] & 1024)	VertexInterp(Threshold, Point[2], Point[6], PointValues[2], PointValues[6], 10, vertlist);
	if (EDGE_TABLE[cubeIndex] & 2048)	VertexInterp(Threshold, Point[3], Point[7], PointValues[3], PointValues[7], 11, vertlist);

	Point.clear();

	for (size_t i = 0; TRIANGLE_TABLE[cubeIndex][i] != -1; i += 3)
	{
		size_t v[3];

		// Make Faces and Vertex ------------------------------------------------------------

		for (int j = 0; j < 3; j++)
		{
			if (mesh.Vertex.count(vertlist[TRIANGLE_TABLE[cubeIndex][i + j]])) 
			{
				v[j] = mesh.Vertex[vertlist[TRIANGLE_TABLE[cubeIndex][i + j]]];
			}
			else
			{
				v[j] = mesh.Vertex.size() + 1;
				//mesh.Vertex[vertlist[TRIANGLE_TABLE[cubeIndex][i + j]]] = v[j];
				mesh.push_Vertex(vertlist[TRIANGLE_TABLE[cubeIndex][i + j]], v[j]);
				mesh.VertexNormal.push_back({ { 0,0,1 } , 0 });
			}
		}
		//mesh.Faces.push_back(vector<size_t>{ v[0], v[1], v[2] });
		mesh.push_Faces(vector<size_t>{ v[0], v[1], v[2] });

		// Make Faces and Vertex ------------------------------------------------------------

		// CalculateNormal ------------------------------------------------------------------

		Coor p1 = vertlist[TRIANGLE_TABLE[cubeIndex][i]];
		Coor p2 = vertlist[TRIANGLE_TABLE[cubeIndex][i + 1]];
		Coor p3 = vertlist[TRIANGLE_TABLE[cubeIndex][i + 2]];

		// a = p2 - p1; b = p3 - p1
		Coor a = { p2.x - p1.x, p2.y - p1.y, p2.z - p1.z };
		Coor b = { p3.x - p1.x, p3.y - p1.y, p3.z - p1.z };

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
			if ((SizeTemp = mesh.VertexNormal[v[j] - 1].second) == 0)
				//mesh.VertexNormal[v[j] - 1] = { N, SizeTemp + 1 };
				mesh.insert_VertexNormal(v[j] - 1, N, SizeTemp + 1);
			else
			{
				// Mean of Vertex Normal
				Coor temp = mesh.VertexNormal[v[j] - 1].first;
				N = {	((SizeTemp * temp.x) + (N.x)) / (SizeTemp + 1)	,
						((SizeTemp * temp.y) + (N.y)) / (SizeTemp + 1)	,
						((SizeTemp * temp.z) + (N.z)) / (SizeTemp + 1)
					};

				// Normalization
				N.normalization();

				//mesh.VertexNormal[v[j] - 1] = { N, SizeTemp + 1 };
				mesh.insert_VertexNormal(v[j] - 1, N, SizeTemp + 1);
			}
		}

		// CalculateNormal ------------------------------------------------------------------
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
	
		//#pragma omp parallel for
	for (int z = 0; z < (ZSize - 1); z++)
	{
		YCoor = 0.0;
		for (int y = 0; y < (YSize - 1); y++)
		{
			XCoor = 0.0;
			for (int x = 0; x < (XSize - 1); x++)
			{
				vector<long double> PointValues	{
													Pixel_Array[z * xyWidth + y * xWidth + x],					// 0
													Pixel_Array[z * xyWidth + (y + 1) * xWidth + x],			// 1
													Pixel_Array[z * xyWidth + (y + 1) * xWidth + (x + 1)],		// 2
													Pixel_Array[z * xyWidth + y * xWidth + (x + 1)],			// 3
													Pixel_Array[(z + 1) * xyWidth + y * xWidth + x],			// 4
													Pixel_Array[(z + 1) * xyWidth + (y + 1) * xWidth + x],		// 5
													Pixel_Array[(z + 1) * xyWidth + (y + 1) * xWidth + (x + 1)],// 6
													Pixel_Array[(z + 1) * xyWidth + y * xWidth + (x + 1)]		// 7
												};

				Polygonise(PointValues, XCoor, YCoor, ZCoor, XCoor + XDist, YCoor + YDist, ZCoor + ZDist, Threshold, mesh);

				XCoor += XDist;
				PointValues.erase(PointValues.begin(), PointValues.end());
			}
			YCoor += YDist;
		}
		ZCoor += ZDist;
	}

	mesh.done();
}

size_t MakeOBJ(const char * name, Mesh & mesh)
{
	cout << "Make OBJ" << endl;

	errno = 0;

	ofstream Obj(name);

	ifstream Vertex("vCache");
	Obj << Vertex.rdbuf();
	Vertex.close();

	// Write Vertex Normal
	ifstream VertexNormal("vnCache");
	Obj << VertexNormal.rdbuf();
	VertexNormal.close();
	
	// Write Faces
	ifstream Faces("fCache");
	Obj << Faces.rdbuf();
	Faces.close();

	Obj.close();

	remove("vCache");
	remove("fCache");
	remove("vnCache");

	return (size_t)errno;
}
