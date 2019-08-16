#include "MarchingCube.h"

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

				size_t cubeIndex = 0;
				int shifter = 1;

				for (int i = 0; i < PointValues.size(); shifter <<= 1, i++)
					if (PointValues[i] < Threshold) cubeIndex |= shifter;
				//if (PointValues[i] > 0.0) cubeIndex |= shifter;

				if (EDGE_TABLE[cubeIndex] != 0)
				{
					vector<Coor> vertlist(12);

					vector<Coor> Point{
						{XCoor, YCoor, ZCoor},							// 0
						{XCoor, YCoor + YDist, ZCoor},					// 1
						{XCoor + XDist, YCoor + YDist, ZCoor},			// 2
						{XCoor + XDist, YCoor, ZCoor},					// 3
						{XCoor, YCoor, ZCoor + ZDist},					// 4
						{XCoor, YCoor + YDist, ZCoor + ZDist},			// 5
						{XCoor + XDist, YCoor + YDist, ZCoor + ZDist},	// 6
						{XCoor + XDist, YCoor, ZCoor + ZDist},			// 7
					};

					auto VertexInterp = [](Coor p1, Coor p2) -> Coor
					{
						return { (p1.x + p2.x) / 2,
									(p1.y + p2.y) / 2,
									(p1.z + p2.z) / 2 };
					};

					if (EDGE_TABLE[cubeIndex] & 1)		vertlist[0] = VertexInterp(Point[0], Point[1]);
					if (EDGE_TABLE[cubeIndex] & 2)		vertlist[1] = VertexInterp(Point[1], Point[2]);
					if (EDGE_TABLE[cubeIndex] & 4)		vertlist[2] = VertexInterp(Point[2], Point[3]);
					if (EDGE_TABLE[cubeIndex] & 8)		vertlist[3] = VertexInterp(Point[3], Point[0]);
					if (EDGE_TABLE[cubeIndex] & 16)		vertlist[4] = VertexInterp(Point[4], Point[5]);
					if (EDGE_TABLE[cubeIndex] & 32)		vertlist[5] = VertexInterp(Point[5], Point[6]);
					if (EDGE_TABLE[cubeIndex] & 64)		vertlist[6] = VertexInterp(Point[6], Point[7]);
					if (EDGE_TABLE[cubeIndex] & 128)	vertlist[7] = VertexInterp(Point[7], Point[4]);
					if (EDGE_TABLE[cubeIndex] & 256)	vertlist[8] = VertexInterp(Point[0], Point[4]);
					if (EDGE_TABLE[cubeIndex] & 512)	vertlist[9] = VertexInterp(Point[1], Point[5]);
					if (EDGE_TABLE[cubeIndex] & 1024)	vertlist[10] = VertexInterp(Point[2], Point[6]);
					if (EDGE_TABLE[cubeIndex] & 2048)	vertlist[11] = VertexInterp(Point[3], Point[7]);

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
								mesh.Vertex[vertlist[TRIANGLE_TABLE[cubeIndex][i + j]]] = v[j];
								mesh.VertexNormal.push_back({ { 0,0,1 } , 0 });
							}
						}
						mesh.push_Faces(vector<size_t>{ v[0], v[1], v[2] });
						//mesh.Faces.push_back(vector<size_t>{ v[0], v[1], v[2] });

						// Make Faces and Vertex ------------------------------------------------------------

						// Make Adj List --------------------------------------------------------------------

						for (int i = 0; i < 3; i++)
						{
							if ((v[i] - 1) == mesh.Adj.size()) mesh.Adj.push_back(set<size_t>());
							mesh.Adj[(v[i] - 1)].insert(v[(i + 1) % 3] - 1);
							mesh.Adj[(v[i] - 1)].insert(v[(i + 2) % 3] - 1);
						}

						// Make Adj List --------------------------------------------------------------------

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
								mesh.VertexNormal[v[j] - 1] = { N, SizeTemp + 1 };
							else
							{
								// Mean of Vertex Normal
								Coor temp = mesh.VertexNormal[v[j] - 1].first;
								N = { ((SizeTemp * temp.x) + (N.x)) / (SizeTemp + 1)	,
										((SizeTemp * temp.y) + (N.y)) / (SizeTemp + 1)	,
										((SizeTemp * temp.z) + (N.z)) / (SizeTemp + 1)
								};

								// Normalization
								N.normalization();

								mesh.VertexNormal[v[j] - 1] = { N, SizeTemp + 1 };
							}
						}

						// CalculateNormal ------------------------------------------------------------------
					}
				}

				XCoor += XDist;
				PointValues.erase(PointValues.begin(), PointValues.end());
			}
			YCoor += YDist;
		}
		ZCoor += ZDist;
	}

	mesh.done();
}

void Smoothing(Mesh & mesh, size_t rounds)
{
	map<size_t, Coor> new_verts;
	vector<pair<Coor, size_t> > new_norms(mesh.newVertex.size(), { {0.0, 0.0, 0.0} , 0 });

	cout << "Smoothing Start" << endl;

	for (size_t i = 0; i < rounds; i++)
	{
		#pragma omp parallel for
		for (long long int vert = 0; vert < mesh.newVertex.size(); vert++)
		{
			Coor new_vert = mesh.newVertex[vert];
			Coor new_norm = mesh.VertexNormal[vert].first;
			auto& neis = mesh.Adj[vert];

			for (auto nei : neis)
			{
				new_vert =	{
								new_vert.x + mesh.newVertex[nei].x,
								new_vert.y + mesh.newVertex[nei].y,
								new_vert.z + mesh.newVertex[nei].z,
							};

				new_norm =	{
								new_norm.x + mesh.VertexNormal[nei].first.x,
								new_norm.y + mesh.VertexNormal[nei].first.y,
								new_norm.z + mesh.VertexNormal[nei].first.z,
							};
			}

			new_verts[vert] =	{
									new_vert.x / (neis.size() + 1),
									new_vert.y / (neis.size() + 1),
									new_vert.z / (neis.size() + 1),
								};
			
			new_norms[vert].first = {
										new_norm.x / (neis.size() + 1),
										new_norm.y / (neis.size() + 1),
										new_norm.z / (neis.size() + 1),
									};
		}

		Swap< vector<pair<Coor, size_t> > >(mesh.VertexNormal, new_norms);
		Swap< map<size_t, Coor> >(mesh.newVertex, new_verts);
	}

	Swap< vector<pair<Coor, size_t> > >(mesh.VertexNormal, new_norms);
	Swap< map<size_t, Coor> >(mesh.newVertex, new_verts);
	mesh.VertexNormal = new_norms;
	mesh.newVertex = new_verts;
}


size_t MakeOBJ(const char * name, Mesh & mesh)
{
	cout << "Make OBJ" << endl;

	errno = 0;

	#pragma omp parallel sections
	{
		#pragma omp section
		{
			FILE *vCache = fopen("vCache", "w");
			fprintf(vCache, "\n# List of geometric vertices\n");

			for (auto v : mesh.newVertex)
				fprintf(vCache, "v %.3Lf %.3Lf %.3Lf\n", v.second.x, v.second.y, v.second.z);

			fclose(vCache);
		}

		#pragma omp section
		{
			FILE *vnCache = fopen("vnCache", "w");
			fprintf(vnCache, "\n# List of vertex normals\n");

			for (size_t i = 0; i < mesh.VertexNormal.size(); i++)
				fprintf(vnCache, "vn %.3Lf %.3Lf %.3Lf\n", mesh.VertexNormal[i].first.x, mesh.VertexNormal[i].first.y, mesh.VertexNormal[i].first.z);

			fclose(vnCache);
		}

		/*#pragma omp section
		{
			FILE *fCache = fopen("fCache", "w");
			fprintf(fCache, "\n# Polygonal face element\n");

			for (size_t i = 0; i < mesh.Faces.size(); i++)
				fprintf(fCache, "f %lld//%lld %lld//%lld %lld//%lld\n", mesh.Faces[i][0], mesh.Faces[i][0], mesh.Faces[i][1], mesh.Faces[i][1], mesh.Faces[i][2], mesh.Faces[i][2]);

			fclose(fCache);
		}*/
	}

	#pragma omp barrier

	mesh.remove();

	ofstream Obj(name);
	Obj << "# OBJ from MarchigCube by 3D BAMAG" << endl;

	ifstream vCache("vCache");
	Obj << vCache.rdbuf() << endl;
	vCache.close();

	ifstream vnCache("vnCache");
	Obj << vnCache.rdbuf() << endl;
	vnCache.close();

	ifstream fCache("fCache");
	Obj << fCache.rdbuf() << endl;
	fCache.close();

	Obj.close();

	remove("vCache");
	remove("vnCache");
	remove("fCache");

	return (size_t)errno;
}