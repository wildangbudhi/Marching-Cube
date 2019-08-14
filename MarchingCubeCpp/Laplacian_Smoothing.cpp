#include "Laplacian_Smoothing.h"

template<typename T>
void Swap(T& a, T&b)
{
	T tmp = a;
	a = b;
	b = tmp;
}

void Smoothing(Mesh & mesh, size_t rounds)
{
	//vector<Coor> new_verts(mesh.Vertex.size());
	//vector<pair<Coor, size_t> > new_norms(mesh.Vertex.size(), { {0.0, 0.0, 0.0} , 0 });

	//cout << "Smoothing Start" << endl;

	//for (size_t i = 0; i < rounds; i++)
	//{
	//	#pragma omp parallel for
	//	for (long long int vert = 0; vert < mesh.Vertex.size(); vert++)
	//	{
	//		// x
	//		double new_vert = mesh.Vertex[vert].x;
	//		double new_norm = mesh.VertexNormal[vert].first.x;
	//		auto& neis = Adj[vert];

	//		for (auto nei : neis)
	//		{
	//			new_vert += mesh.Vertex[nei].x;
	//			new_norm += mesh.VertexNormal[nei].first.x;
	//		}

	//		new_verts[vert].x = new_vert / (neis.size() + 1);
	//		new_norms[vert].first.x = new_norm / (neis.size() + 1);

	//		// y
	//		new_vert = mesh.Vertex[vert].y;
	//		new_norm = mesh.VertexNormal[vert].first.y;

	//		for (auto nei : neis)
	//		{
	//			new_vert += mesh.Vertex[nei].y;
	//			new_norm += mesh.VertexNormal[nei].first.y;
	//		}

	//		new_verts[vert].y = new_vert / (neis.size() + 1);
	//		new_norms[vert].first.y = new_norm / (neis.size() + 1);

	//		// z
	//		new_vert = mesh.Vertex[vert].z;
	//		new_norm = mesh.VertexNormal[vert].first.z;

	//		for (auto nei : neis)
	//		{
	//			new_vert += mesh.Vertex[nei].z;
	//			new_norm += mesh.VertexNormal[nei].first.z;
	//		}

	//		new_verts[vert].z = new_vert / (neis.size() + 1);
	//		new_norms[vert].first.z = new_norm / (neis.size() + 1);
	//	}

	//	Swap< vector<pair<Coor, size_t> > >(mesh.VertexNormal, new_norms);
	//	Swap< vector<Coor> >(mesh.Vertex, new_verts);
	//}

	//Swap< vector<pair<Coor, size_t> > >(mesh.VertexNormal, new_norms);
	//Swap< vector<Coor> >(mesh.Vertex, new_verts);
}
