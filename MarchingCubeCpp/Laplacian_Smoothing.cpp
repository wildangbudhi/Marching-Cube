#include "Laplacian_Smoothing.h"

template<typename T>
void Swap(T& a, T&b)
{
	T tmp = a;
	a = b;
	b = tmp;
}

void AdjList(Mesh & mesh, vector<vector<size_t>> & Adj)
{
	Adj.resize(mesh.Vertex.size(), vector<size_t>());

	#pragma omp parallel for
	for (long long int i = 0; i < mesh.Faces.size(); i++)
	{
		size_t a = mesh.Faces[i][0];
		size_t b = mesh.Faces[i][1];
		size_t c = mesh.Faces[i][2];

		Adj[a].push_back(b);
		Adj[a].push_back(c);
		Adj[b].push_back(a);
		Adj[b].push_back(c);
		Adj[c].push_back(a);
		Adj[c].push_back(b);
	}
}

void Smoothing(Mesh & mesh, size_t rounds)
{
	vector<vector<size_t>> Adj;
	AdjList(mesh, Adj);

	vector<Coor> new_verts(mesh.Vertex.size());
	vector<Coor> new_norms(mesh.Vertex.size());

	for (size_t i = 0; i < rounds; i++)
	{
		#pragma omp parallel for
		for (long long int vert = 0; vert < mesh.Vertex.size(); vert++)
		{
			// x
			double new_vert = mesh.Vertex[vert].x;
			double new_norm = mesh.VertexNormal[vert].x;
			auto& neis = Adj[vert];

			for (auto nei : neis)
			{
				new_vert += mesh.Vertex[nei].x;
				new_norm += mesh.VertexNormal[nei].x;
			}

			new_verts[vert].x = new_vert / (neis.size() + 1);
			new_norms[vert].x = new_norm / (neis.size() + 1);

			// y
			new_vert = mesh.Vertex[vert].y;
			new_norm = mesh.VertexNormal[vert].y;
			neis = Adj[vert];

			for (auto nei : neis)
			{
				new_vert += mesh.Vertex[nei].x;
				new_norm += mesh.VertexNormal[nei].x;
			}

			new_verts[vert].y = new_vert / (neis.size() + 1);
			new_norms[vert].y = new_norm / (neis.size() + 1);

			// z
			new_vert = mesh.Vertex[vert].z;
			new_norm = mesh.VertexNormal[vert].z;
			neis = Adj[vert];

			for (auto nei : neis)
			{
				new_vert += mesh.Vertex[nei].z;
				new_norm += mesh.VertexNormal[nei].z;
			}

			new_verts[vert].z = new_vert / (neis.size() + 1);
			new_norms[vert].z = new_norm / (neis.size() + 1);
		}

		Swap< vector<Coor> >(mesh.VertexNormal, new_norms);
		Swap< vector<Coor> >(mesh.Vertex, new_verts);
	}

	Swap< vector<Coor> >(mesh.VertexNormal, new_norms);
	Swap< vector<Coor> >(mesh.Vertex, new_verts);
}
