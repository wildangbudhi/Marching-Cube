#include "Laplacian_Smoothing.h"

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

void Smoothing(Mesh & mesh)
{
	vector<vector<size_t>> Adj;
	AdjList(mesh, Adj);
}
