#pragma once

#include "TRIANGULAR_SURFACE.h"
#include "LEVELSET_2D.h"
#include "MARCHING_SQUARES_ALGORITHM_TABLE.h"

const static float node_lookup[4][2] = {{-1, 1}, {-1, -1}, {1, 1}, {1, -1},};
static const int Node_of_Edges[4][2] = {{0, 2}, {1, 3}, {1, 0}, {3, 2},};
static const int power_of_two[15] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};

class MARCHING_SQUARE
{
public:
	VERTEX* edgevertices[4];

public: // Constructor and Destructor
	MARCHING_SQUARE(void)
	{
		for (int i = 0; i < 4; i++)
		{
			edgevertices[i] = NULL;
		}
	}
};

class MARCHING_SQUARES_ALGORITHM
{
public: // Essential Data
	MULTITHREADING*				multithreading;
	ARRAY<TRIANGULAR_SURFACE*>	surfaces;
	ARRAY<VERTEX*>				vertices;
	ARRAY<EDGE*>			    edges;

	ARRAY<VERTEX**>				x_edge_vertices, y_edge_vertices, x_edge_vertices_j0;

	ARRAY<int>					vertex_start_indices, vertex_end_indices;
	ARRAY<int>					edge_start_indices, edge_end_indices;

	T							dx, dy;

	T							half_dx, half_dy, half_dz;
	T							quater_dx;

	int							num_threads;

	ARRAY<VT>					min_sizes, max_sizes;
	VT							global_translation;
	int							resolution_x;
	ARRAY<int>					resolution_y;

	VT							deviation;

	T							isocontour;
	
public: // Constructor and Destructor
	MARCHING_SQUARES_ALGORITHM(void)
	{}

	MARCHING_SQUARES_ALGORITHM(MULTITHREADING* multithreading_input, GRID_STRUCTURE_2D& grid, T grid_scale = (T)1, const T x_dev_input = (T)0, const T y_dev_input = (T)0)
		: isocontour((T)0)
	{
		Initialize(multithreading_input, grid, grid_scale, x_dev_input, y_dev_input);
	}

	~MARCHING_SQUARES_ALGORITHM(void)
	{
		if(surfaces.length != 0)
		{
			for(int i = 0; i < surfaces.length; i++)
				DELETE_POINTER(surfaces[i]);
		}
	}

public: // Initialization Function
	void Initialize(MULTITHREADING* multithreading_input, GRID_STRUCTURE_2D& grid, T grid_scale = (T)1, const T x_dev_input = (T)0, const T y_dev_input = (T)0)
	{
		multithreading = multithreading_input;

		int grid_max_res = MIN(grid.i_res, grid.j_res);

		if(multithreading->num_threads < grid_max_res)
			num_threads = multithreading->num_threads;
		else
			num_threads = grid_max_res;

		if(surfaces.length != 0)
		{
			for(int i = 0; i < surfaces.length; i++)
				DELETE_POINTER(surfaces[i]);
		}

		surfaces.Initialize(num_threads, 0);

		vertex_start_indices.Initialize(num_threads);
		vertex_end_indices.Initialize(num_threads);
		edge_start_indices.Initialize(num_threads);
		edge_end_indices.Initialize(num_threads);

		x_edge_vertices.Initialize(num_threads, 0);
		y_edge_vertices.Initialize(num_threads, 0);
		x_edge_vertices_j0.Initialize(num_threads, 0);
		
		ARRAY<GRID_STRUCTURE_2D> partial_grids;

		grid.SplitInYDirection(num_threads, partial_grids);

		dx = grid.dx/grid_scale;
		dy = grid.dy/grid_scale;

		half_dx = grid.dx*(T)0.5;
		half_dy = grid.dy*(T)0.5;

		quater_dx = MAX(dx, dy)*4;

		resolution_x = grid.i_res*(int)grid_scale;

		min_sizes.Initialize(num_threads);
		max_sizes.Initialize(num_threads);
		resolution_y.Initialize(num_threads);

		for(int i = 0; i < num_threads; i++)
		{
			resolution_y[i] = partial_grids[i].j_res;

			min_sizes[i].x = partial_grids[i].min[0];
			max_sizes[i].x = partial_grids[i].max[0];
			min_sizes[i].y = partial_grids[i].min[1];
			max_sizes[i].y = partial_grids[i].max[1];
			
			surfaces[i] = new TRIANGULAR_SURFACE();
		}

		deviation.x = x_dev_input;
		deviation.y = y_dev_input;

		global_translation = VT();

		isocontour = (T)0;
	}

public: // Member Functions
	int index(const int i, const int j) const
	{
		return i + resolution_x*j;
	}

	VT center(const int thread_id, const int i, const int j) const
	{
		VT resolution_offset = VT();
		return VT(min_sizes[thread_id].x + half_dx + dx*(T)i, min_sizes[thread_id].y + half_dy + dy*(T)j) + global_translation;
	}

	VT NodePosition(const int thread_id, const int n, const int i, const int j) const
	{
		VT center_position = center(thread_id, i, j);
		VT node_position;

		node_position.x = center_position.x + (T)node_lookup[n][0]*half_dx;
		node_position.y = center_position.y + (T)node_lookup[n][1]*half_dy;

		return node_position + deviation;
	}

	// Need to be fixed after implementation
	int index_edge_x(int i, int j) const
	{
		return i + resolution_x*j;
	}

	int index_edge_y(int i) const
	{
		return i;
	}

	VERTEX* GetEdgeVertex(int thread_id, int e, int i) const
	{
		if (e == 0)
		{
			return x_edge_vertices[thread_id][index_edge_x(i, 1)];
		}
		if (e == 1)
		{
			return x_edge_vertices[thread_id][index_edge_x(i, 0)];
		}
		if (e == 2)
		{
			return y_edge_vertices[thread_id][index_edge_y(i)];
		}
		if (e == 3)
		{
			return y_edge_vertices[thread_id][index_edge_y(i + 1)];
		}

		cout << "Null getEdgeVertex " << e << " " << i << endl;
		exit(1);
	}

	void SetEdgeVertex(int thread_id, int e, int i, VERTEX* vertex) const
	{
		if (e == 0)
		{
			x_edge_vertices[thread_id][index_edge_x(i, 1)] = vertex;
			return;
		}
		if (e == 1)
		{
			x_edge_vertices[thread_id][index_edge_x(i, 0)] = vertex;
			return;
		}
		if (e == 2)
		{
			y_edge_vertices[thread_id][index_edge_y(i)] = vertex;
			return;
		}
		if (e == 3)
		{
			y_edge_vertices[thread_id][index_edge_y(i + 1)] = vertex;
			return;
		}

		cout << "Null setEdgeVertex " << e << " " << i << endl;
		exit(1);
	}

	void PolygonizeUseFillGhostCellsFromThread(int& thread_id, LEVELSET_2D* levelset)
	{
		levelset->FillGhostCellsFrom(thread_id, levelset->phi, false);
	}

	void PolygonizeThread(int& thread_id, LEVELSET_OBJECT* levelset)
	{
		VERTEX **x_edge_vertices_old, **x_edge_vertices_temp;
		
		surfaces[thread_id]->Reset();

		DELETE_ARRAY(x_edge_vertices[thread_id]);

		int num_01 = resolution_x;
		int num_02 = resolution_x*2;
		int num_03 = resolution_x + 1;

		x_edge_vertices[thread_id] = new VERTEX* [num_02];
		y_edge_vertices[thread_id] = new VERTEX* [num_03];

		x_edge_vertices_old = new VERTEX* [num_02];
		
		x_edge_vertices_j0[thread_id] = new VERTEX* [num_01];

		for (int i = 0; i < num_01; i++)
		{
			x_edge_vertices_j0[thread_id][i] = 0;
		}

		for (int i = 0; i < num_02; i++)
		{
			x_edge_vertices_old[i] = 0;
		}

		for (int k = 0; k < resolution_y[thread_id]; k++)
		{
			for (int i = 0; i < num_02; i++)
			{
				x_edge_vertices[thread_id][i] = 0;
			}
			for (int i = 0; i < num_03; i++)
			{
				y_edge_vertices[thread_id][i] = 0;
			}

			for (int i = 0; i < resolution_x; i++)
			{
				x_edge_vertices[thread_id][index_edge_x(i, 0)] = x_edge_vertices_old[index_edge_x(i, 1)];
			}

			// Set Thread Boundary Edge
			if (resolution_y[thread_id] - 1 == k)
			{
				multithreading->Sync(thread_id, num_threads);
				if (thread_id != num_threads - 1)
				{
					for (int i = 0; i < resolution_x; i++)
					{
						x_edge_vertices[thread_id][index_edge_x(i, 1)] = x_edge_vertices_j0[thread_id + 1][index_edge_x(i, 0)];
					}
				}
			}
			
			for (int i = 0; i < resolution_x; i++)
			{
				T phi[4];
				bool is_skip = false;

				// Sample 4 node values at corners of a square
				for (int node = 0; node < 4; node++)
				{
					phi[node] = (*levelset)(NodePosition(thread_id, node, i, k));
					
					if (phi[node] == (T)0)
					{
						phi[node] = -(T)1e-8;				// To avoid making zero area triangles
					}

					if (phi[node] >= quater_dx)
					{
						is_skip = true;
						break;
					}
				}

				if (is_skip == true)
				{
					continue;
				}

				// Find square index
				int square_index(0);
				if (phi[0] <= (T)0)
				{
					square_index |= 1;
				}
				if (phi[1] <= (T)0)
				{
					square_index |= 2;
				}
				if (phi[2] <= (T)0)
				{
					square_index |= 4;
				}
				if (phi[3] <= (T)0)
				{
					square_index |= 8;
				}
				
				if (MSTABLES::EdgeTable[square_index] == 0)
				{
					continue;
				}

				// Make Vertices
				for (int edge = 0; edge < 4; edge++)
				{
					if (GetEdgeVertex(thread_id, edge, i) != 0)
					{
						continue;
					}

					if (MSTABLES::EdgeTable[square_index] && power_of_two[edge])
					{
						// Node indices of this edge
						const int n0 = Node_of_Edges[edge][0];
						const int n1 = Node_of_Edges[edge][1];

						// Phi values of edge nodes
						const T phi0 = phi[Node_of_Edges[edge][0]];
						const T phi1 = phi[Node_of_Edges[edge][1]];
						
						//if (phi0*phi1 < (T)0)
						{
							VT vertex_position = (abs(phi1)*NodePosition(thread_id, n0, i, k) + abs(phi0)*NodePosition(thread_id, n1, i, k))/(abs(phi0) + abs(phi1));
							VT vertex_normal = levelset->UnitNormal(vertex_position);

							VERTEX* vertex = new VERTEX(vertex_position.values, vertex_normal.values);

							SetEdgeVertex(thread_id, edge, i, vertex);

							surfaces[thread_id]->AddVertex(vertex);
						}
					}
				}
				
				for (int l = 0; l < 2; l++)
				{
					if (MSTABLES::LineTable[square_index][l*2 + 0] != -1)
					{
						surfaces[thread_id]->AddEdge(GetEdgeVertex(thread_id, MSTABLES::LineTable[square_index][l*2 + 0], i), GetEdgeVertex(thread_id, MSTABLES::LineTable[square_index][l*2 + 1], i));
					}
				}
			}// End of resolution_x loop
			
			if (k == 0)
			{
				for (int i = 0; i < resolution_x; i++)
				{
					x_edge_vertices_j0[thread_id][index_edge_x(i, 0)] = x_edge_vertices[thread_id][index_edge_x(i, 0)];
				}
			}

			x_edge_vertices_temp = x_edge_vertices_old;
			x_edge_vertices_old = x_edge_vertices[thread_id];
			x_edge_vertices[thread_id] = x_edge_vertices_temp;
			x_edge_vertices_temp = 0;
		}// End of resolution_y loop
		multithreading->Sync(thread_id, num_threads);

		//for (int i = 0; i < num_threads; i++)
		//{
		//	vertex_start_indices[thread_id] = 0;
		//	vertex_end_indices[thread_id] = 0;
		//	edge_start_indices[thread_id] = 0;
		//	edge_end_indices[thread_id] = 0;
		//}

		//for (int i = 0; i < thread_id; i++)
		//{
		//	vertex_start_indices[thread_id] += (int)surfaces[i]->vertices.size();
		//	edge_start_indices[thread_id] += (int)surfaces[i]->edges.size();
		//}
		//
		//for (int i = 0; i <= thread_id; i++)
		//{
		//	vertex_end_indices[thread_id] += (int)surfaces[i]->vertices.size();
		//	edge_end_indices[thread_id] += (int)surfaces[i]->edges.size();
		//}
		//multithreading->Sync(thread_id, num_threads);

		//if (thread_id == 0)
		//{
		//	vertices.Initialize(vertex_end_indices[num_threads - 1]);
		//	edges.Initialize(edge_end_indices[num_threads - 1]);

		//	// multithreading->SplitDomainIndex1D(0, triangles_end_indices[num_threads-1])
		//	const int j_start = 0;
		//	const int j_res = edge_end_indices[num_threads - 1];
		//	const int j_end = j_start + j_res - 1;
		//	const int quotient = j_res/num_threads;
		//	const int remainder = j_res%num_threads;

		//	int j_start_p = j_start;

		//	for(int i = 0; i < num_threads; i++)
		//	{
		//		int j_depth = i < remainder ? (quotient + 1) : quotient;

		//		multithreading->start_ix_1D[i] = j_start_p;
		//		multithreading->end_ix_1D[i] = j_start_p + j_depth - 1;
		//		
		//		j_start_p += j_depth;
		//	}
		//}
		//multithreading->Sync(thread_id, num_threads);

		//list<EDGE*>::iterator itr_edge;
		//itr_edge = surfaces[thread_id]->edges.begin();
		//for (int i = edge_start_indices[thread_id]; i < edge_end_indices[thread_id]; i++)
		//{
		//	edges[i] = (*itr_edge);
		//	itr_edge++;
		//}
		//multithreading->Sync(thread_id, num_threads);

		//int counter = 0;
		//for(int i = vertex_start_indices[thread_id]; i < vertex_end_indices[thread_id]; i++)
		//{
		//	vertices[i] = surfaces[thread_id]->vertices[counter];
		//	surfaces[thread_id]->vertices[counter]->index = i + 1;
		//	counter++;
		//}
		//multithreading->Sync(thread_id, num_threads);

		//if(thread_id == 0)
		//{
		//	// multithreading->SplitDomainIndex1D(0, vertex_end_indices[num_threads - 1])
		//	const int j_start = 0;
		//	const int j_res = vertex_end_indices[num_threads - 1];
		//	const int j_end = j_start + j_res - 1;
		//	const int quotient = j_res/num_threads;
		//	const int remainder = j_res%num_threads;

		//	int j_start_p = j_start;

		//	for(int i = 0; i < num_threads; i++)
		//	{
		//		int j_depth = i < remainder ? (quotient + 1) : quotient;

		//		multithreading->start_ix_1D[i] = j_start_p;
		//		multithreading->end_ix_1D[i] = j_start_p + j_depth - 1;

		//		j_start_p += j_depth;
		//	}
		//}
		//multithreading->Sync(thread_id, num_threads);

		//BEGIN_1D_ITERATION
		//{
		//	vertices[p]->DetermineNormal();
		//}
		//multithreading->Sync(thread_id, num_threads);}

		//DELETE_ARRAY(x_edge_vertices[thread_id]);
		//DELETE_ARRAY(y_edge_vertices[thread_id]);
		//DELETE_ARRAY(x_edge_vertices_old);
		//DELETE_ARRAY(x_edge_vertices_j0[thread_id]);
	}

	void Polygonize(LEVELSET_OBJECT& levelset)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
			multithreading->thread_list[thread_id] = new MULTITHREADING::THREAD(&MARCHING_SQUARES_ALGORITHM::PolygonizeThread, this, thread_id, &levelset);

		// multithreading->JoinAll()
		for(int i = 0; i < num_threads; i++) multithreading->thread_list[i]->join();

		// DeleteAllThread()
		for(int i = 0; i < num_threads; i++)
		{
			if(multithreading->thread_list[i] != 0)
			{
				delete multithreading->thread_list[i];
				multithreading->thread_list[i] = 0;
			}
		}
		multithreading->num_of_waiting_threads = 0;
	}

	void Polygonize(LEVELSET_2D& levelset, bool use_fill_ghost_cells_from)
	{
		if(use_fill_ghost_cells_from == true)
		{
			for(int thread_id = 0; thread_id < num_threads; thread_id++)
				multithreading->thread_list[thread_id] = new MULTITHREADING::THREAD(&MARCHING_SQUARES_ALGORITHM::PolygonizeUseFillGhostCellsFromThread, this, thread_id, &levelset);

			for(int i = 0; i < num_threads; i++) multithreading->thread_list[i]->join();

			for(int i = 0; i < num_threads; i++)
			{
				if(multithreading->thread_list[i] != 0)
				{
					delete multithreading->thread_list[i];
					multithreading->thread_list[i] = 0;
				}
			}
			multithreading->num_of_waiting_threads = 0;
		}

		for(int thread_id = 0; thread_id < num_threads; thread_id++)
			multithreading->thread_list[thread_id] = new MULTITHREADING::THREAD(&MARCHING_SQUARES_ALGORITHM::PolygonizeThread, this, thread_id, &levelset);

		for(int i = 0; i < num_threads; i++) multithreading->thread_list[i]->join();

		for(int i = 0; i < num_threads; i++)
		{
			if(multithreading->thread_list[i] != 0)
			{
				delete multithreading->thread_list[i];
				multithreading->thread_list[i] = 0;
			}
		}
		multithreading->num_of_waiting_threads = 0;
	}
};