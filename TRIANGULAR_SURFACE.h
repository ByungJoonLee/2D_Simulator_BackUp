#pragma once

#include "COMMON_DEFINITION.h"
#include <list>
#include <vector>
#include "VECTOR_3D.h"
#include "ARRAY.h"
#include "ARRAY_VECTOR2.h"

// Triangular Geometry Components
struct VERTEX;
struct EDGE;
struct TRIANGLE;
struct HOLE;

class TRIANGULAR_SURFACE;

//////////////////////////////////////////////////////////////////////// 
//								VERTEX								 //
///////////////////////////////////////////////////////////////////////

struct VERTEX
{
public: // Essential Data
	T x[2];								// Position of vertices
	T n[2];								// Normal vector
	T s[2];								// Tension stress
	T stress;							// Determinant of tension stress
	T deviation[2];
	T normal_deviation;
	T velocity[2];
	T curvature_normal[2];
	T normalizer;

	int sample_direction;
	int index;

	bool feature;
	bool is_mc_vertex;
	bool is_boundary;
	
	list<TRIANGLE*> triangles;			// Triangles those includes this vertex

public: // Constructors and Destructor
	VERTEX(void);
	VERTEX(T* x);
	VERTEX(T* x, T* n);
	~VERTEX();
	
public: // Member Functions
	void SetPosition(T* x);
	void SetNormal(T* n);
	T*	 GetPosition();
	T*	 GetNormal();
	void AddTriangle(TRIANGLE* triangle);
	void DelTriangle(TRIANGLE* triangle);
	void DetermineNormal();
	void DrawNormal();
	void DrawCurvatureNormal();
	void DrawNeighborConnectivity();
	void DrawDeviation();
	void DrawVelocity();
	void Replace(VERTEX* v_after);
	void DetCurvatureNormal();

	VT GetVTPosition() 
	{
		return VT(x[0], x[1]);
	}

	VT GetVTNormal()
	{
		return VT(n[0], n[1]);
	}
};

//////////////////////////////////////////////////////////////////////// 
//								EDGE								 //
///////////////////////////////////////////////////////////////////////
struct EDGE
{
public: // Essential Data
	vector<VERTEX*> vertices;

public: // Constructor and Destructor
	EDGE(void);
	EDGE(VERTEX* v0, VERTEX* v1);
	~EDGE();

public: // Member Functions
	void AddVertex(VERTEX* v0, VERTEX* v1);
	void AddTriangle(TRIANGLE* triangle);
	bool IsSame(EDGE* edge);
	void Draw();
};

//////////////////////////////////////////////////////////////////////// 
//							   TRIANGLE								 //
///////////////////////////////////////////////////////////////////////
struct TRIANGLE
{
public: // Essential Data
	T n[2];
	VT uv[3];			// Texture Coordinate
	T area;
	VERTEX* vertices[3];
	VERTEX* edge_vertex[3];
	TRIANGLE* triangles[3];
	bool is_old;
	bool wrong;

	VI vertex_indices;

public: // Constructors and Destructor
	TRIANGLE(void);
	TRIANGLE(VERTEX* v0, VERTEX* v1, VERTEX* v2);
	~TRIANGLE(void);

public: // Member Functions
	void	DelTriangle(TRIANGLE* triangle);
	void	DetermineNormal();
	T*		GetNormal();
	void	SetNormal(T* n);
	void	DrawNormal();
	void	Draw();
	void	DrawBack();
	void	DrawEdges();
	void	DrawCenter();
	void	DrawNeighborConnectivity();
	void	CorrectCCW();
	int		CountFeatureVertex();
	void	FindIntersection(int direction, T* p, T* p_intersection, T* uv);
	bool	IsInside(VERTEX* v);
	void	Flip();
	void	Flip(VERTEX* v0, VERTEX* v1);
	VERTEX*	FindAnotherVertex(VERTEX* v0, VERTEX* v1);
	int		GetNeighborIndex(TRIANGLE* triangle);
	int		GetVertexIndex(VERTEX* v);
	void	ChkNeighborConnectivity();
	T		GetOppositeEdgeLength(VERTEX* v);
	void	AddLocalCurvatureNormal();

public: // Member Functions defined in this header file
	VT GetVTNormal()
	{
		return VT(n[0], n[1]);
	}

	VT GetVTCenter()
	{
		VT p0(vertices[0]->x[0], vertices[0]->x[1]);
		VT p1(vertices[1]->x[0], vertices[1]->x[1]);
		VT p2(vertices[2]->x[0], vertices[2]->x[1]);

		return (p0 + p1 + p2)/(T)3;
	}

	T SignedDistance(const VT& location) const
	{
		VT closest_pt = ClosestPoint(location);
		VT normal = VT(n[0], n[1]);

		return DotProduct(normal, location - closest_pt);
	}

	VT ClosestPoint(const VT& location) const
	{
		VT v0 = VT(vertices[0]->x[0], vertices[0]->x[1]);
		VT v1 = VT(vertices[1]->x[0], vertices[1]->x[1]);
		VT v2 = VT(vertices[2]->x[0], vertices[2]->x[1]);

		return ClosestPointFromTriangle(location, v0, v1, v2);
	}

	VT ClosestPoint(const VT& location, bool& on_triangle) const
	{
		VT v0 = VT(vertices[0]->x[0], vertices[0]->x[1]);
		VT v1 = VT(vertices[1]->x[0], vertices[1]->x[1]);
		VT v2 = VT(vertices[2]->x[0], vertices[2]->x[1]);

		return ClosestPointFromTriangle(location, v0, v1, v2, on_triangle);
	}

	static VECTOR_3D<T> BaryCentricCoordinates(const VT& location, const VT& x1, const VT& x2, const VT& x3)
	{
		VT u = x2 - x1, v = x3 - x1, w = location - x1;
		
		T u_dot_u = DotProduct(u, u), v_dot_v = DotProduct(v, v), u_dot_v = DotProduct(u, v),
		  u_dot_w = DotProduct(u, w), v_dot_w = DotProduct(v, w);

		T denominator = u_dot_u*v_dot_v - POW2(u_dot_v);
		T one_over_denominator;

		if(abs(denominator) > (T)1e-16)
		{
			one_over_denominator = 1/denominator;
		}
		else
		{
			one_over_denominator = (T)1e16;
		}

		T a = (v_dot_v*u_dot_w - u_dot_v*v_dot_w)*one_over_denominator, b = (u_dot_u*v_dot_w - u_dot_v*u_dot_w)*one_over_denominator;

		return VECTOR_3D<T>(1-a-b, a, b);
	}

	static VT ClosestPointFromTriangle(const VT& location, const VT& x1, const VT& x2, const VT& x3)
	{
		VT p = location;

		VECTOR_3D<T> weights = BaryCentricCoordinates(p, x1, x2, x3);

		if(weights.y < (T)0 && weights.z < (T)0)
		{
			return x1;
		}
		else if(weights.x < (T)0 && weights.z < (T)0)
		{
			return x2;
		}
		else if(weights.x < (T)0 && weights.y < (T)0)
		{
			return x3;
		}

		if(weights.x < (T)0)
		{
			return ClosestPointFromLine(p, x2, x3);
		}
		else if(weights.y < (T)0)
		{
			return ClosestPointFromLine(p, x1, x3);
		}
		else if(weights.z < (T)0)
		{
			return ClosestPointFromLine(p, x1, x2);
		}

		return p;
	}

	static VT ClosestPointFromTriangle(const VT& location, const VT& x1, const VT& x2, const VT& x3, bool& on_triangle)
	{
		VT p = location;

		VECTOR_3D<T> weights = BaryCentricCoordinates(p, x1, x2, x3);

		on_triangle = false;
		if(weights.y < (T)0 && weights.z < (T)0)
		{
			return x1;
		}
		else if(weights.x < (T)0 && weights.z < (T)0)
		{
			return x2;
		}
		else if(weights.x < (T)0 && weights.y < (T)0)
		{
			return x3;
		}

		if(weights.x < (T)0)
		{
			return ClosestPointFromLine(p, x2, x3);
		}
		else if(weights.y < (T)0)
		{
			return ClosestPointFromLine(p, x1, x3);
		}
		else if(weights.z < (T)0)
		{
			return ClosestPointFromLine(p, x1, x2);
		}

		on_triangle = true;
		return p;
	}
	
	static T SignedDistanceFromTriangle(const VT& location, const VT& x1, const VT& x2, const VT& x3)
	{
		VT closest_pt = ClosestPointFromTriangle(location, x1, x2, x3);
		T min_dist = (location - closest_pt).Magnitude();

		VECTOR_3D<T> weights = BaryCentricCoordinates(location, x1, x2, x3);
		if(weights.x >= 0 && weights.y >= 0 && weights.z >= 0) return -min_dist;
		else return min_dist;
	}

	static T SignedDistanceFromTriangle(const VT& location, const VT& x1, const VT& x2, const VT& x3, bool& on_triangle)
	{
		VT closest_pt = ClosestPointFromTriangle(location, x1, x2, x3, on_triangle);
		T min_dist = (location - closest_pt).Magnitude();

		VECTOR_3D<T> weights = BaryCentricCoordinates(location, x1, x2, x3);
		if(weights.x >= 0 && weights.y >= 0 && weights.z >= 0) return -min_dist;
		else return min_dist;
	}

	static VT ClosestPointFromLine(const VT& location, const VT& x1, const VT& x2)
	{
		T p = DotProduct(location - x1, x2 - x1)/DotProduct(x2 - x1, x2 - x1);

		if(p < (T)0)
		{
			return x1;
		}
		else if(p > (T)0)
		{
			return x2;
		}
		
		return x1 + (x2 - x1)*p;
	}

	/*static bool IntersectTriangle(const VT& i0, const VT& i1, const TRIANGLE& triangle, VT& intersection_position)
	{
		VT p0(triangle.vertices[0]->x);
		VT p1(triangle.vertices[1]->x);
		VT p2(triangle.vertices[2]->x);

		VT normal(triangle.n);
		normal.Normalize();

		VT q0 = i0 - p0;
		VT q1 = i1 - p0;

		T d0 = DotProduct(normal, q0);
		T d1 = DotProduct(normal, q1);

		VT closest_position = i0 - normal*d0;

		VT line_direction = i1 - i0;
		line_direction.Normalize();

		intersection_position = i0 + line_direction*DotProduct(line_direction, closest_position - i0);

		VT weights = BaryCentricCoordinates(intersection_position, p0, p1, p2);

		if(weights.x >= (T)0 && weights.y >= (T)0 && weights.z >= (T)0 && d0*d1 < (T)0) return true;

		return false;
	}

	static bool RayThruTriangle(const TRIANGLE& triangle, const VT& R1, const VT& R2, VT& intersection_position_out, const T& offset)
	{
		VT P1(triangle.vertices[0]->x);
		VT P2(triangle.vertices[1]->x);
		VT P3(triangle.vertices[2]->x);
		
		VT normal(triangle.n);
		normal.Normalize();

		VT intersect_position;

		T dist1 = DotProduct((R1 - P1), normal);
		T dist2 = DotProduct((R2 - P1), normal);

		if(dist1*dist2 >= (T)0 || dist1 == dist2) return false;

		intersect_position  = R1 + (R2 - R1)*(-dist1/(dist2 - dist1));

		VT vTest;

		vTest = CrossProduct(normal, P2 - P1);
		vTest.Normalize();
		if(DotProduct(vTest, intersect_position - P1) < 0) return false;

		vTest = CrossProduct(normal, P3 - P2);
		vTest.Normalize();
		if(DotProduct(vTest, intersect_position - P2) < 0) return false;

		vTest = CrossProduct(normal, P1 - P3);
		vTest.Normalize();
		if(DotProduct(vTest, intersect_position - P3) < 0) return false;

		intersection_position_out = intersect_position;

		return true;
	}*/
};

//////////////////////////////////////////////////////////////////////// 
//							    HOLE								 //
///////////////////////////////////////////////////////////////////////			
struct HOLE
{
public: // Essential Data
	list<VERTEX*> vertices;
	list<EDGE*> edges;
		
public: // Constructor and Destructor
	HOLE(void);
	~HOLE(void);

public: // Member Functions
	EDGE* AddEdge(VERTEX* v0, VERTEX* v1);
	EDGE* AddEdge(EDGE* edge);
	void Draw();
};

//////////////////////////////////////////////////////////////////////// 
//						   Mesh Manager								 //
///////////////////////////////////////////////////////////////////////			
class TRIANGULAR_SURFACE
{
public: // Essential Data
	vector<VERTEX*> vertices;
	list<EDGE*>		edges;
	list<TRIANGLE*> triangles;
	list<HOLE*>		holes;

	//BOX bounding_box;

public: // Constructor and Destructor
	TRIANGULAR_SURFACE(void);
	virtual ~TRIANGULAR_SURFACE(void);

public: // Member Functions
	EDGE*		AddEdge(VERTEX* v0, VERTEX* v1);
	TRIANGLE*	AddTriangle(int v0, int v1, int v2);
	TRIANGLE*	AddTriangle(VERTEX* v0, VERTEX* v1, VERTEX* v2);
	VERTEX*		AddVertex(T* x);
	VERTEX*		AddVertex(T* x, T* n);
	VERTEX*		AddVertex(VERTEX* vertex);
	void		DelAllTriangles();
	void		DelTriangle(TRIANGLE* );
	void		DetCurvatureNormals();
	void		DetVertexNormals();
	void		DetFaceNormals();
	void		AverageDuplexPositionNormal(TRIANGULAR_SURFACE* neighbor, float dx);
	void		DetermineNormalDeviation();
	void		DetTextureCoordinates(T* xy, T* uv, T s);
	void		DetTextureCoordinates(T x, T y, T u, T v, T s);
	void		ChkBoundary();
	void		ChkBoundaryVertices(T* size);
	void		ChkTrianglesNeighborConnectivity();
	VERTEX*		ChkNearestTextureCoordinate(T u, T v);
	void		DrawEdges();
	void		DrawFaceNormals();
	void		DrawHoles();
	void		DrawHoles(int index);
	void		DrawTriangles(const bool& draw_front = true);
	void		DrawVertices();
	void		DrawVertexNormals();
	void		DrawVertexDeviation();
	void		DrawTrianglesNeighborConnectivity();
	void		DrawTrianglesCenter();
	void		DrawVerticesNeighborConnectivity();
	void		DrawCurvatureNormal();
	void		DrawVertexVelocity();
	void		Filtering(T lamda);
	void		Filtering2(T lamda);
	//void		FilteringTangent(T lambda);
	//void		FilteringMinimumVariation(T lambda);
	//void		FilteringDeviation(T lambda);
	//void		Collapse(TRIANGLE* triangle);
	//void		Collapse(VERTEX* v0, VERTEX* v1);
	//void		FilpEdges();
	//void		CorrectCCW();
	//void		RemoveSmallTriangles(T threshold);
	//void		RemoveSmallEdges(T threshold);
	//void		RemoveLongTriangles(T ratio);
	//void		RemoveVertexDeviation();
	//void		UpdateBoundingBox();

	//void		Subdivision();
	//void		Scale(const T& scale);
	//void		Translate(T x, T y, T z);
	//void		FillHoles();
	//void		RemoveNonManifold();
	void		Reset();

	//void		MinMaxPosition(VT& min, VT& max);

	//// File I/O Functions
	//void		Read(const char* filename);
	//void		ReadOBJ(const char* filename, const T angle = (T)0, const VECTOR_3D<T> axis = VECTOR_3D<T>(1,1,1), const VECTOR_3D<T> scale = VECTOR_3D<T>(1,1,1), const VT translation = VT());
	//void		ReadSMF(const char* filename);
	//void		Write(const char* filename);
	//void		WriteOBJ(const char* filename);
	////void		WriteOBJ(const char* filename, VT& position, QUATERNION& quat);
	//void		WritePOLY(const char* filename);
	//static void	WriteOBJ(ARRAY<TRIANGULAR_SURFACE*>& arr, const char* filename);
	//// Realflow bin mesh file exporter
	//static void WriteBinMeshFile(ARRAY<TRIANGULAR_SURFACE*>& arr, int i_frame, const char* filename);
	//void		GetGeometryInfo(vector<VT>& vertices, vector<VI>& faces);
};
