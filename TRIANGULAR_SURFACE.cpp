#include "TRIANGULAR_SURFACE.h"

#include "GL/glut.h"
#include "GL/gl.h"

#include <math.h>
#include <algorithm>

#include <string.h>

#include "io.h"
#include <boost/format.hpp>
#include <direct.h>

#define TRAVERSE_VERTICES	vector<VERTEX*>::iterator itr_vertex;	for(itr_vertex = vertices.begin(); itr_vertex != vertices.end(); itr_vertex++)
#define TRAVERSE_EDGES		list<EDGE*>::iterator itr_edge;			for(itr_edge = edges.begin(); itr_edge != edges.end(); itr_edge++)
#define TRAVERSE_TRIANGLES  list<TRIANGLE*>::iterator itr_triangle; for(itr_triangle = triangles.begin(); itr_triangle != triangles.end(); itr_triangle++)

//////////////////////////////////////////////////////////////////////// 
//				 VERTEX Constructors and Destructor					 //
///////////////////////////////////////////////////////////////////////
VERTEX::VERTEX(void)
{
	feature = false;
	is_mc_vertex = false;
	is_boundary = false;
	x[0] = x[1] = 0;
	ARRAY_VECTOR2::set<T>(curvature_normal, (T)0);
	ARRAY_VECTOR2::set<T>(deviation, (T)0);
	ARRAY_VECTOR2::set<T>(velocity, (T)0);
	normalizer = (T)0;
	normal_deviation = (T)0;
}

VERTEX::VERTEX(T* x)
{
	feature = false;
	is_mc_vertex = false;
	is_boundary = false;
	SetPosition(x);
	ARRAY_VECTOR2::set<T>(curvature_normal, (T)0);
	ARRAY_VECTOR2::set<T>(deviation, (T)0);
	ARRAY_VECTOR2::set<T>(velocity, (T)0);
	normalizer = (T)0;
	normal_deviation = (T)0;
}

VERTEX::VERTEX(T* x, T* n)
{
	feature = false;
	is_mc_vertex = false;
	is_boundary = false;
	SetPosition(x);
	SetNormal(n);
	ARRAY_VECTOR2::set<T>(curvature_normal, (T)0);
	ARRAY_VECTOR2::set<T>(deviation, (T)0);
	ARRAY_VECTOR2::set<T>(velocity, (T)0);
	normalizer = (T)0;
	normal_deviation = (T)0;
}

VERTEX::~VERTEX(void)
{}

//////////////////////////////////////////////////////////////////////// 
//					VERTEX Member Functions						     //
///////////////////////////////////////////////////////////////////////
void VERTEX::SetPosition(T* x_input)
{
	ARRAY_VECTOR2::set<T>(x, x_input);
}

void VERTEX::SetNormal(T* n_input)
{
	ARRAY_VECTOR2::set<T>(n, n_input);
}

T* VERTEX::GetPosition()
{
	return x;
}

T* VERTEX::GetNormal()
{
	return n;
}

void VERTEX::AddTriangle(TRIANGLE* triangle)
{
	triangles.push_back(triangle);
}

void VERTEX::DelTriangle(TRIANGLE* triangle)
{
	triangles.remove(triangle);
}

void VERTEX::DetermineNormal()
{
	if((int)triangles.size() == 0)
	{
		n[0] = (T)0;
		n[1] = (T)0;
		n[2] = (T)0;
		return;
	}

	ARRAY_VECTOR2::set<T>(n, (T)0);

	TRAVERSE_TRIANGLES
	{
		ARRAY_VECTOR2::add<T>(n, (*itr_triangle)->GetNormal(), n);
	}

	ARRAY_VECTOR2::normalize<T>(VERTEX::n);

	ARRAY_VECTOR2::set<T>(curvature_normal, (T)0);			// Initialize curvature

	normalizer = (T)0;
}

void VERTEX::DrawNormal()
{
	glPushMatrix();
	glTranslatef(x[0], x[1], 0.0f);
	glBegin(GL_LINES);
		#ifdef USE_FLOAT_T
			glVertex2f(0.0f, 0.0f);
			glVertex2fv((float*)GetNormal());
		#else
			glVertex2d(0.0, 0.0);
			glVertex2dv(GetNormal());
		#endif
	glEnd();
	glPopMatrix();
}

void VERTEX::DrawCurvatureNormal()
{
	glPushMatrix();
	glTranslatef(x[0], x[1], 0.0f);
	T n[2];
	ARRAY_VECTOR2::set<T>(n, curvature_normal);
	ARRAY_VECTOR2::mul<T>(n, (T)1);
	glBegin(GL_LINES);
		#ifdef USE_FLOAT_T
			glVertex2f(0.0f, 0.0f);
			glVertex2fv((float*)n);
		#else
			glVertex2d(0.0, 0.0);
			glVertex2dv(n);
		#endif
	glEnd();
	glPopMatrix();
}

void VERTEX::DrawNeighborConnectivity()
{
	TRAVERSE_TRIANGLES
	{
		T v_center[3];
		ARRAY_VECTOR2::set<T>(v_center, (*itr_triangle)->vertices[0]->GetPosition());
		ARRAY_VECTOR2::add<T>(v_center, (*itr_triangle)->vertices[1]->GetPosition(), v_center);
		ARRAY_VECTOR2::add<T>(v_center, (*itr_triangle)->vertices[2]->GetPosition(), v_center);
		ARRAY_VECTOR2::div<T>(v_center, (T)3);

		glBegin(GL_LINES);
			#ifdef USE_FLOAT_T
				glVertex2fv((float*)GetPosition());
				glVertex2fv((float*)v_center);
			#else
				glVertex2dv(GetPosition());
				glVertex2dv(v_center);
			#endif
		glEnd();
	}
}

void VERTEX::DrawDeviation()
{
	if((int)triangles.size() == 0)
	{
		return;
	}
	glPushMatrix();
	glTranslatef(x[0], x[1], 0.0f);
	glBegin(GL_LINES);
		#ifdef USE_FLOAT_T
			glVertex2f(0.0f, 0.0f);
			glVertex2fv((float*)deviation);
		#else 
			glVertex2d(0.0, 0.0);
			glVertex2dv(deviation);
		#endif
	glEnd();
	glPopMatrix();
}

void VERTEX::DrawVelocity()
{
	glPushMatrix();
	glTranslatef(x[0], x[1], 0.0f);
	glBegin(GL_LINE);
		#ifdef USE_FLOAT_T
			glVertex2f(0.0f, 0.0f);
			glVertex2fv((float*)velocity);
		#else 
			glVertex2d(0.0, 0.0);
			glVertex2dv(velocity);
		#endif
	glEnd();
	glPopMatrix();
}

// After defining the Triangle structure, check this one more
void VERTEX::Replace(VERTEX* v_after)
{
	if(v_after == this)
	{
		cout << "Want to replace to the same vertex?" << endl;
		return;
	}

	cout << "Vertex replacing";
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;
		for(int i = 0; i < 3; i++)
		{
			if(triangle->vertices[i] == this)
			{
				triangle->vertices[i] = v_after;
				v_after->AddTriangle(triangle);
			}
		}
	}
	cout << "End" << endl;
	v_after->triangles.unique();
}

void VERTEX::DetCurvatureNormal()
{
	T A(0);
	TRAVERSE_TRIANGLES
	{
		A += (*itr_triangle)->area;
	}
	ARRAY_VECTOR2::div<T>(curvature_normal, (T)4*A);
}

//////////////////////////////////////////////////////////////////////// 
//					EDGE Constructor and Destructor  				 //
///////////////////////////////////////////////////////////////////////
EDGE::EDGE(void)
{}

EDGE::EDGE(VERTEX* v0, VERTEX* v1)
{
	vertices.push_back(v0);
	vertices.push_back(v1);
}

EDGE::~EDGE(void)
{}

//////////////////////////////////////////////////////////////////////// 
//					   EDGE Member Functions						 //
///////////////////////////////////////////////////////////////////////
void EDGE::AddVertex(VERTEX* v0, VERTEX* v1)
{
	vertices.push_back(v0);
	vertices.push_back(v1);
}

void EDGE::AddTriangle(TRIANGLE* triangle)
{
	cout << "EDGE::AddTriangle" << endl;
	exit(1);
}

bool EDGE::IsSame(EDGE* edge)
{
	/*VERTEX *v00, *v01, *v10, *v11;
	vector<VERTEX*>::iterator itr;

	itr = EDGE::vertices.begin();
	v00 = *itr;
	itr++;
	v01 = *itr;

	itr = edge->vertices.begin();
	v10 = *itr;
	itr++;
	v11 = *itr;

	if((v00 == v10) && (v01 == v11)) return true;
	else if((v00 == v11) && (v01 == v10)) return true;
	else return false;*/
}

void EDGE::Draw()
{
	
	/*glBegin(GL_LINES);
		
	glEnd();*/
	glBegin(GL_LINES);
		TRAVERSE_VERTICES
		{
			#ifdef USE_FLOAT_T
				glVertex2fv((float*)(*itr_vertex)->GetPosition());
			#else
				glVertex2dv((*itr_vertex)->GetPosition());
			#endif
		}
	glEnd();
}

//////////////////////////////////////////////////////////////////////// 
//	              TRIANGLE Constructor and Destructor  			     //
///////////////////////////////////////////////////////////////////////
TRIANGLE::TRIANGLE(void)
{}

TRIANGLE::TRIANGLE(VERTEX* v0, VERTEX* v1, VERTEX* v2)
{
	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v2;

	edge_vertex[0] = NULL;
	edge_vertex[1] = NULL;
	edge_vertex[2] = NULL;

	is_old = false;
	wrong = false;
}

TRIANGLE::~TRIANGLE(void)
{}

//////////////////////////////////////////////////////////////////////// 
//				      TRIANGLE Member Functions						 //
///////////////////////////////////////////////////////////////////////
void TRIANGLE::DelTriangle(TRIANGLE* triangle)
{
	for(int i = 0; i < 3; i++)
	{
		if(triangles[i] == triangle)
		{
			triangles[i] = NULL;
			return;
		}
	}
}

void TRIANGLE::DetermineNormal()
{
	T l0[2];
	T l1[2];
	ARRAY_VECTOR2::sub<T>(TRIANGLE::vertices[2]->GetPosition(), TRIANGLE::vertices[0]->GetPosition(), l0);
	ARRAY_VECTOR2::sub<T>(TRIANGLE::vertices[1]->GetPosition(), TRIANGLE::vertices[0]->GetPosition(), l1);
	ARRAY_VECTOR2::cross<T>(l1, l0, TRIANGLE::n);
	TRIANGLE::area = ARRAY_VECTOR2::det<T>(TRIANGLE::n)/2.0f;

	ARRAY_VECTOR2::normalize<T>(n);
}

T* TRIANGLE::GetNormal()
{
	return n;
}

void TRIANGLE::SetNormal(T* n_input)
{
	ARRAY_VECTOR2::set<T>(TRIANGLE::n, n_input);
}

void TRIANGLE::DrawNormal()
{
	T p[2];
	ARRAY_VECTOR2::set<T>(p, (T)0);
	ARRAY_VECTOR2::add<T>(p, TRIANGLE::vertices[0]->GetPosition(), p);
	ARRAY_VECTOR2::add<T>(p, TRIANGLE::vertices[1]->GetPosition(), p);
	ARRAY_VECTOR2::add<T>(p, TRIANGLE::vertices[2]->GetPosition(), p);
	ARRAY_VECTOR2::div<T>(p, 3.0f);

	glPushMatrix();
	glTranslatef(p[0], p[1], 0.0f);
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0, 0.0f);
		glVertex3fv((float*)TRIANGLE::GetNormal());
	glEnd();
	glPopMatrix();
}

void TRIANGLE::Draw()
{
	glBegin(GL_TRIANGLES);
	{
		#ifdef USE_FLOAT_T
			glNormal3fv(vertices[0]->GetNormal());
			glTexCoord2f(uv[0].x, uv[0].y);
			glVertex3fv(vertices[0]->GetPosition());
			glNormal3fv(vertices[1]->GetNormal());
			glTexCoord2f(uv[1].x, uv[1].y);
			glVertex3fv(vertices[1]->GetPosition());
			glNormal3fv(vertices[2]->GetNormal());
			glTexCoord2f(uv[2].x, uv[2].y);
			glVertex3fv(vertices[2]->GetPosition());
		#else
			glNormal3dv(vertices[0]->GetNormal());
			glTexCoord2d(uv[0].x, uv[0].y);
			glVertex3dv(vertices[0]->GetPosition());
			glNormal3dv(vertices[1]->GetNormal());
			glTexCoord2d(uv[1].x, uv[1].y);
			glVertex3dv(vertices[1]->GetPosition());
			glNormal3dv(vertices[2]->GetNormal());
			glTexCoord2d(uv[2].x, uv[2].y);
			glVertex3dv(vertices[2]->GetPosition());
		#endif
	}
	glEnd();
}

void TRIANGLE::DrawBack()
{
	glBegin(GL_TRIANGLES);
	{
		#ifdef USE_FLOAT_T
			glNormal3f(-vertices[0]->GetNormal()[0], -vertices[0]->GetNormal()[1], -vertices[0]->GetNormal()[2]);
			glTexCoord2f(uv[0].x, uv[0].y);
			glVertex3fv((float*)vertices[0]->GetPosition());
			glNormal3f(-vertices[2]->GetNormal()[0], -vertices[2]->GetNormal()[1], -vertices[2]->GetNormal()[2]);
			glTexCoord2f(uv[2].x, uv[2].y);
			glVertex3fv((float*)vertices[2]->GetPosition());
			glNormal3f(-vertices[1]->GetNormal()[0], -vertices[1]->GetNormal()[1], -vertices[1]->GetNormal()[2]);
			glTexCoord2f(uv[1].x, uv[1].y);
			glVertex3fv((float*)vertices[1]->GetPosition());
		#else
			glNormal3d(-vertices[0]->GetNormal()[0], -vertices[0]->GetNormal()[1], -vertices[0]->GetNormal()[2]);
			glTexCoord2d(uv[0].x, uv[0].y);
			glVertex3dv((double*)vertices[0]->GetPosition());
			glNormal3d(-vertices[2]->GetNormal()[0], -vertices[2]->GetNormal()[1], -vertices[2]->GetNormal()[2]);
			glTexCoord2d(uv[2].x, uv[2].y);
			glVertex3dv((double*)vertices[2]->GetPosition());
			glNormal3d(-vertices[1]->GetNormal()[0], -vertices[1]->GetNormal()[1], -vertices[1]->GetNormal()[2]);
			glTexCoord2d(uv[1].x, uv[1].y);
			glVertex3dv((double*)vertices[1]->GetPosition());
		#endif
	}
	glEnd();
}

void TRIANGLE::DrawEdges()
{
	glBegin(GL_LINES);
	{
		#ifdef USE_FLOAT_T
			glVertex2fv(TRIANGLE::vertices[1]->GetPosition());
			glVertex2fv(TRIANGLE::vertices[2]->GetPosition());
		#else
			glVertex2dv(TRIANGLE::vertices[1]->GetPosition());
			glVertex2dv(TRIANGLE::vertices[2]->GetPosition());
		#endif
	}
	glEnd();

	glBegin(GL_LINES);
	{
		#ifdef USE_FLOAT_T
			glVertex2fv(TRIANGLE::vertices[2]->GetPosition());
			glVertex2fv(TRIANGLE::vertices[0]->GetPosition());
		#else
			glVertex2dv(TRIANGLE::vertices[2]->GetPosition());
			glVertex2dv(TRIANGLE::vertices[0]->GetPosition());
		#endif
	}
	glEnd();
	
	glBegin(GL_LINES);
	{
		#ifdef USE_FLOAT_T
			glVertex2fv(TRIANGLE::vertices[0]->GetPosition());
			glVertex2fv(TRIANGLE::vertices[1]->GetPosition());
		#else
			glVertex2dv(TRIANGLE::vertices[0]->GetPosition());
			glVertex2dv(TRIANGLE::vertices[1]->GetPosition());
		#endif
	}
	glEnd();
}

void TRIANGLE::DrawCenter()
{
	T center[2];
	ARRAY_VECTOR2::set<T>(center, (T)0);
	ARRAY_VECTOR2::add<T>(center, vertices[0]->GetPosition(), center);
	ARRAY_VECTOR2::add<T>(center, vertices[1]->GetPosition(), center);
	ARRAY_VECTOR2::add<T>(center, vertices[2]->GetPosition(), center);
	ARRAY_VECTOR2::div<T>(center, 3.0f);

	glBegin(GL_POINTS);
		glVertex2fv((float*)center);
	glEnd();
}

void TRIANGLE::DrawNeighborConnectivity()
{
	T center[2];
	ARRAY_VECTOR2::set<T>(center, (T)0);
	ARRAY_VECTOR2::add<T>(center, vertices[0]->GetPosition(), center);
	ARRAY_VECTOR2::add<T>(center, vertices[1]->GetPosition(), center);
	ARRAY_VECTOR2::add<T>(center, vertices[2]->GetPosition(), center);
	ARRAY_VECTOR2::div<T>(center, 3.0f);

	for(int i = 0; i < 3; i++)
	{
		if(TRIANGLE::triangles[i] == NULL)
		{
			cout << "No triangles are connected" << endl;
			continue;
		}
		T neighbor[2];
		
		ARRAY_VECTOR2::set<T>(neighbor, (T)0);
		ARRAY_VECTOR2::add<T>(neighbor, TRIANGLE::triangles[i]->vertices[0]->GetPosition(), neighbor);
		ARRAY_VECTOR2::add<T>(neighbor, TRIANGLE::triangles[i]->vertices[1]->GetPosition(), neighbor);
		ARRAY_VECTOR2::add<T>(neighbor, TRIANGLE::triangles[i]->vertices[2]->GetPosition(), neighbor);
		ARRAY_VECTOR2::div<T>(neighbor, 3.0f);

		ARRAY_VECTOR2::add<T>(center, neighbor, neighbor);
		ARRAY_VECTOR2::div<T>(neighbor, 2.0f);

		glBegin(GL_LINES);
			glVertex2fv((float*)center);
			glVertex2fv((float*)neighbor);
		glEnd();
	}
}

// Need to understand this when you have a time
void TRIANGLE::CorrectCCW()
{
	T l0[2];
	T l1[2];

	ARRAY_VECTOR2::sub<T>(TRIANGLE::vertices[2]->GetPosition(), TRIANGLE::vertices[0]->GetPosition(), l0);
	ARRAY_VECTOR2::sub<T>(TRIANGLE::vertices[1]->GetPosition(), TRIANGLE::vertices[0]->GetPosition(), l1);

	T c[3];
	ARRAY_VECTOR2::cross<T>(l1, l0, c);				// Face normal vector

	TRIANGLE::area = ARRAY_VECTOR2::det(c)/2.0f;	// Calculate triangle area

	T n[3];
	ARRAY_VECTOR2::add<T>(TRIANGLE::vertices[0]->GetNormal(), TRIANGLE::vertices[1]->GetNormal(), n);
	ARRAY_VECTOR2::add<T>(TRIANGLE::vertices[2]->GetNormal(), n, n);

	if((c[0]*n[0] + c[1]*n[1] + c[2]*n[2]) < (T)0)
	{
		VERTEX* temp = TRIANGLE::vertices[1];
		TRIANGLE::vertices[1] = TRIANGLE::vertices[2];
		TRIANGLE::vertices[2] = temp;
	}
}

int TRIANGLE::CountFeatureVertex()
{
	int number = 0;
	for(int i = 0; i < 3; i++)
	{
		if(vertices[i]->feature == true)
		{
			number++;
		}
	}
	return number;
}

void TRIANGLE::FindIntersection(int direction, T* p_input, T* p_intersection, T* uv)
{
	T l0[2];
	T l1[2];

	ARRAY_VECTOR2::sub<T>(TRIANGLE::vertices[1]->x, TRIANGLE::vertices[0]->x, l0);
	ARRAY_VECTOR2::sub<T>(TRIANGLE::vertices[2]->x, TRIANGLE::vertices[0]->x, l1);

	T l0x = l0[(direction+1)%3];
	T l0y = l0[(direction+2)%3];
	T l1x = l0[(direction+1)%3];
	T l1y = l0[(direction+2)%3];

	T p[2];
	p[0] = p_input[(direction+1)%3] - vertices[0]->x[(direction+1)%3];
	p[0] = p_input[(direction+2)%3] - vertices[0]->x[(direction+2)%3];

	T det = l0x*l1y - l1x*l0y;
	det = 1.0f/det;
	uv[0] = det*(l1y*p[0] - l1x*p[1]);
	uv[1] = det*(-l0y*p[0] + l0x*p[1]);

	ARRAY_VECTOR2::set<T>(p_intersection, vertices[0]->x);

	p_intersection[0] += l0[0]*uv[0] + l1[0]*uv[1];
	p_intersection[1] += l0[1]*uv[0] + l1[1]*uv[1];
	
	return;
}

bool TRIANGLE::IsInside(VERTEX* v)
{
	if(vertices[0] == v)
	{
		return true;
	}
	if(vertices[1] == v)
	{
		return true;
	}
	if(vertices[2] == v)
	{
		return true;
	}
	return false;
}

void TRIANGLE::Flip()
{}

void TRIANGLE::Flip(VERTEX* v0, VERTEX* v1)
{}

VERTEX* TRIANGLE::FindAnotherVertex(VERTEX* v0, VERTEX* v1)
{
	for(int i = 0; i < 3; i++)
	{
		if(vertices[i] != v0 && vertices[i] != v1)
		{
			return vertices[i];
		}
	}
	
	return NULL;
}

int TRIANGLE::GetNeighborIndex(TRIANGLE* triangle)
{
	for(int i = 0; i < 3; i++)
	{
		if(triangles[i] == triangle)
		{
			return i;
		}
	}
	cout << "TRIANGLE::GetNeighborIndex error!" << endl;
	return 0;
}

int TRIANGLE::GetVertexIndex(VERTEX* v)
{
	for(int i = 0; i < 3; i++)
	{
		if(vertices[i] == v)
		{
			return i;
		}
	}
	cout << "TRIANGLE::GetVertexIndex error!" << endl;
	return -1;
}

void TRIANGLE::ChkNeighborConnectivity()
{
	for(int i = 0; i < 3; i++)
	{
		triangles[i] = NULL;
		{
			list<TRIANGLE*>::iterator itr;
			vertices[(i+1)%3]->DelTriangle(this);
			vertices[(i+2)%3]->DelTriangle(this);
			for(itr = vertices[(i+1)%3]->triangles.begin(); itr != vertices[(i+1)%3]->triangles.end(); itr++)
			{
				if(find(vertices[(i+2)%3]->triangles.begin(), vertices[(i+2)%3]->triangles.end(), *itr) != vertices[(i+2)%3]->triangles.end())
				{
					triangles[i] = *itr;
					break;
				}
			}
			vertices[(i+1)%3]->AddTriangle(this);
			vertices[(i+2)%3]->AddTriangle(this);
		}
	}

	for(int i = 0; i < 3; i++)
	{
		if(triangles[i] == NULL)
		{
			continue;
		}
		if(triangles[i] == triangles[(i+1)%3])
		{
			wrong = true;
			triangles[i]->wrong = true;
			triangles[(i+1)%3]->wrong = true;
		}
		if(triangles[i] == triangles[(i+2)%3])
		{
			wrong = true;
			triangles[i]->wrong = true;
			triangles[(i+2)%3]->wrong = true;
		}
		if(vertices[i] == vertices[(i+1)%3] || vertices[i] == vertices[(i+2)%3])
		{
			cout << "Triangle vertices duplication" << endl;
		}
	}
}

T TRIANGLE::GetOppositeEdgeLength(VERTEX* v)
{
	int i;
	for(i = 0;i < 3; i++)
	{
		if(TRIANGLE::vertices[i] == v)
		{
			break;
		}
	}

	T length[2];
	ARRAY_VECTOR2::set<T>(length, vertices[(i+1)%3]->GetPosition());
	ARRAY_VECTOR2::sub<T>(length, vertices[(i+2)%3]->GetPosition(), length);

	return ARRAY_VECTOR2::det(length);
}

void TRIANGLE::AddLocalCurvatureNormal()
{
	for(int i = 0; i < 3; i++)
	{
		T l0[2], l1[2];
		ARRAY_VECTOR2::set<T>(l0, vertices[(i+1)%3]->GetPosition());
		ARRAY_VECTOR2::sub<T>(l0, vertices[i]->GetPosition(), l0);
		ARRAY_VECTOR2::set<T>(l1, vertices[(i+2)%3]->GetPosition());
		ARRAY_VECTOR2::sub<T>(l1, vertices[i]->GetPosition(), l1);

		T length[2];
		ARRAY_VECTOR2::set<T>(length, vertices[(i+1)%3]->GetPosition());
		ARRAY_VECTOR2::sub<T>(length, vertices[(i+2)%3]->GetPosition(), length);

		T weight = ARRAY_VECTOR2::det(length);

		ARRAY_VECTOR2::mul<T>(l0, weight);
		ARRAY_VECTOR2::mul<T>(l1, weight);
		ARRAY_VECTOR2::add<T>(vertices[i]->curvature_normal, l0, vertices[i]->curvature_normal);
		ARRAY_VECTOR2::add<T>(vertices[i]->curvature_normal, l1, vertices[i]->curvature_normal);
		vertices[i]->normalizer += weight;
	}
}

//////////////////////////////////////////////////////////////////////// 
//	                HOLE Constructor and Destructor  			     //
///////////////////////////////////////////////////////////////////////
HOLE::HOLE(void)
{}

HOLE::~HOLE(void)
{}

//////////////////////////////////////////////////////////////////////// 
//	                       HOLE Member Functions  			         //
///////////////////////////////////////////////////////////////////////
EDGE* HOLE::AddEdge(VERTEX* v0, VERTEX* v1)
{
	EDGE* edge = new EDGE(v0, v1);
	HOLE::edges.push_back(edge);
	return edge;
}

EDGE* HOLE::AddEdge(EDGE* edge)
{
	HOLE::edges.push_back(edge);
	return edge;
}

void HOLE::Draw()
{
	TRAVERSE_EDGES
	{
		(*itr_edge)->Draw();
	}
}

//////////////////////////////////////////////////////////////////////// 
//	          Mesh Manager Constructor and Destructor  			     //
///////////////////////////////////////////////////////////////////////	
TRIANGULAR_SURFACE::TRIANGULAR_SURFACE(void)
{}

TRIANGULAR_SURFACE::~TRIANGULAR_SURFACE(void)
{
	Reset();
}

//////////////////////////////////////////////////////////////////////// 
//				     Mesh Manager Member Functions  			     //
///////////////////////////////////////////////////////////////////////	
EDGE* TRIANGULAR_SURFACE::AddEdge(VERTEX* v0, VERTEX* v1)
{
	EDGE* edge = new EDGE(v0, v1);
	edges.push_back(edge);
	
	return edge;
}

TRIANGLE* TRIANGULAR_SURFACE::AddTriangle(int v0, int v1, int v2)
{
	return TRIANGULAR_SURFACE::AddTriangle(TRIANGULAR_SURFACE::vertices[v0], TRIANGULAR_SURFACE::vertices[v1], TRIANGULAR_SURFACE::vertices[v2]);
}

TRIANGLE* TRIANGULAR_SURFACE::AddTriangle(VERTEX* v0, VERTEX* v1, VERTEX* v2)
{
	VERTEX* v[3] = {v0, v1, v2};
	TRIANGLE* triangle = new TRIANGLE(v[0], v[1], v[2]);
	triangles.push_back(triangle);
	
	v[0]->AddTriangle(triangle);
	v[1]->AddTriangle(triangle);
	v[2]->AddTriangle(triangle);

	return triangle;
}

VERTEX* TRIANGULAR_SURFACE::AddVertex(T* x)
{
	VERTEX* vertex = new VERTEX(x);
	vertices.push_back(vertex);
	
	return vertex;
}

VERTEX* TRIANGULAR_SURFACE::AddVertex(T* x, T* n)
{
	VERTEX* vertex = new VERTEX(x, n);
	vertices.push_back(vertex);

	return vertex;
}

VERTEX* TRIANGULAR_SURFACE::AddVertex(VERTEX* vertex)
{
	vertices.push_back(vertex);

	return vertex;
}

void TRIANGULAR_SURFACE::DelAllTriangles()
{
	list<TRIANGLE*>::iterator itr_triangle = TRIANGULAR_SURFACE::triangles.begin();
	
	while(itr_triangle != TRIANGULAR_SURFACE::triangles.end())
	{
		TRIANGLE* triangle = *itr_triangle;
		for(int i = 0; i < 3; i++)
		{
			triangle->vertices[i]->DelTriangle(triangle);
			if(triangle->triangles[i] != NULL)
			{
				triangle->triangles[i]->DelTriangle(triangle);
				triangle->triangles[i] = NULL;
			}
		}

		delete triangle;

		triangles.erase(itr_triangle);

		itr_triangle = triangles.begin();
	}
}

void TRIANGULAR_SURFACE::DelTriangle(TRIANGLE* triangle)
{
	triangles.remove(triangle);

	for(int i = 0; i < 3; i++)
	{
		triangle->vertices[i]->DelTriangle(triangle);
		if(triangle->triangles[i] != NULL)
		{
			triangle->triangles[i]->DelTriangle(triangle);
		}
	}

	delete triangle;
}

void TRIANGULAR_SURFACE::DetCurvatureNormals()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->AddLocalCurvatureNormal();
	}
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DetCurvatureNormal();
	}
}

void TRIANGULAR_SURFACE::DetVertexNormals()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DetermineNormal();
	}
}

void TRIANGULAR_SURFACE::DetFaceNormals()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DetermineNormal();
	}
}

void TRIANGULAR_SURFACE::AverageDuplexPositionNormal(TRIANGULAR_SURFACE* neighbor, float dx)
{
	T critical_value = dx/4.0;
	for(unsigned int i = 0; i < vertices.size(); i++)
	{
		VERTEX* i_vertex = vertices[i];

		for(unsigned int j = 0; j < neighbor->vertices.size(); j++)
		{
			VERTEX* j_vertex = neighbor->vertices[j];
			float dis[2] = {i_vertex->x[0] - j_vertex->x[0], i_vertex->x[1] - j_vertex->x[1]};
			if(abs(dis[0]) < critical_value && abs(dis[1]) < critical_value)
			{
				float avrNormal[3] = {(i_vertex->n[0] + j_vertex->n[0])/2, (i_vertex->n[1] + j_vertex->n[1])/2, (i_vertex->n[2] + j_vertex->n[2])/2 };
				i_vertex->n[0] = avrNormal[0]; i_vertex->n[1] = avrNormal[1]; i_vertex->n[2] = avrNormal[2];
				j_vertex->n[0] = avrNormal[0]; j_vertex->n[1] = avrNormal[1]; j_vertex->n[2] = avrNormal[2];
			}
		}
	}
}

void TRIANGULAR_SURFACE::DetermineNormalDeviation()
{
	TRAVERSE_VERTICES
	{
		T detdeviation = ARRAY_VECTOR2::det<T>((*itr_vertex)->deviation);
		if(ARRAY_VECTOR2::dot<T>((*itr_vertex)->n, (*itr_vertex)->deviation) < (T)0)
		{
			detdeviation *= -(T)1;
		}

		(*itr_vertex)->normal_deviation = detdeviation;
		(*itr_vertex)->deviation[0] = (*itr_vertex)->normal_deviation*(*itr_vertex)->n[0];
		(*itr_vertex)->deviation[1] = (*itr_vertex)->normal_deviation*(*itr_vertex)->n[1];
	}
}

void TRIANGULAR_SURFACE::DetTextureCoordinates(T* xy, T* uv, T s)
{
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;

		for(int i = 0; i < 3; i++)
		{
			VT& uv = triangle->uv[i];
			VERTEX* vertex = triangle->vertices[i];
			uv.x = (vertex->GetPosition()[0] - xy[0])/s*uv.x;
			uv.y = (vertex->GetPosition()[1] - xy[1])/s*uv.y;
			if(uv.x > (T)1)
			{
				uv.x = (T)1;
			}
			if(uv.y > (T)1)
			{
				uv.y = (T)1;
			}
			if(uv.x < (T)0)
			{
				uv.x = (T)0;
			}
			if(uv.y < (T)0)
			{
				uv.y = (T)0;
			}
		}
	}
}

void TRIANGULAR_SURFACE::DetTextureCoordinates(T x, T y, T u, T v, T s)
{
	T xy[2] = {x, y};
	T uv[2] = {u, v};

	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;

		for(int i = 0; i < 3; i++)
		{
			VT& uv = triangle->uv[i];
			VERTEX* vertex = triangle->vertices[i];
			uv.x = (vertex->GetPosition()[0] - xy[0])/s*uv.x;
			uv.y = (vertex->GetPosition()[1] - xy[1])/s*uv.y;
			if(uv.x > (T)1)
			{
				uv.x = (T)1;
			}
			if(uv.y > (T)1)
			{
				uv.y = (T)1;
			}
			if(uv.x < (T)0)
			{
				uv.x = (T)0;
			}
			if(uv.y < (T)0)
			{
				uv.y = (T)0;
			}
		}
	}
}

void TRIANGULAR_SURFACE::ChkBoundary()
{
	cout << "# Boundary Checking Started" << endl;

	TRAVERSE_TRIANGLES
	{
		if((*itr_triangle)->triangles[0] == NULL)
		{
			AddEdge((*itr_triangle)->vertices[0], (*itr_triangle)->vertices[1]);
		}
		if((*itr_triangle)->triangles[1] == NULL)
		{
			AddEdge((*itr_triangle)->vertices[1], (*itr_triangle)->vertices[2]);
		}
		if((*itr_triangle)->triangles[2] == NULL)
		{
			AddEdge((*itr_triangle)->vertices[2], (*itr_triangle)->vertices[0]);
		}
	}

	cout << "# Finished" << endl;
	cout << "# Number of Edges = " << edges.size() << endl;
	
	list<EDGE*>::iterator itr_edge = edges.begin();

	while(itr_edge != edges.end())
	{
		cout << "# Number of Edges = " << edges.size() << endl;
		HOLE* hole = new HOLE;
		holes.push_back(hole);
		hole->AddEdge((*itr_edge));
		VERTEX* v0 = (*itr_edge)->vertices[0];
		edges.erase(itr_edge);
		cout << "# Number of Edges = " << edges.size() << endl;
		list<EDGE*>::iterator itr_edge = edges.begin();
		while(itr_edge != edges.end())
		{
			if((*itr_edge)->vertices[0] == v0 || (*itr_edge)->vertices[1] == v0)
			{
				hole->AddEdge((*itr_edge));
				if(v0 == (*itr_edge)->vertices[0])
				{
					v0 = (*itr_edge)->vertices[1];
				}
				else
				{
					v0 = (*itr_edge)->vertices[0];
				}
				edges.erase(itr_edge);
				itr_edge = edges.begin();
			}
			else
			{
				itr_edge++;
			}
		}
		itr_edge = edges.begin();
	}
	cout << "# Finished" << endl;
}

void TRIANGULAR_SURFACE::ChkBoundaryVertices(T* size)
{
	TRAVERSE_VERTICES
	{
		for(int i = 0; i < 3; i++)
		{
			if((*itr_vertex)->GetPosition()[i] <= (T)0 || (*itr_vertex)->GetPosition()[i] >= (T)size[i])
			{
				(*itr_vertex)->is_boundary = true;
			}
		}
	}
}

void TRIANGULAR_SURFACE::ChkTrianglesNeighborConnectivity()
{
	{
		TRAVERSE_TRIANGLES
		{
			TRIANGLE* triangle = *itr_triangle;
			triangle->triangles[0] = NULL;
			triangle->triangles[1] = NULL;
			triangle->triangles[2] = NULL;
		}
	}
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;
		triangle->ChkNeighborConnectivity();
	}
}

VERTEX* TRIANGULAR_SURFACE::ChkNearestTextureCoordinate(T u, T v)
{
	VERTEX* nearest = NULL;
	T distance = (T)100000000;
	TRAVERSE_TRIANGLES
	{
		TRIANGLE* triangle = *itr_triangle;
		VERTEX** tri_vertices = triangle->vertices;
		for(int i = 0; i < 3; i++)
		{
			T uv[2] = {triangle->uv[i].x, triangle->uv[i].y};
			T dis = (uv[0] - u)*(uv[0] - u) + (uv[1] - v)*(uv[1] - v);
			if(dis < distance)
			{
				distance = dis;
				nearest = tri_vertices[i];
			}
		}
	}
	return nearest;
}

void TRIANGULAR_SURFACE::DrawEdges()
{
	//for_each(triangles.begin(), triangles.end(), mem_fun(&TRIANGLE::DrawEdges));
	for_each(edges.begin(), edges.end(), mem_fun(&EDGE::Draw));
}

void TRIANGULAR_SURFACE::DrawFaceNormals()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawNormal();
	}
}

void TRIANGULAR_SURFACE::DrawHoles()
{
	list<HOLE*>::iterator itr_hole;
	for(itr_hole = holes.begin(); itr_hole != holes.end(); itr_hole++)
	{
		(*itr_hole)->Draw();
	}
}

void TRIANGULAR_SURFACE::DrawHoles(int index)
{
	list<HOLE*>::iterator itr_hole = holes.begin();
	for(int i = 0; i < index; i++)
	{
		itr_hole++;
	}
	(*itr_hole)->Draw();
}

void TRIANGULAR_SURFACE::DrawTriangles(const bool& draw_front)
{
	if(draw_front)
	{
		TRAVERSE_TRIANGLES
		{
			(*itr_triangle)->Draw();
		}
	}
	else
	{
		TRAVERSE_TRIANGLES
		{
			(*itr_triangle)->DrawBack();
		}
	}
}

void TRIANGULAR_SURFACE::DrawVertices()
{
	glBegin(GL_POINTS);
	for(int i = 0; i < (int)vertices.size(); i++)
	{
		//glVertex2fv((float*)vertices[i]->x);
		glVertex2f((float)vertices[i]->x[0], (float)vertices[i]->x[1]);
	}
	glEnd();
}

void TRIANGULAR_SURFACE::DrawVertexNormals()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawNormal();
	}
}

void TRIANGULAR_SURFACE::DrawVertexDeviation()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawDeviation();
	}
}

void TRIANGULAR_SURFACE::DrawTrianglesNeighborConnectivity()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawNeighborConnectivity();
	}
}

void TRIANGULAR_SURFACE::DrawTrianglesCenter()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawCenter();
	}
}

void TRIANGULAR_SURFACE::DrawVerticesNeighborConnectivity()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawNeighborConnectivity();
	}
}

void TRIANGULAR_SURFACE::DrawCurvatureNormal()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawCurvatureNormal();
	}
}

void TRIANGULAR_SURFACE::DrawVertexVelocity()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawVelocity();
	}
}

// Need to study about the filtering
void TRIANGULAR_SURFACE::Filtering(T lambda)
{
	TRAVERSE_VERTICES
	{
		if((*itr_vertex)->triangles.size() == 0)
		{
			continue;
		}
		if((*itr_vertex)->is_boundary == true)
		{
			continue;
		}

		list<VERTEX*> v;
		list<VERTEX*>::iterator itr_v;
		list<TRIANGLE*>::iterator itr_triangle;

		for(itr_triangle = (*itr_vertex)->triangles.begin(); itr_triangle != (*itr_vertex)->triangles.end(); itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices[0]);
			v.remove((*itr_triangle)->vertices[1]);
			v.remove((*itr_triangle)->vertices[2]);
			v.push_back((*itr_triangle)->vertices[0]);
			v.push_back((*itr_triangle)->vertices[1]);
			v.push_back((*itr_triangle)->vertices[2]);
		}
		v.remove((*itr_vertex));

		T Li[2] = {(T)0, (T)0};
		T E = (T)0;
		for(itr_v = v.begin(); itr_v != v.end(); itr_v++)
		{
			T dist[2];
			T c = (T)1;
			ARRAY_VECTOR2::sub<T>((*itr_v)->GetPosition(), (*itr_vertex)->GetPosition(), dist);
			ARRAY_VECTOR2::det<T>(dist, &c);
			E += c;
			Li[0] += ((*itr_v)->GetPosition()[0] - (*itr_vertex)->GetPosition()[0])/c;
			Li[1] += ((*itr_v)->GetPosition()[1] - (*itr_vertex)->GetPosition()[1])/c;
		}

		Li[0] = (T)2/E*Li[0];
		Li[1] = (T)2/E*Li[1];

		(*itr_vertex)->deviation[0] -= lambda*Li[0];
		(*itr_vertex)->deviation[1] -= lambda*Li[1];

		(*itr_vertex)->GetPosition()[0] += lambda*Li[0];
		(*itr_vertex)->GetPosition()[1] += lambda*Li[1];
	}
}

void TRIANGULAR_SURFACE::Filtering2(T lambda)
{
	vector<VERTEX*>::iterator itr_vertex;
	for(itr_vertex = vertices.begin(); itr_vertex != vertices.end(); itr_vertex++)
	{
		list<VERTEX*> v;
		list<VERTEX*>::iterator itr_v;
		list<TRIANGLE*>::iterator itr_triangle;
		for(itr_triangle = (*itr_vertex)->triangles.begin(); itr_triangle != (*itr_vertex)->triangles.end(); itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices[0]);
			v.remove((*itr_triangle)->vertices[1]);
			v.remove((*itr_triangle)->vertices[2]);
			v.push_back((*itr_triangle)->vertices[0]);
			v.push_back((*itr_triangle)->vertices[1]);
			v.push_back((*itr_triangle)->vertices[2]);
		}
		v.remove((*itr_vertex));
		
		ARRAY_VECTOR2::set<T>((*itr_vertex)->s, (T)0);
		for(itr_v = v.begin(); itr_v != v.end(); itr_v++)
		{
			T dist[2];
			ARRAY_VECTOR2::sub<T>((*itr_v)->GetPosition(), (*itr_vertex)->GetPosition(), dist);
			ARRAY_VECTOR2::add<T>((*itr_vertex)->s, dist, (*itr_vertex)->s);
		}
		ARRAY_VECTOR2::div<T>((*itr_vertex)->s, (T)v.size());
		ARRAY_VECTOR2::det<T>((*itr_vertex)->s, &(*itr_vertex)->stress);
		if(((*itr_vertex)->s[0]*(*itr_vertex)->n[0] + (*itr_vertex)->s[1]*(*itr_vertex)->n[1] + (*itr_vertex)->s[2]*(*itr_vertex)->n[2]) < (T)0)
		{
			(*itr_vertex)->stress *= -(T)1;
		}
	}

	for(itr_vertex = vertices.begin(); itr_vertex != vertices.end(); itr_vertex++)
	{
		if((*itr_vertex)->triangles.size() == 0)
		{
			continue;
		}

		list<VERTEX*> v;
		list<VERTEX*>::iterator itr_v;
		list<TRIANGLE*>::iterator itr_triangle;
		for(itr_triangle = (*itr_vertex)->triangles.begin(); itr_triangle != (*itr_vertex)->triangles.end(); itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices[0]);
			v.remove((*itr_triangle)->vertices[1]);
			v.remove((*itr_triangle)->vertices[2]);
			v.push_back((*itr_triangle)->vertices[0]);
			v.push_back((*itr_triangle)->vertices[1]);
			v.push_back((*itr_triangle)->vertices[2]);
		}
		v.remove((*itr_vertex));

		T s_mean[3] = {(T)0, (T)0, (T)0};
		for(itr_v = v.begin(); itr_v != v.end(); itr_v++)
		{
			s_mean[0] += (*itr_v)->s[0];
			s_mean[1] += (*itr_v)->s[1];
			s_mean[2] += (*itr_v)->s[2];
		}
		s_mean[0] /= (T)v.size();
		s_mean[1] /= (T)v.size();
		s_mean[2] /= (T)v.size();

		T dx[3];
		ARRAY_VECTOR2::sub<T>((*itr_vertex)->s, s_mean, dx);
		(*itr_vertex)->GetPosition()[0] += lambda*dx[0];
		(*itr_vertex)->GetPosition()[1] += lambda*dx[1];
	}
}

void TRIANGULAR_SURFACE::Reset()
{
	TRAVERSE_VERTICES
	{
		DELETE_POINTER(*itr_vertex);
	}

	TRAVERSE_EDGES
	{
		DELETE_POINTER(*itr_edge);
	}

	TRAVERSE_TRIANGLES
	{
		DELETE_POINTER(*itr_triangle);
	}

	vertices.clear();
	edges.clear();
	triangles.clear();
}