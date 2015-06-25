#pragma once

#include "OPENGL_COMMON.h"
#include "COMMON_DEFINITION.h"
#include "MARCHING_SQUARES_ALGORITHM.h"
#include <Cg/cg.h>
#include <math.h>

class OPENGL_DRIVER
{
public: // Essential Data
	CGcontext				cg_context;
	int						selected_name_id;

	static OPENGL_MATERIAL	material_2Dmode;
	static OPENGL_MATERIAL	material_3Dmode;
	static OPENGL_MATERIAL	material_line;
	
public: // Constructor and Destructor
	OPENGL_DRIVER(void)
	{
		// Default material
		material_3Dmode.color_material = COLOR_MATERIAL_DIFFUSE_AND_AMBIENT;

		material_2Dmode.is_lighting = false;
		material_2Dmode.z_buffer = COMPARISON_NEVER;
		material_2Dmode.color_material = COLOR_MATERIAL_DIFFUSE_AND_AMBIENT;

		material_line.is_lighting = false;
		material_line.color_material = COLOR_MATERIAL_DIFFUSE_AND_AMBIENT;

		cg_context = cgCreateContext();

		SetSelectedNameId(-1);
	}

	~OPENGL_DRIVER(void)
	{}

public: // Initialization Function
	void InitGenericGL()
	{
		// Lighting
		GLfloat data[4] = {0.0f, 0.0f, 0.0f, 0.0f,};
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, data);
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);

		// Anti-alising
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);

		// ETC
		glClearDepth(1.0);

		glDepthFunc(GL_LEQUAL);

		glFrontFace(GL_CW);

		glAlphaFunc(GL_GREATER, 0.0f);

		SetAntialising(ANTIALISING_OFF);
	}

public: // Member Functions
	void SetRenderStatesByMaterial(const OPENGL_MATERIAL& material)
	{
		// Color material
		switch(material.color_material)
		{
		case COLOR_MATERIAL_NONE:
			glDisable(GL_COLOR_MATERIAL);
			break;
		case COLOR_MATERIAL_DIFFUSE:
			glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
			break;
		case COLOR_MATERIAL_AMBIENT:
			glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);
			break;
		case COLOR_MATERIAL_EMISSIVE:
			glColorMaterial(GL_FRONT_AND_BACK, GL_EMISSION);
			break;
		case COLOR_MATERIAL_SPECULAR:
			glColorMaterial(GL_FRONT_AND_BACK, GL_SPECULAR);
			break;
		case COLOR_MATERIAL_DIFFUSE_AND_AMBIENT:
			glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
			break;
		};

		if (material.color_material != COLOR_MATERIAL_NONE)
		{
			glEnable(GL_COLOR_MATERIAL);
		}

		GLfloat color[4];
		if ((material.color_material != COLOR_MATERIAL_AMBIENT) && (material.color_material != COLOR_MATERIAL_DIFFUSE_AND_AMBIENT))
		{
			color[0] = material.ambient_color.GetRed();
			color[1] = material.ambient_color.GetGreen();
			color[2] = material.ambient_color.GetBlue();
			color[3] = material.ambient_color.GetAlpha();
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
		}

		if((material.color_material != COLOR_MATERIAL_DIFFUSE) && (material.color_material != COLOR_MATERIAL_DIFFUSE_AND_AMBIENT))
		{
			color[0] = material.diffuse_color.GetRed();
			color[1] = material.diffuse_color.GetGreen();
			color[2] = material.diffuse_color.GetBlue();
			color[3] = material.diffuse_color.GetAlpha();
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
		}

		if(material.color_material != COLOR_MATERIAL_EMISSIVE)
		{
			color[0] = material.emissive_color.GetRed();
			color[1] = material.emissive_color.GetGreen();
			color[2] = material.emissive_color.GetBlue();
			color[3] = material.emissive_color.GetAlpha();
			glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, color);
		}

		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.shininess);

		GLfloat color2[4] = {0.f, 0.f, 0.f, 1.f};
		if ((material.shininess != 0.0f) && (material.color_material != COLOR_MATERIAL_SPECULAR))
		{
			color2[0] = material.specular_color.GetRed();
			color2[1] = material.specular_color.GetGreen();
			color2[2] = material.specular_color.GetBlue();
			color2[3] = material.specular_color.GetAlpha();
		}

		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color2);

		// Polygon Mode
		switch(material.polygon_mode)
		{
			case POLYGON_SOLID:
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				break;
			case POLYGON_WIREFRAME:
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				break;
			case POLYGON_POINT:
				glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
				break;
		}

		// Shade Mode
		if (material.is_gouraud_shading)
		{
			glShadeModel(GL_SMOOTH);
		}
		else
		{
			glShadeModel(GL_FLAT);
		}

		// Lighting
		if (material.is_lighting)
		{
			glEnable(GL_LIGHTING);
		}
		else
		{
			glDisable(GL_LIGHTING);
		}

		// z buffer
		switch(material.z_buffer)
		{
		case COMPARISON_NEVER:
			glDisable(GL_DEPTH_TEST);
			break;
		case COMPARISON_LESSEQUAL:
			glDepthFunc(GL_LEQUAL);
			break;
		case COMPARISON_EQUAL:
			glDepthFunc(GL_EQUAL);
			break;
		case COMPARISON_LESS:
			glDepthFunc(GL_LESS);
			break;
		case COMPARISON_NOTEQUAL:
			glDepthFunc(GL_NOTEQUAL);
			break;
		case COMPARISON_GREATEREQUAL:
			glDepthFunc(GL_GREATER);
			break;
		case COMPARISON_ALWAYS:
			glDepthFunc(GL_ALWAYS);
			break;
		}

		if (material.z_buffer != COMPARISON_NEVER)
		{
			glEnable(GL_DEPTH_TEST);
		}

		// z write
		if(material.is_z_writable)	
		{
			glDepthMask(GL_TRUE);
		}
		else
		{
			glDepthMask(GL_FALSE);
		}

		// Face culling
		if ((material.is_frontface_culling) && (material.is_backface_culling))
		{
			glCullFace(GL_FRONT_AND_BACK);
			glEnable(GL_CULL_FACE);
		}
		else
		{
			if (material.is_backface_culling)
			{
				glCullFace(GL_BACK);
				glEnable(GL_CULL_FACE);
			}
			else
			{
				if(material.is_frontface_culling)
				{
					glCullFace(GL_FRONT);
					glEnable(GL_CULL_FACE);
				}
				else
				{
					glDisable(GL_CULL_FACE);
				}
			}
		}

		// Normalization 
		if (material.is_normalize_normals)
		{
			glEnable(GL_NORMALIZE);
		}
		else
		{
			glDisable(GL_NORMALIZE);
		}

		// Thickness
		glPointSize(material.thickness);
		glLineWidth(material.thickness);
	}

	void SetDefaultRenderStatesLineMode(OPENGL_COLOR line_color = OPENGL_COLOR(0.0f, 0.0f, 0.0f), GLfloat thickness_input = 1.0f)
	{
		GLfloat temp = material_line.thickness;
		material_line.thickness = thickness_input;
		SetRenderStatesByMaterial(material_line);
		material_line.thickness = temp;
		glColor4f(line_color.GetRed(), line_color.GetGreen(), line_color.GetBlue(), line_color.GetAlpha());
	}

	void SetDefaultRenderStates2DMode()
	{
		SetRenderStatesByMaterial(material_2Dmode);
		glClear(GL_DEPTH_BUFFER_BIT);
	}

	void SetDefaultRenderStates3DMode()
	{
		SetRenderStatesByMaterial(material_3Dmode);
	}

	// Need to review after you study OPEN_GL
	void SetAntialising(ANTIALISING_MODE anti_aliasing)
	{
		// Anti aliasing
		if (GLEW_ARB_multisample && (anti_aliasing && ANTIALISING_MULTISAMPLE))
		{
			if (anti_aliasing & ANTIALISING_ALPHA_TO_COVERAGE)
			{
				glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE_ARB);
			}
			else
			{
				glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE_ARB);
			}

			glEnable(GL_MULTISAMPLE_ARB);
			if (GLEW_NV_multisample_filter_hint)
			{
				glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
			}
		}
		else
		{
			glDisable(GL_MULTISAMPLE_ARB);
		}

		if (anti_aliasing & ANTIALISING_POINT_SMOOTH)
		{
			glEnable(GL_POINT_SMOOTH);
		}
		else
		{
			glDisable(GL_POINT_SMOOTH);
		}

		if (anti_aliasing & ANTIALISING_LINE_SMOOTH)
		{
			glEnable(GL_LINE_SMOOTH);
		}
		else
		{
			glDisable(GL_LINE_SMOOTH);
		}
	}

	void DrawWireTriangles(list<TRIANGLE*>& triangles)
	{
		list<TRIANGLE*>::iterator it;
		for (it = triangles.begin(); it != triangles.end(); it++)
		{
			glBegin(GL_LINES);
			{
				#ifdef USE_FLOAT_T
					glVertex2fv((*it)->vertices[1]->GetPosition());
					glVertex2fv((*it)->vertices[2]->GetPosition());

					glVertex2fv((*it)->vertices[2]->GetPosition());
					glVertex2fv((*it)->vertices[0]->GetPosition());

					glVertex2fv((*it)->vertices[0]->GetPosition());
					glVertex2fv((*it)->vertices[1]->GetPosition());
				#else
					glVertex2dv((*it)->vertices[1]->GetPosition());
					glVertex2dv((*it)->vertices[2]->GetPosition());

					glVertex2dv((*it)->vertices[2]->GetPosition());
					glVertex2dv((*it)->vertices[0]->GetPosition());

					glVertex2dv((*it)->vertices[0]->GetPosition());
					glVertex2dv((*it)->vertices[1]->GetPosition());
				#endif
			}
			glEnd();
		}
	}

	void DrawGradientBackground(const GLfloat* upperColor, const GLfloat* lowerColor)
	{
		GLboolean isDepthTest, isLighting;
		glGetBooleanv(GL_DEPTH_TEST, &isDepthTest);
		glGetBooleanv(GL_LIGHTING, &isLighting);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);

		// glClear(GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
			glLoadIdentity();

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
				glLoadIdentity();

				glBegin(GL_QUADS);
					// Lower color
					glColor3fv(lowerColor);
					glVertex2f(-1.0, -1.0);
					glVertex2f(1.0, -1.0);
					
					// Upper color
					glColor3fv(upperColor);
					glVertex2f(1.0, 1.0);
					glVertex2f(-1.0, 1.0);
				glEnd();
	
			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();

		if (isDepthTest)
		{
			glEnable(GL_DEPTH_TEST);
		}				
		if (isLighting)
		{
			glEnable(GL_LIGHTING);
		}	
	}

	void DrawGradientBackground(GLfloat upperR, GLfloat upperG, GLfloat upperB, GLfloat lowerR, GLfloat lowerG, GLfloat lowerB)
	{
		GLfloat upper[3] = {upperR, upperG, upperB};
		GLfloat lower[3] = {lowerR, lowerG, lowerB};

		DrawGradientBackground(&upper[0], &lower[0]);
	}
	
	void DrawGradientBackground(const GLubyte* upperColor, const GLubyte* lowerColor)
	{
		GLboolean isDepthTest, isLighting;
		glGetBooleanv(GL_DEPTH_TEST, &isDepthTest);
		glGetBooleanv(GL_LIGHTING, &isLighting);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);

		// glClear(GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
			glLoadIdentity();

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
				glLoadIdentity();

				glBegin(GL_QUADS);
					// Lower color
					glColor3ubv(lowerColor);
					glVertex2f(-1.0, -1.0);
					glVertex2f(1.0, -1.0);
					
					// Upper color
					glColor3ubv(upperColor);
					glVertex2f(1.0, 1.0);
					glVertex2f(-1.0, 1.0);
				glEnd();
	
			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();

		if (isDepthTest)
		{
			glEnable(GL_DEPTH_TEST);
		}				
		if (isLighting)
		{
			glEnable(GL_LIGHTING);
		}	
	}

	void DrawGradientBackground(GLubyte upperR, GLubyte upperG, GLubyte upperB, GLubyte lowerR, GLubyte lowerG, GLubyte lowerB)
	{
		GLubyte upper[3] = {upperR, upperG, upperB};
		GLubyte lower[3] = {lowerR, lowerG, lowerB};

		DrawGradientBackground(&upper[0], &lower[0]);
	}

	void DrawWireBox(GLfloat dx, GLfloat dy)
	{
		glBegin(GL_LINE_LOOP);
			glVertex2f(-dx, -dy);
			glVertex2f(-dx,  dy);
			glVertex2f( dx,  dy);
			glVertex2f( dx, -dy);
		glEnd();
	}

	void DrawAxis(GLfloat x_center, GLfloat y_center, GLfloat dx, GLfloat dy)
	{
		glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 0.0f);
			glVertex2f(x_center, -1.1*dy);
			glVertex2f(x_center,  1.1*dy);
		glEnd();

		glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 0.0f);	
			glVertex2f(-1.1*dx , y_center);
			glVertex2f( 1.1*dx , y_center);
		glEnd();

		glDisable(GL_LIGHTING);
		glColor3f(0.0f, 0.0f, 0.0f);
		glRasterPos3f(1.2f + 0.01f, 0.0f, 0.0f);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'x');

		glColor3f(0.0f, 0.0f, 0.0f);
		glRasterPos3f(0.0f, 1.5f + 0.07f, 0.0f);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'y');
	}

	CGcontext GetCGContext()
	{
		return cg_context;
	}

	void SetSelectedNameId(int name_id)
	{
		selected_name_id = name_id;
	}

	int GetSelectedNameId()
	{
		return selected_name_id;
	}
};