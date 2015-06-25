#pragma once

#include "OPENGL_DRIVER.h"
#include <vector>

class OPENGL_LIGHT_MANAGER
{
public: // Essential Data
	GLint							max_lights;
	GLint							num_lights;

	vector<OPENGL_LIGHT_PROPERTY>	lights;

	OPENGL_DRIVER*					driver;

	OPENGL_VEC3						position;
	OPENGL_VEC3						scale;

public: // Construtor and Destructor
	OPENGL_LIGHT_MANAGER(OPENGL_DRIVER* driver_input)
		: driver(driver_input), max_lights(8), num_lights(0)
	{}

	~OPENGL_LIGHT_MANAGER(void)
	{}

public: // Initialization Function
	void InitLights()
	{
		glGetIntegerv(GL_MAX_LIGHTS, &max_lights);
	}

	void EnableLight(int light_index, bool is_enable)
	{
		if ((light_index < 0) || (light_index > max_lights))
		{
			return;
		}

		if ((light_index < 0) || (light_index >= num_lights))
		{
			return;
		}

		int light_index_full = GL_LIGHT0 + light_index;

		if (is_enable)
		{
			glEnable(light_index_full);
		}
		else
		{
			glDisable(light_index_full);
		}
	}

	int AddLight(const OPENGL_LIGHT_PROPERTY& light)
	{
		int light_index = GL_LIGHT0;
		for (light_index = GL_LIGHT0; light_index < GL_LIGHT0 + max_lights; ++light_index)
		{
			if (!glIsEnabled(light_index))
			{
				break;
			}
		}

		if (light_index == GL_LIGHT0 + max_lights)
		{
			return -1;
		}

		ApplyLight(light_index, light);

		glEnable(light_index);				// When added, enable the light

		lights.push_back(light);
		num_lights++ ;

		return light_index;
	}

	void ApplyLight(int light_index, const OPENGL_LIGHT_PROPERTY& light, bool debug_light_draw = false)
	{
		if (light_index >= GL_LIGHT0 + max_lights)
		{
			return;
		}

		GLfloat data[4];
		switch (light.light_type)
		{
		case LIGHT_SPOT:
			{
				data[0] = light.direction.x;
				data[1] = light.direction.y;
				data[2] = light.direction.z;
				data[3] = 0.0f;
				glLightfv(light_index, GL_SPOT_DIRECTION, data);

				data[0] = light.position.x;
				data[1] = light.position.y;
				data[2] = light.position.z;
				data[3] = 1.0f;

				glLightfv(light_index, GL_POSITION, data);
				glLightf(light_index, GL_SPOT_EXPONENT, light.spot_exponent);
				glLightf(light_index, GL_SPOT_CUTOFF, light.spot_cutoff);

				// For debugging
				if (debug_light_draw)
				{
					driver->SetDefaultRenderStates3DMode();
					glColor3f(0.2f, 1.0f, 0.2f);
					glPushMatrix();
						glTranslatef(light.position.x, light.position.y, light.position.z);
						glutSolidSphere(0.05, 10, 10);
					glPopMatrix();

					glBegin(GL_LINES);
						glVertex3f(light.position.x, light.position.y, light.position.z);
						glVertex3f(light.position.x + light.direction.x, light.position.y + light.direction.y, light.position.z + light.direction.z);
					glEnd();
				}
			}
			break;
			
		case LIGHT_POINT:
			{
				data[0] = light.position.x;
				data[1] = light.position.y;
				data[2] = light.position.z;
				data[3] = 1.0f;						// 1.0f for point light
				glLightfv(light_index, GL_POSITION, data);
				glLightf(light_index, GL_SPOT_EXPONENT, 0.0f);
				glLightf(light_index, GL_SPOT_CUTOFF, 180.0f);

				// For debugging
				if (debug_light_draw)
				{
					driver->SetDefaultRenderStates3DMode();
					glColor3f(0.2f, 0.2f, 1.0f);
					glPushMatrix();
						glTranslatef(light.position.x, light.position.y, light.position.z);
						glutSolidSphere(0.05, 10, 10);
					glPopMatrix();
				}
			}
			break;

		case LIGHT_DIRECTIONAL:
			{
				data[0] = -light.direction.x;
				data[1] = -light.direction.y;
				data[2] = -light.direction.z;
				data[3] = 0.0f;						// 0.0f for directional light
				glLightfv(light_index, GL_POSITION, data);
				glLightf(light_index, GL_SPOT_EXPONENT, 0.0f);
				glLightf(light_index, GL_SPOT_CUTOFF, 180.0f);

				// For debugging
				if (debug_light_draw)
				{
					driver->SetDefaultRenderStates3DMode();
					glColor3f(1.0f, 0.2f, 0.2f);
					glutSolidSphere(0.05, 10, 10);
					glBegin(GL_LINES);
						glVertex3f(0.0f, 0.0f, 0.0f);
						glVertex3f(light.direction.x, light.direction.y, light.direction.z);
					glEnd();
				}
			}
			break;
		}

		// Set Diffuse color
		data[0] = light.diffuse_color.r;
		data[1] = light.diffuse_color.g;
		data[2] = light.diffuse_color.b;
		data[3] = light.diffuse_color.a;
		glLightfv(light_index, GL_DIFFUSE, data);

		// Set Specular color
		data[0] = light.specular_color.r;
		data[1] = light.specular_color.g;
		data[2] = light.specular_color.b;
		data[3] = light.specular_color.a;
		glLightfv(light_index, GL_SPECULAR, data);

		// Set Ambient color
		data[0] = light.ambient_color.r;
		data[1] = light.ambient_color.g;
		data[2] = light.ambient_color.b;
		data[3] = light.ambient_color.a;
		glLightfv(light_index, GL_AMBIENT, data);

		// Set Attenuation
		glLightf(light_index, GL_CONSTANT_ATTENUATION, light.constant_attenuation);
		glLightf(light_index, GL_LINEAR_ATTENUATION, light.linear_attenuation);
		glLightf(light_index, GL_QUADRATIC_ATTENUATION, light.quadratic_attenuation);
	}

	void DeleteAllLight()
	{
		for (int i = 0; i < max_lights; ++i)
		{
			glDisable(GL_LIGHT0 + i);
		}

		vector<OPENGL_LIGHT_PROPERTY>().swap(lights);
		num_lights = 0;
	}

	void UpdateLightPosition(bool debug_light_draw = false, bool enable_light = false)
	{
		for (unsigned int i = 0; i < lights.size(); ++i)
		{
			ApplyLight(GL_LIGHT0 + i, lights[i], debug_light_draw);
			if (enable_light)
			{
				glEnable(GL_LIGHT0 + i);
			}
		}
	}

	void ToggleEnableLight(int light_index)
	{
		if ((light_index < 0) || (light_index > max_lights))
		{
			return;
		}

		if ((light_index < 0) || (light_index >= num_lights))
		{
			return;
		}

		int light_index_full = GL_LIGHT0 + light_index;
		if (glIsEnabled(light_index_full))
		{
			glDisable(light_index_full);
			cout << "Disabled Light = " << light_index << endl;
		}
		else
		{
			glEnable(light_index_full);
			cout << "Enabled Light = " << light_index << endl;
		}
	}

	void LoadOpenGLSettings(SCRIPT_BLOCK* script_block, GRID_STRUCTURE_2D& grid, bool default_lights = false)
	{
		DeleteAllLight();

		if (!script_block)
		{
			return;
		}

		if (default_lights)
		{
			OPENGL_LIGHT_PROPERTY light_property1;
			light_property1.light_type = LIGHT_POINT;

			light_property1.ambient_color.Set(0.1f, 0.1f, 0.1f);
			light_property1.diffuse_color.Set(1.0f, 1.0f, 1.0f);
			light_property1.specular_color.Set(1.0f, 1.0f, 1.0f);

			light_property1.constant_attenuation = 2.4;
			light_property1.linear_attenuation = 0.0;
			light_property1.quadratic_attenuation = 0.0;

			light_property1.position.Set(grid.x_max, grid.y_max, 0);
			AddLight(light_property1);

			OPENGL_LIGHT_PROPERTY light_property2;
			light_property2.light_type = LIGHT_POINT;
			
			light_property2.ambient_color.Set(0.1f, 0.1f, 0.1f);
			light_property2.diffuse_color.Set(1.0f, 1.0f, 1.0f);
			light_property2.specular_color.Set(1.0f, 1.0f, 1.0f);

			light_property2.constant_attenuation = 2.4;
			light_property2.linear_attenuation = 0.0;
			light_property2.quadratic_attenuation = 0.0;

			light_property2.position.Set(grid.x_min, grid.y_max, 0);
			AddLight(light_property2);
		}

		vector<SCRIPT_BLOCK*> light_list = script_block->GetBlockList("OPENGL_LIGHT_OBJECT");
		for (int i = 0; i < light_list.size(); i++)
		{
			OPENGL_LIGHT_PROPERTY light_property;
			// Light_type
			const char* light_type = light_list[i]->GetString("light_type", "POINT");
			if (boost::iequals(light_type, "POINT"))
			{
				light_property.light_type = LIGHT_POINT;
			}
			else if (boost::iequals(light_type, "DIRECTIONAL"))
			{
				light_property.light_type = LIGHT_DIRECTIONAL;
			}
			else if (boost::iequals(light_type, "SPOT"))
			{
				light_property.light_type = LIGHT_SPOT;
			}

			// Position
			VECTOR_3D<T> pos = light_list[i]->GetVector3("position", VECTOR_3D<T>((T)0.0, (T)0.0, (T)0.0));
			light_property.position.Set(pos.x, pos.y, pos.z);

			// Direction
			VECTOR_3D<T> dir = light_list[i]->GetVector3("direction", VECTOR_3D<T>((T)0.0, (T)0.0, (T)0.0));
			light_property.direction.Set(dir.x, dir.y, dir.z);

			// Color
			VECTOR_3D<T> ambient = light_list[i]->GetVector3("ambient_color", VECTOR_3D<T>((T)0.0, (T)0.0, (T)0.0));
			light_property.ambient_color.Set((GLfloat) ambient.x, (GLfloat) ambient.y, (GLfloat) ambient.z);
			VECTOR_3D<T> diffuse = light_list[i]->GetVector3("diffuse_color", VECTOR_3D<T>((T)0.0, (T)0.0, (T)0.0));
			light_property.diffuse_color.Set((GLfloat) diffuse.x, (GLfloat) diffuse.y, (GLfloat) diffuse.z);
			VECTOR_3D<T> specular = light_list[i]->GetVector3("specular_color", VECTOR_3D<T>((T)0.0, (T)0.0, (T)0.0));
			light_property.specular_color.Set((GLfloat) specular.x, (GLfloat) specular.y, (GLfloat) specular.z);
			
			// Attenuation
			light_property.constant_attenuation = light_list[i]->GetFloat("constant_attenuation", (T)1.0f);
			light_property.linear_attenuation = light_list[i]->GetFloat("linear_attenuation", (T)0.0f);
			light_property.quadratic_attenuation = light_list[i]->GetFloat("quadratic_attenuation", (T)0.0f);

			// Spot
			light_property.spot_exponent = light_list[i]->GetFloat("spot_exponent", (T)2.0f);
			light_property.spot_cutoff = light_list[i]->GetFloat("spot_cutoff", (T)45.0f);

			AddLight(light_property);
		}
	}
};
