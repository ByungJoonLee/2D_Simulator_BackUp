#pragma once

#include "nvGlutManipulators.h"
#include "GRID_STRUCTURE_2D.h"
#include <Windows.h>
#include <GL/glut.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <io.h>
#include <boost/filesystem.hpp>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>

#define SCALE_PAN_FACTOR		0.001f
#define SCALE_ZOOM_FACTOR		0.0015f
#define SELECT_BUFFER_LENGTH	1024

class OPENGL_GLUT_CAMERA
{
public: // Enumerates
	enum CAMERA_MODE
	{
		MANIPULATION_MODE = 0,
		FILE_MODE,
	};

public: // Essential Data
	nv::GlutExamine			trackball;
	CAMERA_MODE				camera_mode;
	std::string				camera_save_dir;
	float					camera_loaded_matrix[16];
	int						loaded_camera_num;
	double					angle_of_view;
	double					z_near;
	double					z_far;
	int						screen_width;
	int						screen_height;
	float					focus_distance;
	float					scale_pan;
	float					scale_zoom;

	GLuint*					select_buffer;
	int						select_id;

	bool					auto_camera_load;
	int						camera_num;

public: // Constructor and Destructor
	OPENGL_GLUT_CAMERA(void)
		: camera_mode(MANIPULATION_MODE)
		, loaded_camera_num(0)
		, angle_of_view(40.0)
		, z_near(0.1)
		, z_far(1000.0)
		, screen_width(640)
		, screen_height(480)
		, focus_distance(1.0)
		, scale_pan(SCALE_PAN_FACTOR)
		, scale_zoom(SCALE_ZOOM_FACTOR)
		, select_buffer(0)
		, select_id(-1)
		, auto_camera_load(false)
		, camera_num(0)
	{
		trackball.setDollyActivate(GLUT_LEFT_BUTTON, GLUT_ACTIVE_CTRL);
		trackball.setPanActivate(GLUT_LEFT_BUTTON, GLUT_ACTIVE_SHIFT);
		trackball.setTrackballActivate(GLUT_LEFT_BUTTON);

		trackball.setTrackballScale(0.7f);
		trackball.setDollyScale(scale_zoom);
		trackball.setPanScale(scale_pan);
		trackball._examine._pan.z = -5.0f;

		for (int i = 0; i < 16; ++i)
		{
			camera_loaded_matrix[i] = 0.0f;
		}
		camera_loaded_matrix[0] = camera_loaded_matrix[5] = camera_loaded_matrix[10] = camera_loaded_matrix[15] = 1.0f;

		camera_save_dir = GetApplicationDir();
	}

	~OPENGL_GLUT_CAMERA(void)
	{}

public: // Initialization Function
	void Initialize(std::string& script_abs_path, GRID_STRUCTURE_2D& grid)
	{
		SCRIPT_READER script_reader(script_abs_path.c_str());
		SCRIPT_BLOCK script_block = script_reader.FindBlock("SIMULATION_WORLD");

		auto_camera_load = script_block.GetBoolean("auto_camera_load", false);
		camera_num = script_block.GetInteger("camera_num", 0);

		SetCameraNum(camera_num);
		CalcDollyPosition(grid);
	}

public: // Member Functions
	void GlutDisplay()
	{
		glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			trackball.applyTransform();
	}

	void GlutMouse(int button, int state, int x, int y)
	{
		trackball.mouse(button, state, x, y);
	}

	void GlutMotion(int x, int y)
	{
		// Reset Scale
		if (trackball._examine._dActive)
		{
			trackball.setDollyScale(scale_zoom - (SCALE_ZOOM_FACTOR*std::abs(trackball._examine._dolly.z)));
		}
		if (trackball._examine._pActive)
		{
			trackball.setPanScale(scale_pan - (SCALE_PAN_FACTOR*std::abs(trackball._examine._dolly.z)));
		}

		trackball.motion(x, y);

		// Set to manipulation mode
		camera_mode = MANIPULATION_MODE;
	}

	void GlutReshape(int width, int height)
	{
		screen_width = width;
		screen_height = height;

		LoadProjectMatrix();

		glViewport(0, 0, screen_width, screen_height);

		trackball.reshape(width, height);
	}

	void SetCameraMode(CAMERA_MODE mode)
	{
		camera_mode = mode;
	}

	CAMERA_MODE GetCameraMode()
	{
		return camera_mode;
	}

	void SetCameraNum(int num)
	{
		camera_num = num;
	}

	int GetCameraNum()
	{
		return loaded_camera_num - 1;
	}

	void SetCameraFileDir(const char* dir)
	{
		if (dir)
		{
			camera_save_dir = dir;
		}
	}

	void SaveCamera()
	{
		cout << "Save Camera....";

		char buffer[MAX_PATH] = {0, };
		GetCurrentDirectory(MAX_PATH, buffer);
		
		SetCurrentDirectory(camera_save_dir.c_str());			//.cam file is saved under online_simulator.exe
		static int count = 0;
		char filename[MAX_PATH] = {0, };
		
		while (1)
		{
			sprintf_s(filename, "./cam_%d.cam2", count);
			if (_access(filename, 0) == -1)		// file not exist
			{
				break;
			}
			count++;
		}

		ofstream file(filename, std::ios::out | std::ios::trunc);

		if (!file)
		{
			cerr << "can't open file \"" << "./cam.cam2" << "\"" << endl;
			return;
		}

		file << trackball._examine._pan.x << endl;
		file << trackball._examine._pan.y << endl;
		file << trackball._examine._pan.z << endl;

		file << trackball._examine._dolly.x << endl;
		file << trackball._examine._dolly.y << endl;
		file << trackball._examine._dolly.z << endl;

		file << trackball._examine._centroid.x << endl;
		file << trackball._examine._centroid.y << endl;
		file << trackball._examine._centroid.z << endl;

		file << trackball._examine._r.x << endl;
		file << trackball._examine._r.y << endl;
		file << trackball._examine._r.z << endl;
		file << trackball._examine._r.w << endl;

		cout << filename << " success" << endl;

		file.close();

		loaded_camera_num = count + 1;
		camera_mode = FILE_MODE;

		SetCurrentDirectory(buffer);
	}

	void LoadCamera()
	{
		char cam_file_name[MAX_PATH] = {0, };
		sprintf_s(cam_file_name, "&scam_%d.cam2", camera_save_dir.c_str(), loaded_camera_num);

		char buffer[MAX_PATH] = {0, };
		GetCurrentDirectory(MAX_PATH, buffer);

		SetCurrentDirectory(camera_save_dir.c_str());			//	.cam file is loaded under online_simulator.exe
		ifstream file(cam_file_name, std::ios::in);

		if (!file.is_open())
		{
			file.close();
			SetCurrentDirectory(buffer);

			if (loaded_camera_num == 0)
			{
				cerr << cam_file_name << " is not loaded correctly." << endl;
				camera_mode = MANIPULATION_MODE;
			}
			else
			{
				loaded_camera_num = 0;
				LoadCamera();
			}
			return;
		}

		file >> trackball._examine._pan.x;
		file >> trackball._examine._pan.y;
		file >> trackball._examine._pan.z;

		file >> trackball._examine._dolly.x;
		file >> trackball._examine._dolly.y;
		file >> trackball._examine._dolly.z;

		file >> trackball._examine._centroid.x;
		file >> trackball._examine._centroid.y;
		file >> trackball._examine._centroid.z;

		file >> trackball._examine._r.x;
		file >> trackball._examine._r.y;
		file >> trackball._examine._r.z;
		file >> trackball._examine._r.w;

		file.close();

		camera_mode = FILE_MODE;
		cout << cam_file_name << " is loaded successfully." << endl;
		loaded_camera_num++;

		SetCurrentDirectory(buffer);
	}

	void AutoLoadCamera()
	{
		if (auto_camera_load)
		{
			loaded_camera_num = 0;
			LoadCamera();
		}
	}

	void BeginSelect(int x, int y)
	{
		if (select_buffer)
		{
			delete [] select_buffer;
			select_buffer = 0;
		}

		select_buffer = new GLuint[SELECT_BUFFER_LENGTH];
		glSelectBuffer(SELECT_BUFFER_LENGTH, select_buffer);

		static GLint viewport[4];
		viewport[0] = 0;
		viewport[1] = screen_height;
		viewport[2] = screen_width;
		viewport[3] = -screen_height;

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
			glLoadIdentity();
			glRenderMode(GL_SELECT);
			gluPickMatrix(x, y, 2, 2, viewport);
			LoadProjectMatrix(true);
			glInitNames();
	}

	void EndSelect()
	{
		glFlush();

		GLint nbHits = glRenderMode(GL_RENDER);

		if (nbHits <= 0)
		{
			select_id = -1;
		}
		else
		{
			GLuint zMin = select_buffer[1];
			select_id = select_buffer[3];
			for (int i = 1; i < nbHits; ++i)
			{
				if (select_buffer[4*i + 1] < zMin)
				{
					zMin = select_buffer[4*i + 1];
					select_id = select_buffer[4*i + 3];
				}
			}
		}

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();

		if (select_buffer)
		{
			delete [] select_buffer;
			select_buffer = 0;
		}

		glMatrixMode(GL_MODELVIEW);
	}

	int GetSelectedId()
	{
		return select_id;
	}

	void UpAngle()
	{
		angle_of_view += 1.0;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(angle_of_view, (double)screen_width/screen_height, z_near, z_far);
		glMatrixMode(GL_MODELVIEW);
	}

	void DownAngle()
	{
		angle_of_view -= 1.0;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(angle_of_view, (double)screen_width/screen_height, z_near, z_far);
		glMatrixMode(GL_MODELVIEW);
	}

	void SavePBRTCamera()
	{
		GLfloat mv[16] = {0.0, };
		glGetFloatv(GL_MODELVIEW_MATRIX, mv);

		float camera_pos_x = -(mv[12]*mv[0] + mv[13]*mv[1] + mv[14]*mv[2]);
		float camera_pos_y = -(mv[12]*mv[4] + mv[13]*mv[5] + mv[14]*mv[6]);
		float camera_pos_z = -(mv[12]*mv[8] + mv[13]*mv[9] + mv[14]*mv[10]);

		float camera_target_x = camera_pos_x - mv[2];
		float camera_target_y = camera_pos_y - mv[6];
		float camera_target_z = camera_pos_z - mv[10];

		float camera_up_x = mv[1];
		float camera_up_y = mv[5];
		float camera_up_z = mv[9];

		ofstream file("d:/pbrtcamera.txt");
		file << "LookAt " << camera_pos_x << " " << camera_pos_y << " " << camera_pos_z << endl;
		file << camera_target_x << " " << camera_target_y << " " << camera_target_z << endl;
		file << camera_up_x << " " << camera_up_y << " " << camera_up_z << endl;

		file << "[" << angle_of_view << "]" << endl;
		file.close();
	}

	void LoadProjectMatrix(bool for_select = false)
	{
		glMatrixMode(GL_PROJECTION);
		if (!for_select)
		{
			glLoadIdentity();
		}
		gluPerspective(angle_of_view, (double)screen_width/screen_height, z_near, z_far);
	}

	void CalcDollyPosition(GRID_STRUCTURE_2D& grid)
	{
		float diff_x = (grid.x_max - grid.x_min);
		float diff_y = (grid.y_max - grid.y_min);
		//float diff_z = (grid.z_max - grid.z_min);
		float diff_max = (std::max)(diff_x, diff_y);
		//diff_max = (std::max)(diff_max, diff_z);
		focus_distance = diff_max/tan((angle_of_view)*M_PI_2/180.0);

		float cen_x = (grid.x_max + grid.x_min)/2.0f;
		float cen_y = (grid.y_max + grid.y_min)/2.0f;
		//float cen_z = (grid.z_max + grid.z_min)/2.0f;
		/*nv::vec3f cv(cen_x, cen_y, cen_z);
		trackball.setCenterOfRotation(cv);

		trackball._examine._pan.x = -cen_x;
		trackball._examine._pan.y = -cen_y;
		trackball._examine._pan.z = -focus_distance;

		scale_pan = focus_distance*SCALE_PAN_FACTOR;
		scale_zoom = focus_distance*SCALE_ZOOM_FACTOR;*/
	}

	static std::string GetApplicationDir()
	{
		char buffer[MAX_PATH] = {0, };
		::GetModuleFileNameA(NULL, buffer, MAX_PATH);
		boost::filesystem::path p(buffer);
		return p.parent_path().string() + std::string("\\");
	}
};

		
		
