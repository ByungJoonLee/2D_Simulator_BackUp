#pragma once

#include "PROJECT_INFO.h"
#include "OPENGL_WORLD.h"
#include "SIMULATION_WORLD.h"
#include <string>
#include <boost/chrono.hpp>
#include <iomanip>
#include <deque>
#include <numeric>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

struct DISPLAY_INFO;
class SIMULATION_MANAGER;

struct DISPLAY_INFO
{
	string camera_mode;
	bool is_capture_image;
	bool is_capture_move_at_last_frame;
};

class SIMULATION_MANAGER
{
public: // Enumerates
	enum SIMULATION_STATE
	{
		NONE_STATE = 0,
		INIT_STATE,
		RUN_STATE,
		PAUSE_STATE,
	};

public: // Essential Data
	SIMULATION_STATE				state;

	SIMULATION_WORLD*				simulation;
	OPENGL_WORLD*					opengl;
	
	DISPLAY_INFO					dp_info;
	boost::chrono::duration<double> simulation_time;
	boost::chrono::duration<double> total_time;

	double							opengl_fps;

public: // Constructor and Destructor
	SIMULATION_MANAGER(void)
		: state(NONE_STATE),
		simulation(0),
		opengl(0)
	{
		dp_info.camera_mode = "Camera";
		dp_info.is_capture_image = false;
		dp_info.is_capture_move_at_last_frame = false;
		total_time.zero();
        simulation_time.zero();
	}

	~SIMULATION_MANAGER(void)
	{
		DeleteWorlds();
	}

public: // Initialization Function
	bool Initialize(string& script_path)
	{
		if (false == AnalysisScriptPath(script_path))
		{
			return false;
		}

		PROJECT_INFO::app_abs_dir = GetApplicationDir();

		CreateWorlds();

		simulation->script_abs_dir = PROJECT_INFO::script_abs_dir;
		simulation->app_abs_dir = PROJECT_INFO::app_abs_dir;
		simulation->script_base_name = PROJECT_INFO::script_basename;

		ResetSimulation();

		return true;
	}
	
	static std::string GetApplicationDir()
	{
		char buffer[MAX_PATH] = {0, };
		::GetModuleFileNameA(NULL, buffer, MAX_PATH);
		boost::filesystem::path p(buffer);

		return p.parent_path().string() + std::string("\\");
	}

	// This gives us the name of give folder
	bool AnalysisScriptPath(string& script_path)
	{
		boost::filesystem::path p(script_path);
		if (p.is_relative())
		{
			p = boost::filesystem::absolute(p);
		}

		if (!boost::filesystem::exists(p))
		{
			cerr << "Script file does not exist!" << endl;
			return false;
		}

		string script_extension = p.extension().string();			// Include dot ex) ".sim"

		if (!boost::iequals(script_extension, ".sim"))
		{
			cerr << "Script file extension is not 'sim'!" << endl;
			return false;
		}

		PROJECT_INFO::script_abs_path = p.string();
		PROJECT_INFO::script_abs_dir = p.parent_path().string() + string("\\");
		PROJECT_INFO::script_basename = p.stem().string();
		PROJECT_INFO::script_filename = p.filename().string();

		return true;
	}

	void CreateWorlds()
	{
		DeleteWorlds();

		simulation = new SIMULATION_WORLD;
		opengl = new OPENGL_WORLD;
	}

	void DeleteWorlds()
	{
		DELETE_POINTER(simulation);
		DELETE_POINTER(opengl);
	}

	void ResetSimulation()
	{
		simulation->InitializeFromScript(PROJECT_INFO::script_filename.c_str());
	
		if (opengl)
		{
			opengl->Initialize(simulation);
			opengl->Update();
		}

        total_time.zero();
        simulation_time.zero(); 
	}

	void Update2DDisplayInfo()
	{
		if (opengl)
		{
			ostringstream string_stream;
			int hour, minute;
			string_stream << fixed << setprecision(3);
			hour = (int)total_time.count()/3600;
			minute = ((int)total_time.count()/60)%60;
			string_stream << "acc_dt: " << GetSimulationAccumulatedTime() << "s\n";
			string_stream << "sim_dt: " << simulation_time.count() << "s\n";
			string_stream << "tot_dt: " << hour << "h " << minute << "m ";

			string str = string_stream.str();
			opengl->SetInfoText(str);
		}
	}

	void OneStepSimulation()
	{
		boost::chrono::system_clock::time_point start_time = boost::chrono::system_clock::now();

		if (simulation)
		{
			simulation->AdvanceOneFrame();
		}
		simulation_time = boost::chrono::system_clock::now() - start_time;
		total_time += simulation_time;
		
		if (opengl)
		{
			opengl->Update();
		}
	}

	

	void Render()
	{
		static deque<double> time_deque;
		static boost::chrono::system_clock::time_point prev_time;
		static bool is_first = true;

		Update2DDisplayInfo();

		if (opengl)
		{
			opengl->Render();
		}

		if (is_first)
		{
			prev_time = boost::chrono::system_clock::now();
			is_first = false;
			return;
		}

		boost::chrono::system_clock::time_point current_time = boost::chrono::system_clock::now();
		boost::chrono::duration<double> end_time = current_time - prev_time;

		if (time_deque.size() >= 30)
		{
			time_deque.pop_front();
		}

		time_deque.push_back(end_time.count());

		if (!time_deque.empty())
		{
			double sum = std::accumulate(time_deque.begin(), time_deque.end(), 0.0);
			opengl_fps = (double)time_deque.size() / sum;
		}

		prev_time = current_time;
	}

	void Set2DDisplayInfo(DISPLAY_INFO& info)
	{
		dp_info = info;
		Update2DDisplayInfo();
	}

	void Key(unsigned char c)
	{
		if (opengl)
		{
			opengl->Key(c);
		}
	}

	void KeyWithAlt(unsigned char c)
	{
		if (opengl)
		{
			opengl->KeyWithAlt(c);
		}
	}

	void SpecialKey(int key)
	{
		if (opengl)
		{
			opengl->SpecialKey(key);
		}
	}

	double GetSimulationTimeStep()
	{
		if (simulation)
		{
			return simulation->dt;
		}
		
		return 0;
	}

	double GetSimulationAccumulatedTime()
	{
		if (simulation)
		{
			return simulation->accu_dt;
		}

		return 0;
	}

	int GetCurrentFrame()
	{
		if (simulation)
		{
			return simulation->num_current_frame;
		}

		return 0;
	}

	int GetLastFrame()
	{
		if (simulation)
		{
			return simulation->last_frame;
		}

		return 0;
	}

	void FinalizeSimulation()
	{
	}

	SIMULATION_WORLD* GetSimulationWorld()
	{
		return simulation;
	}

	OPENGL_WORLD* GetOpenGLWorld()
	{
		return opengl;
	}
};


