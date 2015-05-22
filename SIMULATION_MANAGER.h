#pragma once

#include "PROJECT_INFO.h"
#include "SIMULATION_WORLD.h"
#include <string>
#include <boost/chrono.hpp>
#include <iomanip>
#include <deque>
#include <numeric>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

class SIMULATION_MANAGER
{
public: // Essential Data
	SIMULATION_WORLD*			simulation;

public: // Constructor and Destructor
	SIMULATION_MANAGER(void)
		: simulation(0)
	{}

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
	}

	void DeleteWorlds()
	{
		DELETE_POINTER(simulation);
	}

	void ResetSimulation()
	{
		simulation->InitializeFromScript(PROJECT_INFO::script_filename.c_str());
	}

	SIMULATION_WORLD* GetSimulationWorld()
	{
		return simulation;
	}
};


