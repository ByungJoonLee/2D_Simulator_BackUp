#pragma once

#include <string>

class PROJECT_INFO
{
public: // Essential Data
	static std::string app_abs_dir;
	static std::string script_abs_path;
	static std::string script_abs_dir;
	static std::string script_basename;
	static std::string script_filename;

	friend class SIMULATION_MANAGER;

public: // Constructor and Destructor
	PROJECT_INFO(void)
	{}

	PROJECT_INFO(const PROJECT_INFO& other)
	{}

public: // Member Functions
	static std::string GetAppDir()
	{
		return app_abs_dir;
	}

	static std::string GetScriptDir()
	{
		return script_abs_dir;
	}

	static std::string GetScriptAbsPath()
	{
		return script_abs_path;
	}

	static std::string GetScriptAbsDir()
	{
		return script_abs_dir;
	}

	static std::string GetScriptBaseName()
	{
		return script_basename;
	}

	static std::string GetScriptFilename()
	{
		return script_filename;
	}
};


		