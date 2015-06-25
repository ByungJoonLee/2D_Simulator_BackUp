#pragma once

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <io.h>
#include <string.h>
#include <fstream>
//#include <windef.h>

#include <Windows.h>

class FILE_UTILITIES
{
public: // Constructor and Destructor
	FILE_UTILITIES(void)
	{}

	~FILE_UTILITIES(void)
	{}

public:	// Member Functions
	int Execute_Process(const std::string& cmd)
	{
		return system(cmd.c_str());
	}

	bool Directory_Exists(const std::string& dirname)
	{
		DWORD attr=GetFileAttributes(dirname.c_str());
		
		return((attr!=-1) && (attr&FILE_ATTRIBUTE_DIRECTORY));
	}

	bool Create_Directory(const std::string& dirname)
	{
		if (!Directory_Exists(dirname))
		{
			CreateDirectory(dirname.c_str(), 0);
			if (GetLastError() == ERROR_ALREADY_EXISTS)
			{
				std::cout << "ERROR_ALREADY_EXISTS : " << dirname.c_str() << std::endl;
			}
			if (GetLastError() == ERROR_PATH_NOT_FOUND)
			{
				std::cout << "ERROR_PATH_NOT_FOUND : " << dirname.c_str() << std::endl;
			}
			if (!Directory_Exists(dirname))
			{
				std::cerr << dirname.c_str() << " is NOT created Successfully." << std::endl;
			}
			std::cout << dirname.c_str() << " is Create Successfully." << std::endl;
		}
		return true;
	}
};


