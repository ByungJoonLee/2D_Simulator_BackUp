#pragma once

#include "OPENGL_OBJECT_BASE.h"
#include "SCRIPT_READER.h"
#include <vector>

class OPENGL_SOLVER_BASE
{
public: // Essential Data
	vector<OPENGL_OBJECT_BASE*> all_objects;

public: // Constructor and Destructor
	OPENGL_SOLVER_BASE(void)
	{}

	virtual ~OPENGL_SOLVER_BASE(void)
	{
		for (vector<OPENGL_OBJECT_BASE*>::iterator it = all_objects.begin(); it != all_objects.end(); ++it)
		{
			delete (*it);
		}
	}

public: // Member Functions
	std::vector<OPENGL_OBJECT_BASE*>& GetObjects()
	{
		return all_objects;
	}

	void AddObject(OPENGL_OBJECT_BASE* obj)
	{
		all_objects.push_back(obj);
	}

public: // Virtual Functions
	virtual void Key(unsigned char c)
	{}

	virtual void KeyWithAlt(unsigned char c)
	{}

	virtual void LoadOpenGLSettings(SCRIPT_BLOCK* script_block) = 0;
};