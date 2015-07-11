#pragma once

#include "OPENGL_OBJECT_BASE.h"

class OPENGL_LEVELSET : public OPENGL_OBJECT_BASE
{
public: // Essential Data
	bool						is_polygonized;
	int							num_threads;
	float						grid_scale;
	MULTITHREADING*				multithreading;
	MARCHING_SQUARES_ALGORITHM  levelset_polygonizer;
	LEVELSET_2D*				levelset_object;

public: // Enumerate
	enum SLEVELSET_DRAW_TYPE
	{
		LEVELSET_DRAW_HIDE			= 0x0000,
		LEVELSET_DRAW_SOLID			= 0x0001,
		LEVELSET_DRAW_WIRE			= 0x0002,
		LEVELSET_DRAW_WIRE_SOLID	= LEVELSET_DRAW_SOLID | LEVELSET_DRAW_WIRE,
	};

public: // Constructor and Destructor
	OPENGL_LEVELSET(const char* display_name, OPENGL_DRIVER* driver, LEVELSET_2D* levelset_object_input, MULTITHREADING* multithreading_input, float grid_scale_input)
		: OPENGL_OBJECT_BASE(display_name, driver), levelset_object(levelset_object_input), num_threads(multithreading_input->num_threads), multithreading(multithreading_input), grid_scale(grid_scale_input), is_polygonized(false)
	{
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_HIDE, "HIDE");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_SOLID, "SOLID");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE, "WIRE");
		RegisterDrawType((int) OPENGL_LEVELSET::LEVELSET_DRAW_WIRE_SOLID, "WIRE_SOLID");
		SetDrawType((int) LEVELSET_DRAW_WIRE);
		is_levelset = true;
	}

	~OPENGL_LEVELSET(void)
	{}

public: // Member Functions
	virtual int NextDrawType()
	{
		switch(GetDrawType())
		{
		case LEVELSET_DRAW_HIDE:
			SetDrawType((int) LEVELSET_DRAW_SOLID);
			break;
		case LEVELSET_DRAW_SOLID:
			SetDrawType((int) LEVELSET_DRAW_WIRE);
			break;
		case LEVELSET_DRAW_WIRE:
			SetDrawType((int) LEVELSET_DRAW_WIRE_SOLID);
			break;
		case LEVELSET_DRAW_WIRE_SOLID:
			SetDrawType((int) LEVELSET_DRAW_HIDE);
			break;
		}

		return GetDrawType();
	}

	virtual int PreviousDrawType()
	{
		switch(GetDrawType())
		{
		case LEVELSET_DRAW_HIDE:
			SetDrawType((int) LEVELSET_DRAW_WIRE_SOLID);
			break;
		case LEVELSET_DRAW_SOLID:
			SetDrawType((int) LEVELSET_DRAW_HIDE);
			break;
		case LEVELSET_DRAW_WIRE:
			SetDrawType((int) LEVELSET_DRAW_SOLID);
			break;
		case LEVELSET_DRAW_WIRE_SOLID:
			SetDrawType((int) LEVELSET_DRAW_WIRE);
			break;
		}

		return GetDrawType();
	}

	virtual void Update()
	{
		if (!levelset_object)
		{
			return;
		}
				
		if (GetDrawType() == LEVELSET_DRAW_HIDE)
		{
			is_polygonized = false;
		}
		else
		{
			levelset_polygonizer.Initialize(multithreading, levelset_object->grid, grid_scale);
			
			levelset_polygonizer.Polygonize(*(levelset_object), true);
			
			is_polygonized = true;
		}
	}

	virtual void Render()
	{
		if (!levelset_object)
		{
			return;
		}

		if (GetDrawType() & LEVELSET_DRAW_SOLID)
		{
			RenderSolid();
		}

		if (GetDrawType() & LEVELSET_DRAW_WIRE)
		{
			RenderWire();
		}
	}

	void RenderSolid()
	{
		if (!is_polygonized)
		{
			Update();
		}

		GetDriver()->SetRenderStatesByMaterial(material);
		for (int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			if (levelset_polygonizer.surfaces[thread_id])
			{
				levelset_polygonizer.surfaces[thread_id]->DrawTriangles(true);
			}
		}
	}

	void RenderWire()
	{
		if (!is_polygonized)
		{
			return;
		}

		GetDriver()->SetDefaultRenderStatesLineMode();
		for (int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			if (levelset_polygonizer.surfaces[thread_id])
			{
				levelset_polygonizer.surfaces[thread_id]->DrawEdges();
				//levelset_polygonizer.surfaces[thread_id]->DrawVertices();
			}
		}
	}

	MARCHING_SQUARES_ALGORITHM* MarchingSquareAlgorithm()
	{
		return& levelset_polygonizer;
	}
};