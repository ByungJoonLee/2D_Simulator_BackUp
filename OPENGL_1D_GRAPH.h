#pragma once

#include "OPENGL_OBJECT_BASE.h"
#include "FIELD_STRUCTURE_1D.h"

typedef void (*COLOR_FUNCTION)(float value, float scale_min, float scale_max);

enum FIELD_1D_COLOR_MODE
{
	FIELD_1D_GREY_SCALE_MODE = 0,
	FIELD_1D_BC_CONDITION_MODE,
	FIELD_1D_LEVELSET_MODE,
	FIELD_1D_PRESSURE_MODE,
	FIELD_1D_RESIDUAL_MODE
};

// Class begining
class OPENGL_1D_GRAPH : public OPENGL_OBJECT_BASE
{
public: // Enumerates
	enum FIELD_1D_DRAW_TYPE
	{
		FIELD_1D_DRAW_HIDE = 0,
		FIELD_1D_DRAW_SHOW
	};

	enum MODE
	{
		FLOAT_MODE = 0,
		INT_MODE
	};

	enum FIELD
	{
		FIELD_X = 0,
		FIELD_Y,
		
		FIELD_NUM
	};

public: // Essential Data
	FIELD_STRUCTURE_1D<T>*		scalar_field_f;
	FIELD_STRUCTURE_1D<int>*	scalar_field_i;
	GRID_STRUCTURE_1D&			grid;

	int							min[FIELD_NUM];
	int							max[FIELD_NUM];
	int							index[FIELD_NUM];

	float						scale_min;
	float						scale_max;
	float						scale_length;

	float						hx;
	float						hy;
	
	MODE						mode;

	bool						is_cal_scale;

	COLOR_FUNCTION				color_pf;
	FIELD_1D_COLOR_MODE		    color_mode;

	int							name_base;

public: // Constructors and Destructor
	OPENGL_1D_GRAPH(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_1D<T>* scalar_field_f_input, FIELD_1D_COLOR_MODE color_mode_input = FIELD_1D_PRESSURE_MODE, COLOR_FUNCTION color_pf_input = 0)
		: OPENGL_OBJECT_BASE(display_name, driver), scalar_field_f(scalar_field_f_input), mode(FLOAT_MODE), grid(scalar_field_f->grid), color_mode(color_mode_input), color_pf(color_pf_input)
	{
		Initialize();
	}

	OPENGL_1D_GRAPH(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_1D<int>* scalar_field_i_input, FIELD_1D_COLOR_MODE color_mode_input = FIELD_1D_GREY_SCALE_MODE, COLOR_FUNCTION color_pf_input = 0)
		: OPENGL_OBJECT_BASE(display_name, driver), scalar_field_i(scalar_field_i_input), mode(FLOAT_MODE), grid(scalar_field_f->grid), color_mode(color_mode_input), color_pf(color_pf_input)
	{
		Initialize();
	}

	~OPENGL_1D_GRAPH(void)
	{}

public: // Initialization Function
	void Initialize()
	{
		count_object_for_name++;
		name_base = count_object_for_name*NAME_BASE;

		RegisterDrawType((int) OPENGL_1D_GRAPH::FIELD_1D_DRAW_HIDE, "HIDE");
		RegisterDrawType((int) OPENGL_1D_GRAPH::FIELD_1D_DRAW_SHOW, "DRAW (-,+)");
		SetDrawType((int) FIELD_1D_DRAW_SHOW);

		for (int i = 0; i < FIELD_NUM; i++)
		{
			min[i] = max[i] = index[i] = 0;
		}

		scale_length = 1.0f;
		is_cal_scale = false;

		hx = grid.dx/2.0f;
		
		switch(mode)
		{
		case FLOAT_MODE:
			InitParam(this, scalar_field_f);
			break;
		case INT_MODE:
			InitParam(this, scalar_field_i);
			break;
		}
		
		is_scalar = true;
	}
	
	void RenderField(bool draw_with_name = false)
	{
		int i(0);
		int count_cell = 0;
		
		// For scaling of graph
		float max_s = 0.0f;
		for (int i = min[FIELD_X]; i <= max[FIELD_X]; i++)
		{
			if (max_s <= scalar_field_f->array_for_this(i))
			{
				max_s = scalar_field_f->array_for_this(i);
			}
		}

		for (int i = min[FIELD_X]; i <= max[FIELD_X]; i++)
		{
			float scalar = 0.0f;
			if (mode == FLOAT_MODE)
			{
				scalar = scalar_field_f->array_for_this(i);
			}
			else if (mode == INT_MODE)
			{
				scalar = scalar_field_i->array_for_this(i);
			}

			//SetColor(scalar);

			/*if (draw_with_name)
			{
				glPushName(name_base + count_cell);
			}*/

			GLfloat x_coor = scalar_field_f->x_min + i*scalar_field_f->dx;
			// Translation for drawing
			GLfloat x_cen = (scalar_field_f->x_min + (scalar_field_f->x_min + scalar_field_f->i_end*scalar_field_f->dx))/(T)2;

			glBegin(GL_POINTS);
				//glVertex3f(ces.x, ces.y, 1.0f);
				glVertex3f(x_coor - x_cen, scalar/max_s - 0.5f, 0.0f);
			glEnd();
			
			if (draw_with_name)
			{
				glPopName();
			}

			count_cell++;
		}
	}

	/*void SetColor(float scalar)
	{
		if (color_pf)
		{
			color_pf(scalar, scale_min, scale_max);
		}
		else
		{
			switch (color_mode)
			{
			case SCALARFIELD_GREY_SCALE_MODE:
				ColorFuncGreyScale(scalar, scale_min, scale_max);
				break;
			case SCALARFIELD_BC_CONDITION_MODE:
				ColorFuncForBoundaryCondition(scalar, scale_min, scale_max);
				break;
			case SCALARFIELD_LEVELSET_MODE:
				ColorFuncForLevelset(scalar, scale_min, scale_max);
				break;
			case SCALARFIELD_PRESSURE_MODE:
				ColorFuncForPressure(scalar, scale_min, scale_max);
				break;
			}
		}
	}*/

	void DrawSelectedGrid(VI& ijk, VT& ces, float scalar, int grid_name)
	{
		if (driver->GetSelectedNameId() == grid_name)
		{
			glColor3f(0.0f, 1.0f, 0.0f);

			glPushMatrix();
				glTranslatef(ces.x, ces.y, 0.0f);
				GetDriver()->DrawWireBox(hx, hy);
			glPopMatrix();
		}
	}

	void LeftValue()
	{
		if (GetDrawType() == FIELD_1D_DRAW_HIDE)
		{
			return;
		}

		FIELD field = FIELD_X;

		switch (GetDrawType())
		{
		case FIELD_1D_DRAW_SHOW:
			field = FIELD_NUM;
			break;
		}

		index[field] -= 1;
		index[field] = CLAMP(index[field], min[field], max[field]);

		/*LOG::*/cout << "cutting plane pos: " << index[field] << endl;
	}

	void RightValue()
	{
		if (GetDrawType() == FIELD_1D_DRAW_HIDE)
		{
			return;
		}

		FIELD field = FIELD_NUM;

		switch (GetDrawType())
		{
		case FIELD_1D_DRAW_SHOW:
			field = FIELD_NUM;
		break;
		}

		index[field] += 1;
		index[field] = CLAMP(index[field], min[field], max[field]);

		/*LOG::*/cout << "cutting plane pos: " << index[field] << endl;
	}

	void UserAction(USER_ACTION_TYPE user_action)
	{
		switch (user_action)
		{
		case ACTION_INCREMENT:
			break;
		case ACTION_DECREMENT:
			break;
		case ACTION_LEFT:
			LeftValue();
			break;
		case ACTION_RIGHT:
			RightValue();
			break;
		}
	}

public: // Virtual Functions
	virtual int NextDrawType()
	{
		switch (GetDrawType())
		{
		case FIELD_1D_DRAW_HIDE:
			SetDrawType((int)FIELD_1D_DRAW_SHOW);
			break;
		case FIELD_1D_DRAW_SHOW:
			SetDrawType((int)FIELD_1D_DRAW_HIDE);
			break;
		}

		return GetDrawType();
	}

	virtual int PreviousDrawType()
	{
		switch (GetDrawType())
		{
		case FIELD_1D_DRAW_HIDE:
			SetDrawType((int)FIELD_1D_DRAW_SHOW);
			break;
		case FIELD_1D_DRAW_SHOW:
			SetDrawType((int)FIELD_1D_DRAW_HIDE);
			break;
		}

		return GetDrawType();
	}

	virtual void Update()
	{
		if (GetDrawType() == FIELD_1D_DRAW_HIDE)
		{
			is_cal_scale = false;
		}
		else
		{
			int i(0);
			int max_ix_i;
			int min_ix_i;

			for (int i = min[FIELD_X]; i <= max[FIELD_X]; i++)
			{
				float scalar = 0.0f;
				if (mode == FLOAT_MODE)
				{
					scalar = scalar_field_f->array_for_this(i);
				}
				else if (mode == INT_MODE)
				{
					scalar = (float) scalar_field_i->array_for_this(i);
				}

				scale_min = MIN(scale_min, scalar);
				scale_max = MAX(scale_max, scalar);

				if (scale_min == scalar)
				{
					min_ix_i = i;
				}
				if (scale_max = scalar)
				{
					max_ix_i = i;
				}

				scale_length = scale_max - scale_min;
				is_cal_scale = true;
			}
		}
	}

	virtual void Render()
	{
		if ((GetDrawType() != FIELD_1D_DRAW_HIDE) && !is_cal_scale)
		{
			Update();
		}

		GetDriver()->SetDefaultRenderStatesLineMode();
		
		switch(GetDrawType())
		{
		case FIELD_1D_DRAW_HIDE:
			break;
		case FIELD_1D_DRAW_SHOW:
			RenderField();
			break;
		}
	}

	virtual void RenderWithName()
	{
		if ((GetDrawType() != FIELD_1D_DRAW_HIDE) && !is_cal_scale)
		{
			Update();
		}

		GetDriver()->SetDefaultRenderStatesLineMode();

		switch (GetDrawType())
		{
		case FIELD_1D_DRAW_HIDE:
			break;
		case FIELD_1D_DRAW_SHOW:
			RenderField(true);
			break;
		}
	}
	
public: // Friend Function
	template<class TT>
	friend static void InitParam(OPENGL_1D_GRAPH* me, FIELD_STRUCTURE_1D<TT>* fu1d);
};

template<class TT>
static void InitParam(OPENGL_1D_GRAPH* me, FIELD_STRUCTURE_1D<TT>* fu1d)
{
	if (fu1d)
	{
		me->min[OPENGL_1D_GRAPH::FIELD_X] = fu1d->i_start;
		me->max[OPENGL_1D_GRAPH::FIELD_X] = fu1d->i_end;
		me->index[OPENGL_1D_GRAPH::FIELD_X] = (me->min[OPENGL_1D_GRAPH::FIELD_X] + me->max[OPENGL_1D_GRAPH::FIELD_X])/2;
	}
}
