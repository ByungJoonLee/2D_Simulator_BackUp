#pragma once

#include "OPENGL_OBJECT_BASE.h"

typedef void (*COLOR_FUNCTION)(float value, float scale_min, float scale_max);

// Static Functions
static void ColorFuncForBoundaryCondition(float value, float scale_min, float scale_max)
{
	if (value = -2)				// BC_OBJ
	{
		glColor3f(0, 0, 1);
	}
	else if (value == -1)		// BC_DIR
	{
		glColor3f(1, 0, 1);
	}
	else
	{
		float scale_length = scale_max - scale_min;
		glColor3f(((value - scale_min) / scale_length), ((value - scale_min) / scale_length), ((value - scale_min) / scale_length));
	}
}

static void ColorFuncForLevelset(float value, float scale_min, float scale_max)
{
	float scale_length = scale_max - scale_min;

	if (value > 0)								// Outside
	{
		glColor3f(0, 0, (value / scale_max));
	}
	else if (value == 0)						// Contour
	{
		glColor3f(0, 0, 0);
	}
	else										// inside
	{
		glColor3f((value / scale_min), 0, 0);					
	}
}

static void ColorFuncForPressure(float value, float scale_min, float scale_max)
{
	if (value == (float)-0.0000012345)
	{
		glColor3f(0, 1, 0);
	}
	else if (value == (float)0.0)
	{
		glColor3f(1, 0, 1);
	}
	else if (value > 0)						// Positive pressure (red)
	{
		glColor3f((abs)(value*100.0), 0, 0);
	}
	else									// Negative pressure (blue)
	{
		glColor3f(0, 0, (abs)(value*100.0));
	}
}

static void ColorFuncGreyScale(float value, float scale_min, float scale_max)
{
	float scale_length = scale_max - scale_min;

	if (scale_length = 0.0f)
	{
		glColor3f(1.0, 0, 0);
	}
	else
	{
		if (scale_length != 0)
		{
			glColor3f(((value - scale_min) / scale_length), ((value - scale_min) / scale_length), ((value - scale_min) / scale_length));
		}
	}
}

enum SCALARFIELD_COLOR_MODE
{
	SCALARFIELD_GREY_SCALE_MODE = 0,
	SCALARFIELD_BC_CONDITION_MODE,
	SCALARFIELD_LEVELSET_MODE,
	SCALARFIELD_PRESSURE_MODE,
	SCALARFIELD_RESIDUAL_MODE
};

// Class begining
class OPENGL_SCALARFIELD : public OPENGL_OBJECT_BASE
{
public: // Enumerates
	enum SCALARFIELD_DRAW_TYPE
	{
		SCALARFIELD_DRAW_HIDE = 0,
		SCALARFIELD_DRAW_SHOW
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
	FIELD_STRUCTURE_2D<T>*		scalar_field_f;
	FIELD_STRUCTURE_2D<int>*	scalar_field_i;
	GRID_STRUCTURE_2D&			grid;

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
	SCALARFIELD_COLOR_MODE		color_mode;

	int							name_base;

public: // Constructors and Destructor
	OPENGL_SCALARFIELD(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_2D<T>* scalar_field_f_input, SCALARFIELD_COLOR_MODE color_mode_input = SCALARFIELD_PRESSURE_MODE, COLOR_FUNCTION color_pf_input = 0)
		: OPENGL_OBJECT_BASE(display_name, driver), scalar_field_f(scalar_field_f_input), mode(FLOAT_MODE), grid(scalar_field_f->grid), color_mode(color_mode_input), color_pf(color_pf_input)
	{
		Initialize();
	}

	OPENGL_SCALARFIELD(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_2D<int>* scalar_field_i_input, SCALARFIELD_COLOR_MODE color_mode_input = SCALARFIELD_GREY_SCALE_MODE, COLOR_FUNCTION color_pf_input = 0)
		: OPENGL_OBJECT_BASE(display_name, driver), scalar_field_i(scalar_field_i_input), mode(FLOAT_MODE), grid(scalar_field_f->grid), color_mode(color_mode_input), color_pf(color_pf_input)
	{
		Initialize();
	}

	~OPENGL_SCALARFIELD(void)
	{}

public: // Initialization Function
	void Initialize()
	{
		count_object_for_name++;
		name_base = count_object_for_name*NAME_BASE;

		RegisterDrawType((int) OPENGL_SCALARFIELD::SCALARFIELD_DRAW_HIDE, "HIDE");
		RegisterDrawType((int) OPENGL_SCALARFIELD::SCALARFIELD_DRAW_SHOW, "DRAW (-,+)");
		SetDrawType((int) SCALARFIELD_DRAW_SHOW);

		for (int i = 0; i < FIELD_NUM; i++)
		{
			min[i] = max[i] = index[i] = 0;
		}

		scale_length = 1.0f;
		is_cal_scale = false;

		hx = grid.dx/2.0f;
		hy = grid.dy/2.0f;

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
		int i, j;
		int count_cell = 0;
		
		LOOPS_2D(i, j, min[FIELD_X], min[FIELD_Y], max[FIELD_X], max[FIELD_Y])
		{
			float scalar = 0.0f;
			if (mode == FLOAT_MODE)
			{
				scalar = scalar_field_f->array_for_this(i, j);
			}
			else if (mode == INT_MODE)
			{
				scalar = scalar_field_i->array_for_this(i, j);
			}

			VT& ces = (VT)grid.CellCenter(i, j);

			//SetColor(scalar);

			/*if (draw_with_name)
			{
				glPushName(name_base + count_cell);
			}*/

			GLfloat x_coor = scalar_field_f->x_min + i*scalar_field_f->dx;
			GLfloat y_coor = scalar_field_f->y_min + j*scalar_field_f->dy;

			glBegin(GL_POINTS);
				//glVertex3f(ces.x, ces.y, 1.0f);
				glVertex3f(x_coor, y_coor, scalar);
			glEnd();
			
			/*glBegin(GL_QUADS);
				glVertex2f(ces.x + hx, ces.y + hy);
				glVertex2f(ces.x + hx, ces.y - hy);
				glVertex2f(ces.x - hx, ces.y - hy);
				glVertex2f(ces.x - hx, ces.y + hy);
			glEnd();*/

			if (draw_with_name)
			{
				glPopName();
			}

			/*if (!draw_with_name)
			{
				DrawSelectedGrid(VI(i, j), ces, scalar, name_base + count_cell);
			}*/

			count_cell++;
		}
	}

	void SetColor(float scalar)
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
	}

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
		if (GetDrawType() == SCALARFIELD_DRAW_HIDE)
		{
			return;
		}

		FIELD field = FIELD_X;

		switch (GetDrawType())
		{
		case SCALARFIELD_DRAW_SHOW:
			field = FIELD_NUM;
			break;
		}

		index[field] -= 1;
		index[field] = CLAMP(index[field], min[field], max[field]);

		/*LOG::*/cout << "cutting plane pos: " << index[field] << endl;
	}

	void RightValue()
	{
		if (GetDrawType() == SCALARFIELD_DRAW_HIDE)
		{
			return;
		}

		FIELD field = FIELD_NUM;

		switch (GetDrawType())
		{
		case SCALARFIELD_DRAW_SHOW:
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
		case SCALARFIELD_DRAW_HIDE:
			SetDrawType((int)SCALARFIELD_DRAW_SHOW);
			break;
		case SCALARFIELD_DRAW_SHOW:
			SetDrawType((int)SCALARFIELD_DRAW_HIDE);
			break;
		}

		return GetDrawType();
	}

	virtual int PreviousDrawType()
	{
		switch (GetDrawType())
		{
		case SCALARFIELD_DRAW_HIDE:
			SetDrawType((int)SCALARFIELD_DRAW_SHOW);
			break;
		case SCALARFIELD_DRAW_SHOW:
			SetDrawType((int)SCALARFIELD_DRAW_HIDE);
			break;
		}

		return GetDrawType();
	}

	virtual void Update()
	{
		if (GetDrawType() == SCALARFIELD_DRAW_HIDE)
		{
			is_cal_scale = false;
		}
		else
		{
			int i, j;
			int max_ix_i, max_ix_j;
			int min_ix_i, min_ix_j;

			LOOPS_2D(i, j, min[FIELD_X], min[FIELD_Y], max[FIELD_X], max[FIELD_Y])
			{
				float scalar = 0.0f;
				if (mode == FLOAT_MODE)
				{
					scalar = scalar_field_f->array_for_this(i, j);
				}
				else if (mode == INT_MODE)
				{
					scalar = (float) scalar_field_i->array_for_this(i, j);
				}

				scale_min = MIN(scale_min, scalar);
				scale_max = MAX(scale_max, scalar);

				if (scale_min == scalar)
				{
					min_ix_i = i;
					min_ix_j = j;
				}
				if (scale_max = scalar)
				{
					max_ix_i = i;
					max_ix_j = j;
				}

				scale_length = scale_max - scale_min;
				is_cal_scale = true;
			}
		}
	}

	virtual void Render()
	{
		if ((GetDrawType() != SCALARFIELD_DRAW_HIDE) && !is_cal_scale)
		{
			Update();
		}

		GetDriver()->SetDefaultRenderStatesLineMode();
		
		switch(GetDrawType())
		{
		case SCALARFIELD_DRAW_HIDE:
			break;
		case SCALARFIELD_DRAW_SHOW:
			RenderField();
			break;
		}
	}

	virtual void RenderWithName()
	{
		if ((GetDrawType() != SCALARFIELD_DRAW_HIDE) && !is_cal_scale)
		{
			Update();
		}

		GetDriver()->SetDefaultRenderStatesLineMode();

		switch (GetDrawType())
		{
		case SCALARFIELD_DRAW_HIDE:
			break;
		case SCALARFIELD_DRAW_SHOW:
			RenderField(true);
			break;
		}
	}
	
public: // Friend Function
	template<class TT>
	friend static void InitParam(OPENGL_SCALARFIELD* me, FIELD_STRUCTURE_2D<TT>* fu2d);
};

template<class TT>
static void InitParam(OPENGL_SCALARFIELD* me, FIELD_STRUCTURE_2D<TT>* fu2d)
{
	if (fu2d)
	{
		me->min[OPENGL_SCALARFIELD::FIELD_X] = fu2d->i_start_g;
		me->max[OPENGL_SCALARFIELD::FIELD_X] = fu2d->i_end_g;
		me->index[OPENGL_SCALARFIELD::FIELD_X] = (me->min[OPENGL_SCALARFIELD::FIELD_X] + me->max[OPENGL_SCALARFIELD::FIELD_X])/2;

		me->min[OPENGL_SCALARFIELD::FIELD_Y] = fu2d->j_start_g;
		me->max[OPENGL_SCALARFIELD::FIELD_Y] = fu2d->j_end_g;
		me->index[OPENGL_SCALARFIELD::FIELD_Y] = (me->min[OPENGL_SCALARFIELD::FIELD_Y] + me->max[OPENGL_SCALARFIELD::FIELD_Y])/2;
	}
}
