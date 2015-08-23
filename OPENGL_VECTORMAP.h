#pragma once

#include "OPENGL_OBJECT_BASE.h"

class OPENGL_VECTORMAP : public OPENGL_OBJECT_BASE
{
public: // Enumerates
	enum FIELD
	{
		FIELD_X = 0,
		FIELD_Y,

		FIELD_NUM,
	};

	enum VECTORFIELD_DRAW_TYPE
	{
		VECTORFIELD_HIDE	= 0,
		VECTORFIELD_DRAW_SHOW  ,
		VECTORFIELD_DRAW_SHOW_X,
		VECTORFIELD_DRAW_SHOW_Y,
	};

public: // Essential Data
	FIELD_STRUCTURE_2D<VT>*			vector_field;
	
	FIELD_STRUCTURE_2D<T>*			vector_field_mac;
	
	FIELD_STRUCTURE_2D<T>*			vector_field_x;
	FIELD_STRUCTURE_2D<T>*			vector_field_y;

	GRID_STRUCTURE_2D&				grid;

	int								min[FIELD_NUM];
	int								max[FIELD_NUM];
	int								index[FIELD_NUM];

	int								name_base;

	float							length_scale;

public: // Constructor and Destructor
	OPENGL_VECTORFIELD(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_2D<VT>* vector_field_input)
		: OPENGL_OBJECT_BASE(display_name, driver), vector_field(vector_field_input), vector_field_mac(0), vector_field_x(0), vector_field_y(0), grid(vector_field_input->grid), length_scale(1.0)
	{
		count_object_for_name++;
		name_base = count_object_for_name*NAME_BASE;

		RegisterDrawType((int) VECTORFIELD_HIDE, "HIDE");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW  , "DRAW (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_X, "DRAW_X (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_Y, "DRAW_Y (-,+,/,*)");
		SetDrawType((int) VECTORFIELD_HIDE);

		for (int i = 0; i < FIELD_NUM; i++)
		{
			min[i] = max[i] = index[i] = 0;
		}

		min[FIELD_X] = vector_field->i_start_g;
		max[FIELD_X] = vector_field->i_end_g;
		index[FIELD_X] = (min[FIELD_X] + max[FIELD_X])/2;

		min[FIELD_Y] = vector_field->j_start_g;
		max[FIELD_Y] = vector_field->j_end_g;
		index[FIELD_Y] = (min[FIELD_Y] + max[FIELD_Y])/2;

		is_velocity = true;
	}

	OPENGL_VECTORFIELD(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_2D<VT>* vector_field_input, FIELD_STRUCTURE_2D<T>* vector_field_mac_input)
		: OPENGL_OBJECT_BASE(display_name, driver), vector_field(vector_field_input), vector_field_mac(vector_field_mac_input), vector_field_x(0), vector_field_y(0), length_scale(1.0), grid(vector_field->grid)
	{
		count_object_for_name++;
		name_base = count_object_for_name*NAME_BASE;
		
		RegisterDrawType((int) VECTORFIELD_HIDE, "HIDE");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW  , "DRAW (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_X, "DRAW_X (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_Y, "DRAW_Y (-,+,/,*)");
		SetDrawType((int) VECTORFIELD_HIDE);

		for (int i = 0; i < FIELD_NUM; i++)
		{
			min[i] = max[i] = index[i] = 0;
		}

		// Speed up variable
		int i_res(vector_field_mac_input->grid.i_res), j_res(vector_field_mac_input->grid.j_res);

		// MAC grid
		// x-component
		if (i_res > j_res)
		{
			min[FIELD_X] = vector_field_mac->i_start_g;
			max[FIELD_X] = vector_field_mac->i_end_g;
			index[FIELD_X] = (min[FIELD_X] + max[FIELD_X])/2;

			min[FIELD_Y] = vector_field_mac->j_start_g;
			max[FIELD_Y] = vector_field_mac->j_end_g;
			index[FIELD_Y] = (min[FIELD_Y] + max[FIELD_Y])/2;

			is_velocity_x = true;
		}		
			
		// y-component
		if (j_res > i_res)
		{
			min[FIELD_X] = vector_field_mac->i_start_g;
			max[FIELD_X] = vector_field_mac->i_end_g;
			index[FIELD_X] = (min[FIELD_X] + max[FIELD_X])/2;

			min[FIELD_Y] = vector_field_mac->j_start_g;
			max[FIELD_Y] = vector_field_mac->j_end_g;
			index[FIELD_Y] = (min[FIELD_Y] + max[FIELD_Y])/2;

			is_velocity_y = true;
		}
	}

	OPENGL_VECTORFIELD(const char* display_name, OPENGL_DRIVER* driver, FIELD_STRUCTURE_2D<VT>* vector_field_input, FIELD_STRUCTURE_2D<T>* vector_field_x_input, FIELD_STRUCTURE_2D<T>* vector_field_y_input)
		: OPENGL_OBJECT_BASE(display_name, driver), vector_field(vector_field_input), vector_field_mac(0), vector_field_x(vector_field_x_input), vector_field_y(vector_field_y_input), length_scale(1.0), grid(vector_field->grid)
	{
		count_object_for_name++;
		name_base = count_object_for_name*NAME_BASE;

		RegisterDrawType((int) VECTORFIELD_HIDE, "HIDE");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW  , "DRAW (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_X, "DRAW_X (-,+,/,*)");
		RegisterDrawType((int) VECTORFIELD_DRAW_SHOW_Y, "DRAW_Y (-,+,/,*)");
		SetDrawType((int) VECTORFIELD_HIDE);

		for (int i = 0; i < FIELD_NUM; i++)
		{
			min[i] = max[i] = index[i] = 0;
		}

		min[FIELD_X] = vector_field->i_start;
		max[FIELD_X] = vector_field->i_end;
		index[FIELD_X] = (min[FIELD_X] + max[FIELD_X])/2;

		min[FIELD_Y] = vector_field->j_start;
		max[FIELD_Y] = vector_field->j_end;
		index[FIELD_Y] = (min[FIELD_Y] + max[FIELD_Y])/2;

		is_velocity = true;
	}

	~OPENGL_VECTORFIELD(void)
	{}

public: // Initialization Function
	void Initialize()
	{}

	void DrawLine(VI& ij, VT& ces, VT& vec, int name_cell, bool draw_with_name)
	{
		/*T mag_vec = vec.Magnitude();
		VT scaled_vec = (T)9/10*vec;
		VT perpen_vec = VT(-vec.y, vec.x);
		perpen_vec.MakeThisUnit();
		VT p_vec = scaled_vec + (T)1/10*perpen_vec;
		VT m_vec = scaled_vec - (T)1/10*perpen_vec;*/

		if (draw_with_name)
		{
			glPushName(name_cell);
		}

		if (!draw_with_name && (driver->GetSelectedNameId() == name_cell))
		{
			glColor3f(0.0f, 1.0f, 0.0f);

			glLineWidth(2.0);
			glBegin(GL_LINES);
				glVertex2f(ces.x, ces.y);
				glVertex2f(ces.x + (vec.x*length_scale), ces.y + (vec.y*length_scale));
			glEnd();
			/*glBegin(GL_LINES);
				glVertex2f(ces.x + (vec.x*length_scale), ces.y + (vec.y*length_scale));
				glVertex2f(ces.x + (p_vec.x*length_scale), ces.y + (p_vec.y*length_scale));
			glEnd();
			glBegin(GL_LINES);
				glVertex2f(ces.x + (vec.x*length_scale), ces.y + (vec.y*length_scale));
				glVertex2f(ces.x + (m_vec.x*length_scale), ces.y + (m_vec.y*length_scale));
			glEnd();*/
			glLineWidth(1.0);
		}
		else
		{
			glColor3f(0.65f, 0.0f, 0.56f);

			glBegin(GL_LINES);
				glVertex2f(ces.x, ces.y);
				glVertex2f(ces.x + (vec.x*length_scale), ces.y + (vec.y*length_scale));
			glEnd();
			/*glBegin(GL_LINES);
				glVertex2f(ces.x + (vec.x*length_scale), ces.y + (vec.y*length_scale));
				glVertex2f(ces.x + (p_vec.x*length_scale), ces.y + (p_vec.x*length_scale));
			glEnd();
			glBegin(GL_LINES);
				glVertex2f(ces.x + (vec.x*length_scale), ces.y + (vec.y*length_scale));
				glVertex2f(ces.x + (m_vec.x*length_scale), ces.y + (m_vec.x*length_scale));
			glEnd();*/
			
		}

		if (draw_with_name)
		{
			glPopName();
		}
	}

	void DrawLine(VI& ij, VT& ces, T& vec_x, T& vec_y, int name_cell, bool draw_with_name)
	{
		if (draw_with_name)
		{
			glPushName(name_cell);
		}

		if (!draw_with_name && (driver->GetSelectedNameId() == name_cell))
		{
			glColor3f(0.0f, 1.0f, 0.0f);

			glLineWidth(2.0);
			glBegin(GL_LINES);
				glVertex2f(ces.x, ces.y);
				glVertex2f(ces.x + (vec_x*length_scale), ces.y + (vec_y*length_scale));
			glEnd();
			glLineWidth(1.0);
		}
		else
		{
			glColor3f(0.65f, 0.0f, 0.56f);

			glBegin(GL_LINES);
				glVertex2f(ces.x, ces.y);
				glVertex2f(ces.x + (vec_x*length_scale), ces.y + (vec_y*length_scale));
			glEnd();
		}

		if (draw_with_name)
		{
			glPopName();
		}
	}

	void RenderXField(bool draw_with_name = false)
	{
		int i, j;
		int count_cell = 0;
		int name_cell = -1;

		LOOPS_2D(i, j, index[FIELD_X], min[FIELD_Y], index[FIELD_X], max[FIELD_Y])
		{
			T& x_component = vector_field_mac->array_for_this(i, j);
			VT& vec_x = VT(x_component, 0);
			VT& ces = (VT)grid.CellCenter(i, j);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j), ces, vec_x, name_cell, draw_with_name);

			count_cell++;
		}
	}

	void RenderYField(bool draw_with_name = false)
	{
		int i, j;
		int count_cell = 0;
		int name_cell = -1;

		LOOPS_2D(i, j, min[FIELD_X], index[FIELD_Y], max[FIELD_X], index[FIELD_Y])
		{
			T& y_component = vector_field_mac->array_for_this(i, j);
			VT& vec_y = VT(0, y_component);
			VT& ces = (VT)grid.CellCenter(i, j);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j), ces, vec_y, name_cell, draw_with_name);

			count_cell++;
		}
	}
	
	void RenderField(bool draw_with_name = false)
	{
		int i, j;
		int count_cell = 0;
		int name_cell = -1;

		for (int i = min[FIELD_X]; i <= max[FIELD_X]; i++)
		{
			for (int j = min[FIELD_Y]; j <= max[FIELD_X]; j++)
			{
				T average_x = (T)0.5*(vector_field_x->array_for_this(i, j) + vector_field_x->array_for_this(i + 1, j));
				T average_y = (T)0.5*(vector_field_y->array_for_this(i, j) + vector_field_y->array_for_this(i, j + 1));
				T mag = sqrt(POW2(average_x) + POW2(average_y));
				
				if (mag >= length_scale)
				{
					length_scale = mag;
				}
			}
		}

		length_scale = (T)sqrt(2)*vector_field->dx/length_scale;
		
		LOOPS_2D(i, j, min[FIELD_X], min[FIELD_Y], max[FIELD_X], max[FIELD_Y])
		{
			T average_x = (T)0.5*(vector_field_x->array_for_this(i, j) + vector_field_x->array_for_this(i + 1, j));
			T average_y = (T)0.5*(vector_field_y->array_for_this(i, j) + vector_field_y->array_for_this(i, j + 1));
			T& x_component = average_x;
			T& y_component = average_y;
			VT& vec = VT(x_component, y_component);
			VT& ces = (VT)grid.CellCenter(i, j);
			name_cell = name_base + count_cell;

			DrawLine(VI(i, j), ces, vec, name_cell, draw_with_name);

			count_cell++;
		}
	}

	void IncrementValue()
	{
		if (GetDrawType() == VECTORFIELD_HIDE)
		{
			return;
		}

		length_scale += 0.1f;

		cout << "Vector field length scale = " << length_scale << endl;
	}

	void DecrementValue()
	{
		if (GetDrawType() == VECTORFIELD_HIDE)
		{
			return;
		}

		length_scale -= 0.1f;
		if (length_scale < 0.1f)
		{
			length_scale = 0.1f;
		}

		cout << "Vector field length scale = " << length_scale << endl;
	}

	void LeftValue()
	{
		if (GetDrawType() == VECTORFIELD_HIDE)
		{
			return;
		}

		FIELD field = FIELD_X;
		switch (GetDrawType())
		{
		case VECTORFIELD_DRAW_SHOW_X:
			field = FIELD_X;
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			field = FIELD_Y;
			break;
		}

		index[field] -= 1;
		index[field] = CLAMP(index[field], min[field], max[field]);
	}

	void RightValue()
	{
		if (GetDrawType() == VECTORFIELD_HIDE)
		{
			return;
		}

		FIELD field = FIELD_X;
		switch (GetDrawType())
		{
		case VECTORFIELD_DRAW_SHOW_X:
			field = FIELD_X;
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			field = FIELD_Y;
			break;
		}

		index[field] += 1;
		index[field] = CLAMP(index[field], min[field], max[field]);
	}

public: // Virtual Functions
	virtual int NextDrawType()
	{
		switch (GetDrawType())
		{
		case VECTORFIELD_HIDE:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_X);
			break;
		case VECTORFIELD_DRAW_SHOW_X:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_Y);
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			SetDrawType((int) VECTORFIELD_HIDE);
			break;
		}

		return GetDrawType();
	}

	virtual int PreviousDrawType()
	{
		switch (GetDrawType())
		{
		case VECTORFIELD_HIDE:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_Y);
			break;
		case VECTORFIELD_DRAW_SHOW_X:
			SetDrawType((int) VECTORFIELD_HIDE);
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			SetDrawType((int) VECTORFIELD_DRAW_SHOW_X);
			break;
		}

		return GetDrawType();
	}

	virtual void Update()
	{}

	virtual void Render()
	{
		GetDriver()->SetDefaultRenderStatesLineMode();
		switch (GetDrawType())
		{
		case VECTORFIELD_HIDE:
			break;
		case VECTORFIELD_DRAW_SHOW:
			RenderField();
			break;
		case VECTORFIELD_DRAW_SHOW_X:
			RenderXField();
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			RenderYField();
			break;
		}
	}

	virtual void RenderWithName()
	{
		GetDriver()->SetDefaultRenderStatesLineMode();
		switch (GetDrawType())
		{
		case VECTORFIELD_HIDE:
			break;
		case VECTORFIELD_DRAW_SHOW_X:
			RenderXField(true);
			break;
		case VECTORFIELD_DRAW_SHOW_Y:
			RenderYField(true);
			break;
		}
	}

	virtual void UserAction(USER_ACTION_TYPE user_action)
	{
		switch (user_action	)
		{
		case OPENGL_OBJECT_BASE::ACTION_INCREMENT:
			IncrementValue();
			break;
		case OPENGL_OBJECT_BASE::ACTION_DECREMENT:
			DecrementValue();
			break;
		case OPENGL_OBJECT_BASE::ACTION_LEFT:
			LeftValue();
			break;
		case OPENGL_OBJECT_BASE::ACTION_RIGHT:
			RightValue();
			break;
		}
	}
};

