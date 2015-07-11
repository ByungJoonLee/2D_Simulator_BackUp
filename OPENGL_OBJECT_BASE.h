#pragma once

#include "OPENGL_DRIVER.h"
#include <string>
#include <map>

#define NAME_BASE 10000000;

class OPENGL_OBJECT_BASE
{
public: // Essential Data
	int					id;
	string				name;
	int					draw_type;
	bool				is_visible;

	OPENGL_ROT4			rotation;
	OPENGL_VEC2			position;
	OPENGL_VEC2			scale;
	OPENGL_VEC2			length;
	OPENGL_VEC2			center;

	OPENGL_MATERIAL		material;

	//BOX					bounding_box;

	OPENGL_DRIVER*		driver;

	map<int, string>	map_draw_types;

	static int			count_object_for_name;

public: // Constructor and Destructor
	OPENGL_OBJECT_BASE(const char* display_name, OPENGL_DRIVER* driver_input)
		: id(GenerateID()), name(display_name), draw_type(0), position(0.0, 0.0), rotation(0.0, 0.0, 0.0, 0.0), scale(1.0, 1.0), driver(driver_input)
	{}

	~OPENGL_OBJECT_BASE(void)
	{}

public: // Member Functions
	enum USER_ACTION_TYPE
	{
		ACTION_INCREMENT = 0,
		ACTION_DECREMENT,

		ACTION_LEFT,
		ACTION_RIGHT,
	};

	static int GenerateID()
	{
		static int genID = 0;
		return genID++;
	}

	int GetID()
	{
		return id;
	}

	OPENGL_DRIVER* GetDriver()
	{
		return driver;
	}

	// For rendering
	void Draw()
	{
		PreDraw();
		Render();
		PostDraw();
	}

	void DrawWithName()
	{
		PreDraw();
		RenderWithName();
		PostDraw();
	}

	void PreDraw()
	{
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
			glScalef(scale.GetX(), scale.GetY(), 0.0f);
			glTranslatef(position.GetX(), position.GetY(), 0.0f);
	}

	void PostDraw()
	{
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
	}

	void RenderWithName()
	{
		glPushName(GetID());
		Render();
		glPopName();
	}
	
	// Set up the basic properties
	void SetName(const char* display_name)
	{
		name = display_name;
	}

	string& GetName()
	{
		return name;
	}

	void SetPosition(const OPENGL_VEC2& position_input)
	{
		position = position_input;
	}

	void SetPosition(GLfloat x, GLfloat y)
	{
		position = OPENGL_VEC2(x, y);
	}
	
	const OPENGL_VEC2& GetPosition() const
	{
		return position;
	}

	void SetRotation(const OPENGL_ROT4& rotation_input)
	{
		rotation = rotation_input;
	}
	
	void SetRotation(GLfloat s, GLfloat x, GLfloat y, GLfloat z)
	{
		rotation = OPENGL_ROT4(s, x, y, z);
	}

	const OPENGL_ROT4& GetRotation() const
	{
		return rotation;
	}

	void SetScale(const OPENGL_VEC2& scale_input)
	{
		scale = scale_input;
	}

	void SetScale(GLfloat x, GLfloat y) 
	{
		scale = OPENGL_VEC2(x, y);
	}

	const OPENGL_VEC2& GetScale() const
	{
		return scale;
	}
	
	void SetLength(const OPENGL_VEC2& length_input)
	{
		length = length_input;
	}

	void SetLength(GLfloat x, GLfloat y)
	{
		length = OPENGL_VEC2(x, y);
	}

	const OPENGL_VEC2& GetLength() const
	{
		return length;
	}

	void SetCenter(const OPENGL_VEC2& center_input)
	{
		center = center_input;
	}

	void SetCenter(GLfloat x, GLfloat y)
	{
		center = OPENGL_VEC2(x, y);
	}

	const OPENGL_VEC2& GetCenter() const
	{
		return center;
	}

	void SetMaterial(const OPENGL_MATERIAL& material_input)
	{
		material = material_input;
	}
	
	const OPENGL_MATERIAL& GetMaterial() const
	{
		return material;
	}

	// Draw type
	int GetDrawType()
	{
		return draw_type;
	}

	void SetDrawType(int draw_type_input)
	{
		draw_type = draw_type_input;
	}

	map<int, string>& GetDrawTypeMap()
	{
		return map_draw_types;
	}

	// Register and Cout Functions
	void RegisterDrawType(int draw_type_input, string type_name)
	{
		map_draw_types[draw_type] = type_name;
	}

	void CoutDrawType()
	{
		cout << GetDrawTypeMap()[GetDrawType()] << endl;
	}
	
	// Virtual functions
	virtual void Render() = 0;
	
	virtual void Update()
	{}

	virtual void RenderExtra()
	{}

	virtual int PreviousDrawType() = 0;
	

	virtual int NextDrawType() = 0;
	
	// Values
	virtual void UserAction(USER_ACTION_TYPE user_action)
	{}

};