#pragma once

#include <GL/glew.h>
#include <GL/glut.h>
#include <string>
#include <iostream>
#include <boost/algorithm/string.hpp>

//////////////////////////////////////////////////////////////////////// 
//						    OPENGL_COLOR							 //
///////////////////////////////////////////////////////////////////////
class OPENGL_COLOR
{
public: // Essential data
	GLfloat r;
	GLfloat g;
	GLfloat b;
	GLfloat a;

public: // Constructors and Destructor
	OPENGL_COLOR(void)
		: r(0.0f), g(0.0f), b(0.0f), a(1.0f)
	{}

	OPENGL_COLOR(GLfloat r_input, GLfloat g_input, GLfloat b_input, GLfloat a_input = 1.0f)
		: r(r_input), g(g_input), b(b_input), a(a_input)
	{}

public: // Member Functions
	GLfloat GetRed(void) const
	{
		return r;
	}

	GLfloat GetGreen(void) const
	{
		return g;
	}

	GLfloat GetBlue(void) const
	{
		return b;
	}

	GLfloat GetAlpha(void) const
	{
		return a;
	}

	void SetRed(GLfloat r_input)
	{
		r = r_input;
	}

	void SetGreen(GLfloat g_input)
	{
		g = g_input;
	}

	void SetBlue(GLfloat b_input)
	{
		b = b_input;
	}

	void SetAlpha(GLfloat a_input)
	{
		a = a_input;
	}

	void Set(GLfloat r_input, GLfloat g_input, GLfloat b_input, GLfloat a_input = 1.0f)
	{
		r = r_input;
		g = g_input;
		b = b_input;
		a = a_input;
	}
};

//////////////////////////////////////////////////////////////////////// 
//						     OPENGL_VEC3							 //
///////////////////////////////////////////////////////////////////////
class OPENGL_VEC3
{
public: // Essential Data
	GLfloat x;
	GLfloat y;
	GLfloat z;

public: // Constructors and Destructor
	OPENGL_VEC3(void)
		: x(0.0f), y(0.0f), z(0.0f)
	{}

	OPENGL_VEC3(GLfloat x_input, GLfloat y_input, GLfloat z_input)
		: x(x_input), y(y_input), z(z_input)
	{}

public: // Operator Overloading
	OPENGL_VEC3& operator+=(const OPENGL_VEC3& p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}

	OPENGL_VEC3& operator-=(const OPENGL_VEC3& p) 
	{
		x -= p.x;
		y -= p.y;
		z -= p.z;
		return *this;
	}

	OPENGL_VEC3& operator*=(GLfloat c)
	{
		x = x*c;
		y = y*c;
		z = z*c;
		return *this;
	}

	OPENGL_VEC3& operator/=(GLfloat c)
	{
		x = x/c;
		y = y/c;
		z = z/c;
		return *this;
	}

public: // Member Functions
	GLfloat GetX(void) const
	{
		return x;
	}

	GLfloat GetY(void) const
	{
		return y;
	}

	GLfloat GetZ(void) const
	{
		return z;
	}

	void SetX(GLfloat x_input)
	{
		x = x_input;
	}

	void SetY(GLfloat y_input)
	{
		y = y_input;
	}

	void SetZ(GLfloat z_input)
	{
		z = z_input;
	}

	void Set(GLfloat x_input, GLfloat y_input, GLfloat z_input)
	{
		x = x_input;
		y = y_input;
		z = z_input;
	}
};

inline const OPENGL_VEC3 operator+(const OPENGL_VEC3& p1, const OPENGL_VEC3& p2)
{
	return OPENGL_VEC3(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

inline const OPENGL_VEC3 operator-(const OPENGL_VEC3& p1, const OPENGL_VEC3& p2)
{
	return OPENGL_VEC3(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

inline const OPENGL_VEC3 operator*(const OPENGL_VEC3& p, GLfloat c)
{
	return OPENGL_VEC3(p.x * c, p.y * c, p.z * c);
}

inline const OPENGL_VEC3 operator*(GLfloat c, const OPENGL_VEC3& p)
{
	return OPENGL_VEC3(p.x * c, p.y * c, p.z * c);
}

inline const OPENGL_VEC3 operator-(const OPENGL_VEC3& p)
{
	return OPENGL_VEC3(-p.x, -p.y, -p.z);
}

inline const OPENGL_VEC3 operator/(const OPENGL_VEC3& p, GLfloat c)
{
	return OPENGL_VEC3(p.x/c, p.y/c, p.z/c);
}

//////////////////////////////////////////////////////////////////////// 
//						      OPENGL_ROT4							 //
///////////////////////////////////////////////////////////////////////
class OPENGL_ROT4
{
public: // Essential Data
	GLfloat s;
	GLfloat x;
	GLfloat y;
	GLfloat z;

public: // Constructors and Destructor
	OPENGL_ROT4(void)
		: s(0.0f), x(0.0f), y(0.0f), z(0.0f)
	{}

	OPENGL_ROT4(GLfloat s_input, GLfloat x_input, GLfloat y_input, GLfloat z_input)
		: s(s_input), x(x_input), y(y_input), z(z_input)
	{}

public: // Member Functions
	GLfloat GetS() const
	{
		return s;
	}

	GLfloat GetX() const
	{
		return x;
	}

	GLfloat GetY() const
	{
		return y;
	}

	GLfloat GetZ() const
	{
		return z;
	}

	void SetS(GLfloat s_input)
	{
		s = s_input;
	}

	void SetX(GLfloat x_input)
	{
		x = x_input;
	}

	void SetY(GLfloat y_input)
	{
		y = y_input;
	}

	void SetZ(GLfloat z_input)
	{
		z = z_input;
	}

	void Set(GLfloat s_input, GLfloat x_input, GLfloat y_input, GLfloat z_input)
	{
		s = s_input;
		x = x_input;
		y = y_input;
		z = z_input;
	}
};

static int OPENGL_ROUND(double d)
{
	return d >= 0.0 ? int(d + 0.5) : int(d - int(d - 1) + 0.5) + int(d - 1);
}

//////////////////////////////////////////////////////////////////////// 
//						    OPENGL_POINT							 //
///////////////////////////////////////////////////////////////////////
class OPENGL_POINT
{
public: // Essential Data
	int xp;
	int yp;

public: // Constructors and Destructor
	OPENGL_POINT(void)
		: xp(0), yp(0)
	{}

	OPENGL_POINT(int x_pos, int y_pos)
		: xp(x_pos), yp(y_pos)
	{}

public: // Operator Overloading
	OPENGL_POINT& operator+=(const OPENGL_POINT& p) 
	{
		xp += p.xp;
		yp += p.yp;
		return *this;
	}

	OPENGL_POINT& operator-=(const OPENGL_POINT& p)
	{
		xp -= p.xp;
		yp -= p.yp;
		return *this;
	}

	OPENGL_POINT& operator*=(double c)
	{
		xp = OPENGL_ROUND(xp*c);
		yp = OPENGL_ROUND(yp*c);
		return *this;
	}

	OPENGL_POINT& operator/=(double c)
	{
		xp = OPENGL_ROUND(xp/c);
		yp = OPENGL_ROUND(yp/c);
		return *this;
	}
	
	friend inline bool operator==(const OPENGL_POINT&, const OPENGL_POINT&);
	friend inline bool operator!=(const OPENGL_POINT&, const OPENGL_POINT&);
	friend inline const OPENGL_POINT operator+(const OPENGL_POINT&, const OPENGL_POINT&);
	friend inline const OPENGL_POINT operator-(const OPENGL_POINT&, const OPENGL_POINT&);
	friend inline const OPENGL_POINT operator*(const OPENGL_POINT&, double);
	friend inline const OPENGL_POINT operator*(double, const OPENGL_POINT&);
	friend inline const OPENGL_POINT operator-(const OPENGL_POINT&);
	friend inline const OPENGL_POINT operator/(const OPENGL_POINT&, double);

public: // Member Functions
	bool IsNull() const
	{
		return xp == 0 && yp == 0;
	}

	int X() const
	{
		return xp;
	}

	int Y() const
	{
		return yp;
	}

	void SetX(int x)
	{
		xp = x;
	}

	void SetY(int y)
	{
		yp = y;
	}

	int ManhattanLength() const
	{
		return abs(X()) + abs(Y());
	}

	int& RX()
	{
		return xp;
	}

	int& RY()
	{
		return yp;
	}
};

inline bool operator==(const OPENGL_POINT& p1, const OPENGL_POINT& p2)
{
	return p1.xp == p2.xp && p1.yp == p2.yp;
}

inline bool operator!=(const OPENGL_POINT& p1, const OPENGL_POINT& p2)
{
	return p1.xp != p2.xp || p1.yp != p2.yp;
}

inline const OPENGL_POINT operator+(const OPENGL_POINT& p1, const OPENGL_POINT& p2)
{
	return OPENGL_POINT(p1.xp + p2.xp, p1.yp + p2.yp);
}

inline const OPENGL_POINT operator-(const OPENGL_POINT& p1, const OPENGL_POINT& p2)
{
	return OPENGL_POINT(p1.xp - p2.xp, p1.yp - p2.yp);
}

inline const OPENGL_POINT operator*(const OPENGL_POINT& p, double c)
{
	return OPENGL_POINT(OPENGL_ROUND(p.xp*c), OPENGL_ROUND(p.yp*c));
}

inline const OPENGL_POINT operator*(double c, const OPENGL_POINT& p)
{
	return OPENGL_POINT(OPENGL_ROUND(p.xp*c), OPENGL_ROUND(p.yp*c));
}

inline const OPENGL_POINT operator-(const OPENGL_POINT& p)
{
	return OPENGL_POINT(-p.xp, -p.yp);
}

inline const OPENGL_POINT operator/(const OPENGL_POINT& p, double c) 
{
	return OPENGL_POINT(OPENGL_ROUND(p.xp/c), OPENGL_ROUND(p.yp/c));
}

//////////////////////////////////////////////////////////////////////// 
//						    OPENGL_SIZE  							 //
///////////////////////////////////////////////////////////////////////
class OPENGL_SIZE
{
public: // Essential Data
	int wd;
	int ht;

public: // Constructors and Destructor
	OPENGL_SIZE(void)
		: wd(-1), ht(-1)
	{}

	OPENGL_SIZE(int w, int h)
		: wd(w), ht(h)
	{}

public: // Member Functions
	bool IsNull() const
	{
		return wd == 0 && ht == 0;
	}

	bool IsEmpty() const
	{
		return wd < 1 || ht < 1;
	}

	bool IsValid() const
	{
		return wd >= 0 && ht >= 0;
	}

	int Width() const
	{
		return wd;
	}

	int Height() const
	{
		return ht;
	}

	void SetWidth(int w)
	{
		wd = w;
	}

	void SetHeight(int h)
	{
		ht = h;
	}

	void Transpose()
	{
		int tmp = wd;
		wd = ht;
		ht = tmp;
	}

	OPENGL_SIZE ExpandedTo(const OPENGL_SIZE& othersize) const
	{
		return OPENGL_SIZE(std::max(wd, othersize.wd), std::max(ht, othersize.ht));
	}

	OPENGL_SIZE BoundedTo(const OPENGL_SIZE& othersize) const
	{
		return OPENGL_SIZE(std::min(wd, othersize.wd), std::min(ht, othersize.ht));
	}

	int& RWidth()
	{
		return wd;
	}

	int& RHeight()
	{
		return ht;
	}

	OPENGL_SIZE& operator+=(const OPENGL_SIZE& s)
	{
		wd += s.wd;
		ht += s.ht;
		return *this;
	}

	OPENGL_SIZE& operator-=(const OPENGL_SIZE& s)
	{
		wd -= s.wd;
		ht -= s.ht;
		return *this;
	}

	OPENGL_SIZE& operator*(double c)
	{
		wd = OPENGL_ROUND(wd*c);
		ht = OPENGL_ROUND(ht*c);
		return *this;
	}

	OPENGL_SIZE& operator/(double c)
	{
		wd = OPENGL_ROUND(wd/c);
		ht = OPENGL_ROUND(ht/c);
		return *this;
	}

	friend inline bool operator==(const OPENGL_SIZE&, const OPENGL_SIZE&);
	friend inline bool operator!=(const OPENGL_SIZE&, const OPENGL_SIZE&);
	friend inline const OPENGL_SIZE operator+(const OPENGL_SIZE&, const OPENGL_SIZE&);
	friend inline const OPENGL_SIZE operator-(const OPENGL_SIZE&, const OPENGL_SIZE&);
	friend inline const OPENGL_SIZE operator*(const OPENGL_SIZE&, double);
	friend inline const OPENGL_SIZE operator*(double, const OPENGL_SIZE&);
	friend inline const OPENGL_SIZE operator/(const OPENGL_SIZE&, double);
};

inline bool operator==(const OPENGL_SIZE& s1, const OPENGL_SIZE& s2)
{
	return s1.wd == s2.wd && s1.ht == s2.ht;
}

inline bool operator!=(const OPENGL_SIZE& s1, const OPENGL_SIZE& s2)
{
	return s1.wd != s2.wd || s1.ht != s2.ht;
}

inline const OPENGL_SIZE operator+(const OPENGL_SIZE& s1, const OPENGL_SIZE& s2)
{
	return OPENGL_SIZE(s1.wd + s2.wd, s1.ht + s2.ht);
}

inline const OPENGL_SIZE operator-(const OPENGL_SIZE& s1, const OPENGL_SIZE& s2)
{
	return OPENGL_SIZE(s1.wd - s2.wd, s1.ht - s2.ht);
}

inline const OPENGL_SIZE operator*(const OPENGL_SIZE& s, double d)
{
	return OPENGL_SIZE(OPENGL_ROUND(s.wd*d), OPENGL_ROUND(s.ht*d));
}

inline const OPENGL_SIZE operator*(double d, const OPENGL_SIZE& s)
{
	return OPENGL_SIZE(OPENGL_ROUND(s.wd*d), OPENGL_ROUND(s.ht*d));
}
inline const OPENGL_SIZE operator/(const OPENGL_SIZE& s, double d)
{
	return OPENGL_SIZE(OPENGL_ROUND(s.wd/d), OPENGL_ROUND(s.ht/d));
}

//////////////////////////////////////////////////////////////////////// 
//						    ENUMERATE       						 //
///////////////////////////////////////////////////////////////////////
enum LIGHT_TYPE
{
	LIGHT_POINT = 0,
	LIGHT_DIRECTIONAL,
	LIGHT_SPOT,
};

enum POLYGON_MODE
{
	POLYGON_SOLID = 0,
	POLYGON_WIREFRAME,
	POLYGON_POINT,
};

enum ANTIALISING_MODE
{
	ANTIALISING_OFF				= 0x0000,

	ANTIALISING_MULTISAMPLE		= 0x0001,
	ANTIALISING_LINE_SMOOTH		= 0x0002,
	ANTIALISING_POINT_SMOOTH	= 0x0004,

	ANTIALISING_ALL				= ANTIALISING_MULTISAMPLE | ANTIALISING_LINE_SMOOTH | ANTIALISING_POINT_SMOOTH,

	ANTIALISING_ALPHA_TO_COVERAGE = 0x0010,
};

enum COLOR_MATERIAL_MODE
{
	COLOR_MATERIAL_NONE =0,

	COLOR_MATERIAL_DIFFUSE,
	COLOR_MATERIAL_AMBIENT,
	COLOR_MATERIAL_EMISSIVE,
	COLOR_MATERIAL_SPECULAR,
	COLOR_MATERIAL_DIFFUSE_AND_AMBIENT
};

enum COMPARISON_FUNC
{
	COMPARISON_NEVER =0,

	COMPARISON_LESSEQUAL,		// <=
	COMPARISON_EQUAL,			// ==
	COMPARISON_LESS,			// <
	COMPARISON_NOTEQUAL,		// !=
	COMPARISON_GREATEREQUAL,	// >=
	COMPARISON_GREATER,			// >
	COMPARISON_ALWAYS,			// always true
};

static std::string MatNames[] =
{
	// 1
	"BRASS",
	"BRONZE",
	"POLISHED_BRONZE",
	"CHROME",
	"COPPER",

	// 6
	"POLISHED_COPPER",
	"GOLD",
	"POLISHED_GOLD",
	"TIN",
	"SILVER",

	// 11
	"POLISHED_SILVER",
	"EMERALD",
	"JADE",
	"OBSIDIAN",
	"PEARL",

	// 16
	"RUBY",
	"TURQUOISE",
	"BLACK_PLASTIC",
	"CYAN_PLASTIC",
	"GREEN_PLASTIC",

	// 21
	"RED_PLASTIC",
	"WHITE_PLASTIC",
	"YELLOW_PLASTIC",
	"BLACK_RUBBER",
	"CYAN_RUBBER",

	// 26
	"GREEN_RUBBER",
	"RED_RUBBER",
	"WHITE_RUBBER",
	"YELLOW_RUBBER",
	"BRIGHT_WHITE",

	// 31
	"LESS_BRIGHT_WHITE",
	"WARM_WHITE",
	"COOL_WHITE",
	"CYAN_PLASTIC_TRANSPARENT",
};

static GLfloat Material[][11] = 
{
	// 1
	//BRASS
	0.329412f, 0.223529f, 0.027451f,
	0.780392f, 0.568627f, 0.113725f, 
	0.992157f, 0.941176f, 0.807843f,
	27.9f, 1.0f,

	//BRONZE
	0.2125f, 0.1275f, 0.054f,
	0.714f, 0.4284f, 0.18144f, 
	0.393548f, 0.271906f, 0.166721f, 
	25.6f, 1.0f,

	// POLISHED_BRONZE
	0.25f, 0.148f, 0.06475f,
	0.4f, 0.2368f, 0.1036f,
	0.774597f, 0.458561f, 0.200621f,
	76.8f, 1.0f,

	//CHROME
	0.25f, 0.25f, 0.25f,
	0.4f, 0.4f, 0.4f, 
	0.774597f, 0.774597f, 0.774597f, 
	76.8f, 1.0f,

	//COPPER
	0.19125f, 0.0735f, 0.0225f,
	0.7038f, 0.27048f, 0.0828f, 
	0.256777f, 0.137622f, 0.086014f, 
	12.8f, 1.0f,

	// 6
	//POLISHED_COPPER
	0.2295f, 0.08825f, 0.0275f,
	0.5508f, 0.2118f, 0.066f,
	0.580594f, 0.223257f, 0.0695701f,
	51.2f, 1.0f,

	//GOLD
	0.24725f, 0.1995f, 0.0745f,
	0.75164f, 0.60648f, 0.22648f, 
	0.628281f, 0.555802f, 0.366065f, 
	51.2f, 1.0f,

	//POLISHED_GOLD
	0.24725f, 0.2245f, 0.0645f,
	0.34615f, 0.3143f, 0.0903f,
	0.797357f, 0.723991f, 0.208006f,
	83.2f, 1.0f,

	//TIN
	0.105882f, 0.058824f, 0.113725f,
	0.427451f, 0.470588f, 0.541176f,
	0.333333f, 0.333333f, 0.521569f,
	9.84615f, 1.0f,

	//SILVER
	0.19225f, 0.19225f, 0.19225f,
	0.50754f, 0.50754f, 0.50754f, 
	0.508273f, 0.508273f, 0.508273f, 
	51.2f, 1.0f,

	// 11
	//POLISHED_SILVER
	0.23125f, 0.23125f, 0.23125f,
	0.2775f, 0.2775f, 0.2775f,
	0.773911f, 0.773911f, 0.773911f,
	89.6f, 1.0f,

	//EMERALD
	0.0215f, 0.1745f, 0.0215f,
	0.07568f, 0.61424f, 0.07568f, 
	0.633f, 0.727811f, 0.633f, 
	76.8f, 1.0f,

	//JADE
	0.135f, 0.2225f, 0.1575f,
	0.54f, 0.89f, 0.63f, 
	0.316228f, 0.316228f, 0.316228f, 
	12.8f, 1.0f,

	//OBSIDIAN
	0.05375f, 0.05f, 0.06625f,
	0.18275f, 0.17f, 0.22525f, 
	0.332741f, 0.328634f, 0.346435f, 
	38.4f, 1.0f,

	//PEARL
	0.25f, 0.20725f, 0.20725f,
	1.0f, 0.829f, 0.829f, 
	0.296648f, 0.296648f, 0.296648f, 
	11.264f, 1.0f,

	// 16
	//RUBY
	0.1745f, 0.01175f, 0.01175f,
	0.61424f, 0.04136f, 0.04136f, 
	0.727811f, 0.626959f, 0.626959f, 
	76.8f, 1.0f,

	//TURQUOISE
	0.1f, 0.18725f, 0.1745f,
	0.396f, 0.74151f, 0.69102f, 
	0.297254f, 0.30829f, 0.306678f, 
	12.8f, 1.0f,

	//BLACK PLASTIC
	0.0f, 0.0f, 0.0f, 
	0.01f, 0.01f, 0.01f,
	0.50f, 0.50f, 0.50f, 
	32.0f, 1.0f,

	//CYAN PLASTIC
	0.0f, 0.1f, 0.06f, 
	0.0f, 0.50980392f, 0.50980392f,
	0.50196078f, 0.50196078f, 0.50196078f, 
	32.0f, 1.0f,

	//GREEN PLASTIC
	0.0f, 0.0f, 0.0f,
	0.1f, 0.35f, 0.1f, 
	0.45f, 0.55f, 0.45f, 
	32.0f, 1.0f,

	// 21
	//RED PLASTIC
	0.0f, 0.0f, 0.0f, 
	0.5f, 0.0f, 0.0f,
	0.7f, 0.6f, 0.6f, 
	32.0f, 1.0f,

	//WHITE PLASTIC
	0.0f, 0.0f, 0.0f, 
	0.55f, 0.55f, 0.55f,
	0.70f, 0.70f, 0.70f, 
	32.0f, 1.0f,

	//YELLOW PLASTIC
	0.0f, 0.0f, 0.0f, 
	0.5f, 0.5f, 0.0f,
	0.60f, 0.60f, 0.50f, 
	32.0f, 1.0f,

	//BLACK RUBBER
	0.02f, 0.02f, 0.02f, 
	0.01f, 0.01f, 0.01f,
	0.4f, 0.4f, 0.4f, 
	10.0f, 1.0f,

	//CYAN RUBBER
	0.0f, 0.05f, 0.05f, 
	0.4f, 0.5f, 0.5f,
	0.04f, 0.7f, 0.7f, 
	10.0f, 1.0f,

	// 26
	//GREEN RUBBER
	0.0f, 0.05f, 0.0f, 
	0.4f, 0.5f, 0.4f,
	0.04f, 0.7f, 0.04f, 
	10.0f, 1.0f,

	//RED RUBBER
	0.05f, 0.0f, 0.0f, 
	0.5f, 0.4f, 0.4f,
	0.7f, 0.04f, 0.04f, 
	10.0f, 1.0f,

	//WHITE RUBBER
	0.05f, 0.05f, 0.05f, 
	0.5f, 0.5f, 0.5f,
	0.7f, 0.7f, 0.7f, 
	10.0f, 1.0f,

	//YELLOW RUBBER
	0.05f, 0.05f, 0.0f, 
	0.5f, 0.5f, 0.4f,
	0.7f, 0.7f, 0.04f, 
	10.0f, 1.0f,

	//BRIGHT WHITE
	0.2f, 0.2f, 0.2f,
	1.0f, 1.0f, 1.0f, 
	0.8f, 0.8f, 0.8f, 
	51.2f, 1.0f,

	// 31
	//LESS BRIGHT WHITE
	0.2f, 0.2f, 0.2f,
	0.8f, 0.8f, 0.8f, 
	0.5f, 0.5f, 0.5f, 
	44.8f, 1.0f,

	//WARMISH WHITE
	0.3f, 0.2f, 0.2f,
	1.0f, 0.9f, 0.8f, 
	0.4f, 0.2f, 0.2f, 
	44.8f, 1.0f,

	//COOLISH WHITE
	0.2f, 0.2f, 0.3f,
	0.8f, 0.9f, 1.0f, 
	0.2f, 0.2f, 0.4f, 
	44.8f, 1.0f,

	//CYAN PLASTIC TRANSPARENT
	0.0f, 0.1f, 0.06f, 
	0.0f, 0.50980392f, 0.50980392f,
	0.50196078f, 0.50196078f, 0.50196078f, 
	32.0f, 0.2f,
};

struct OPENGL_MATERIAL
{
	OPENGL_COLOR			ambient_color;
	OPENGL_COLOR			diffuse_color;
	OPENGL_COLOR			specular_color;
	OPENGL_COLOR			emissive_color;
	GLfloat					shininess;

	GLfloat					thickness; // point and line thickness
	GLboolean				is_gouraud_shading;
	GLboolean				is_lighting;
	GLboolean				is_backface_culling;
	GLboolean				is_frontface_culling;
	GLboolean				is_normalize_normals;
	GLboolean				is_z_writable;

	POLYGON_MODE			polygon_mode;
	COLOR_MATERIAL_MODE		color_material;
	COMPARISON_FUNC			z_buffer;

	OPENGL_MATERIAL()
		: ambient_color(0.0f, 0.0f, 0.0f)
		, diffuse_color(0.0f, 0.0f, 0.0f)
		, specular_color(0.0f, 0.0f, 0.0f)
		, emissive_color(0.0f, 0.0f, 0.0f)
		, shininess(0.0f)
		, thickness(1.0f)
		, is_gouraud_shading(true)
		, is_lighting(true)
		, is_backface_culling(false)
		, is_frontface_culling(false)
		, is_normalize_normals(false)
		, is_z_writable(true)
		, polygon_mode(POLYGON_SOLID)
		, color_material(COLOR_MATERIAL_NONE)
		, z_buffer(COMPARISON_LESSEQUAL)
	{
		SetMaterialByName("CHROME");
	}

	void SetMaterialByName(const char* color_name)
	{
		int num_material = sizeof(MatNames)/sizeof(MatNames[0]);
		for(int i = 0; i < num_material; ++i)
		{
			if(boost::iequals(MatNames[i], color_name))
			{
				GLfloat *pMat = Material[i];
				ambient_color.Set(pMat[0], pMat[1], pMat[2], pMat[10]);
				diffuse_color.Set(pMat[3], pMat[4], pMat[5], pMat[10]);
				specular_color.Set(pMat[6], pMat[7], pMat[8], pMat[10]);
				shininess = pMat[9];
				break;
			}
		}
	}
};

struct OPENGL_LIGHT_PROPERTY
{
	LIGHT_TYPE		light_type;

	OPENGL_COLOR	ambient_color;
	OPENGL_COLOR	diffuse_color;
	OPENGL_COLOR	specular_color;
	
	GLfloat			constant_attenuation;
	GLfloat			linear_attenuation;
	GLfloat			quadratic_attenuation;

	GLfloat			spot_exponent;
	GLfloat			spot_cutoff;

	OPENGL_VEC3		position;
	OPENGL_VEC3		direction;

	OPENGL_LIGHT_PROPERTY()
		: light_type(LIGHT_POINT)
		, ambient_color(0.2f, 0.2f, 0.2f)
		, diffuse_color(1.0f, 1.0f, 1.0f)
		, specular_color(1.0f, 1.0f, 1.0f)
		, constant_attenuation(1.0f)
		, linear_attenuation(0.0f)
		, quadratic_attenuation(0.0f)
		, spot_exponent(2.0f)
		, spot_cutoff(45.0f)
		, position(0.0f, 0.0f, 0.0f)
		, direction(0.0f, 0.0f, 0.0f)
	{}
};








