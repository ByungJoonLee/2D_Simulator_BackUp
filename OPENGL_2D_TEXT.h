/********************************************************************
	OPENGL_2D_TEXT.h
*********************************************************************/

#pragma once

#include "OPENGL_COMMON.h"
#include <FTGL/ftgl.h>

enum OPENGL_FONT_TYPE
{
	BITMAP_FONT = 0,
	PIXMAP_FONT,
};

enum FONT_NAME
{
	FONT_ARIAL = 0,
	FONT_CONSOLAS,
	FONT_COURIER_NEW,
	FONT_GEORGIA,
	FONT_IMPACT,
	FONT_PLATINO,
	FONT_TAHOMA,
	FONT_TIMES_NEW_ROMAN,
	FONT_VERDANA,
};

enum FONT_STYLE
{
	STYLE_REGULAR		= 0x0000,
	STYLE_BOLD			= 0x0001,
	STYLE_ITALIC		= 0x0002,
	STYLE_BOLD_ITALIC	= STYLE_BOLD | STYLE_ITALIC,
};

enum TEXT_BLOCK_LAYOUT
{
	TEXT_BLOCK_TOP = 0,
	TEXT_BLOCK_LEFT,
	TEXT_BLOCK_RIGHT,
	TEXT_BLOCK_BOTTOM,
};

//////////////////////////////////////////////////////////////////////////
// OPENGL_2D_BASE
class OPENGL_2D_BASE
{
public:
	OPENGL_2D_BASE(FONT_NAME font_name
		, int font_size
		, FONT_STYLE font_style = STYLE_REGULAR
		, OPENGL_FONT_TYPE opengl_font_type = PIXMAP_FONT);
	~OPENGL_2D_BASE();

	OPENGL_VEC3 ProjecteVec3(OPENGL_VEC3 vec);

	FTFont*	GetFont() { return ft_font_; }
	OPENGL_FONT_TYPE GetOpenglFontType() { return openg_font_type_; }
	FONT_NAME GetFontName() { return font_name_; }
	FONT_STYLE GetFontStyle() { return font_style_; }
	void SetFontSize(int size);
	int GetfontSize() { return font_size_; }

protected:
	OPENGL_SIZE GetScreenSize();
	OPENGL_SIZE PushScreenCoordinateMatrix();
	void PopScreenCoordinateMatrix();

	FTFont				*ft_font_;
	OPENGL_FONT_TYPE	openg_font_type_;
	FONT_NAME			font_name_;
	FONT_STYLE			font_style_;
	int					font_size_;
};

//////////////////////////////////////////////////////////////////////////
// OPENGL_2D_LINE_TEXT
class OPENGL_2D_LINE_TEXT : public OPENGL_2D_BASE
{
public:
	OPENGL_2D_LINE_TEXT(FONT_NAME font_name
		, int font_size
		, FONT_STYLE font_style = STYLE_REGULAR
		, OPENGL_FONT_TYPE opengl_font_type = PIXMAP_FONT);
	~OPENGL_2D_LINE_TEXT();

	void PrintLineText(int pos_x, int pos_y, int num_line, OPENGL_COLOR text_color, const char *fmt, ...);
	void PrintLineText(int pos_x, int pos_y, int num_line, OPENGL_COLOR text_color, bool is_bg, OPENGL_COLOR bg_color, const char *fmt, ...);
};

//////////////////////////////////////////////////////////////////////////
// OPENGL_2D_BLOCK_TEXT
class OPENGL_2D_BLOCK_TEXT : public OPENGL_2D_BASE
{
public:
	OPENGL_2D_BLOCK_TEXT(FONT_NAME font_name
		, int font_size
		, FONT_STYLE font_style = STYLE_REGULAR
		, OPENGL_FONT_TYPE opengl_font_type = PIXMAP_FONT
		, TEXT_BLOCK_LAYOUT block_layout = TEXT_BLOCK_RIGHT);
	~OPENGL_2D_BLOCK_TEXT();

	void Init();
	void PrintLineText(int num_line, OPENGL_COLOR text_color, const char *fmt, ...);
	void PrintLineText(int num_line, OPENGL_COLOR text_color, bool is_bg, OPENGL_COLOR bg_color, const char *fmt, ...);

	void PrintBlockText(OPENGL_COLOR text_color, const char *fmt, ...);
	void PrintBlockText(OPENGL_COLOR text_color, bool is_bg, OPENGL_COLOR bg_color, const char *fmt, ...);

	void SetLeftRightLength(int lrLength);
private:
	void ChangeLayout();

	FTSimpleLayout		*font_layout_;
	OPENGL_POINT		pos_;
	OPENGL_SIZE			margin_;
	TEXT_BLOCK_LAYOUT	block_layout_;
	int					left_right_length_;
};

//////////////////////////////////////////////////////////////////////////
// OPENGL_2D_MENU_TEXT
class OPENGL_2D_MENU_TEXT : public OPENGL_2D_BASE
{
public:
	OPENGL_2D_MENU_TEXT(FONT_NAME font_name
		, int font_size
		, FONT_STYLE font_style = STYLE_REGULAR
		, OPENGL_FONT_TYPE opengl_font_type = PIXMAP_FONT);
	~OPENGL_2D_MENU_TEXT();

	void Init();
	void PrintMenu(const char* str, OPENGL_COLOR text_color, int hightlight_line, OPENGL_COLOR highlight_color, bool is_bg, OPENGL_COLOR bg_color);
	void PrintSubMenu(OPENGL_POINT pos, const char* str, OPENGL_COLOR text_color, int hightlight_line, OPENGL_COLOR highlight_color, bool is_bg, OPENGL_COLOR bg_color);
	void SetLineLength(int line_length);

	OPENGL_POINT GetSubMenuPos(int selected_line);

private:
	void ResetPosition();

	FTSimpleLayout		*font_layout_;
	OPENGL_POINT		pos_;
	OPENGL_SIZE			margin_;
	OPENGL_SIZE			offset_;
};