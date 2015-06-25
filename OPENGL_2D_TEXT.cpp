/********************************************************************
	OPENGL_2D_TEXT.cpp
*********************************************************************/

#include "OPENGL_2D_TEXT.h"
#include <iostream>
#include <string>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////////
// OPENGL_2D_BASE
OPENGL_2D_BASE::OPENGL_2D_BASE(FONT_NAME font_name
	, int font_size
	, FONT_STYLE font_style
	, OPENGL_FONT_TYPE opengl_font_type)
	: ft_font_(0)
	, font_name_(font_name)
	, font_style_(font_style)
	, font_size_(font_size)
	, openg_font_type_(opengl_font_type)
{
	// find folder dir
	char * pPath = 0;
	pPath = getenv("SystemRoot");
	std::string fontPath;
	if (pPath != 0)
		fontPath = pPath;
	else
		fontPath = "C:/Windows";
	fontPath += "/fonts/";

	// font name
	std::string ttf_name;
	switch (font_name_)
	{
	case FONT_ARIAL:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "arial.ttf";		break;
			case STYLE_BOLD:		ttf_name = "arialdb.ttf";	break;
			case STYLE_ITALIC:		ttf_name = "ariali.ttf";	break;
			case STYLE_BOLD_ITALIC: ttf_name = "arialbi.ttf";	break;
			}
		}
		break;

	case FONT_CONSOLAS:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "Consola.ttf";		break;
			case STYLE_BOLD:		ttf_name = "Consolab.ttf";		break;
			case STYLE_ITALIC:		ttf_name = "Consolai.ttf";		break;
			case STYLE_BOLD_ITALIC:	ttf_name = "Consolaz.ttf";		break;
			}
		}
		break;

	case FONT_COURIER_NEW:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "cour.ttf";		break;
			case STYLE_BOLD:		ttf_name = "courb.ttf";		break;
			case STYLE_ITALIC:		ttf_name = "couri.ttf";		break;
			case STYLE_BOLD_ITALIC:	ttf_name = "courbi.ttf";	break;
			}
		}
		break;

	case FONT_GEORGIA:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "georgia.ttf";		break;
			case STYLE_BOLD:		ttf_name = "georgiab.ttf";		break;
			case STYLE_ITALIC:		ttf_name = "georgiai.ttf";		break;
			case STYLE_BOLD_ITALIC:	ttf_name = "georgiaz.ttf";		break;
			}
		}
		break;

	case FONT_IMPACT:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "impact.ttf";		break;
			case STYLE_BOLD:		ttf_name = "impact.ttf";		break;
			case STYLE_ITALIC:		ttf_name = "impact.ttf";		break;
			case STYLE_BOLD_ITALIC:	ttf_name = "impact.ttf";		break;
			}
		}
		break;

	case FONT_PLATINO:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "pala.ttf";		break;
			case STYLE_BOLD:		ttf_name = "palab.ttf";		break;
			case STYLE_ITALIC:		ttf_name = "palai.ttf";		break;
			case STYLE_BOLD_ITALIC:	ttf_name = "palabi.ttf";	break;
			}
		}
		break;

	case FONT_TAHOMA:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "tahoma.ttf";		break;
			case STYLE_BOLD:		ttf_name = "tahomabd.ttf";		break;
			case STYLE_ITALIC:		ttf_name = "tahoma.ttf";		break;
			case STYLE_BOLD_ITALIC:	ttf_name = "tahomabd.ttf";		break;
			}
		}
		break;

	case FONT_TIMES_NEW_ROMAN:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "times.ttf";		break;
			case STYLE_BOLD:		ttf_name = "timesb.ttf";	break;
			case STYLE_ITALIC:		ttf_name = "timesi.ttf";	break;
			case STYLE_BOLD_ITALIC:	ttf_name = "timesbi.ttf";	break;
			}
		}
		break;

	case FONT_VERDANA:
		{
			switch (font_style) 
			{
			case STYLE_REGULAR:		ttf_name = "verdana.ttf";		break;
			case STYLE_BOLD:		ttf_name = "verdanab.ttf";		break;
			case STYLE_ITALIC:		ttf_name = "verdanai.ttf";		break;
			case STYLE_BOLD_ITALIC:	ttf_name = "verdanaz.ttf";		break;
			}
		}
		break;
	}
	fontPath += ttf_name;

	// opengl font type
	switch(openg_font_type_)
	{
	case BITMAP_FONT: ft_font_ = new FTBitmapFont(fontPath.c_str());	break;
	case PIXMAP_FONT: ft_font_ = new FTPixmapFont(fontPath.c_str());	break;
	default: ft_font_ = new FTPixmapFont(fontPath.c_str());	break;
	}

	SetFontSize(font_size_);

	if(ft_font_->Error())
	{
		std::cout << "Failed to open font " << fontPath << std::endl;
		delete ft_font_;
		ft_font_ = 0;
	}
}

OPENGL_2D_BASE::~OPENGL_2D_BASE()
{
	if(ft_font_)
		delete ft_font_;
}

OPENGL_SIZE OPENGL_2D_BASE::GetScreenSize()
{
	GLint	viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	return OPENGL_SIZE(viewport[2] - viewport[0], viewport[3] - viewport[1]);
}

OPENGL_SIZE OPENGL_2D_BASE::PushScreenCoordinateMatrix() 
{
	GLint	viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(viewport[0],viewport[2],viewport[1],viewport[3]);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	return OPENGL_SIZE(viewport[2] - viewport[0], viewport[3] - viewport[1]);
}

void OPENGL_2D_BASE::PopScreenCoordinateMatrix() 
{
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void OPENGL_2D_BASE::SetFontSize(int size)
{
	if(size < 2)
		size = 2;

	font_size_ = size;

	if(ft_font_)
		ft_font_->FaceSize(font_size_);
}

OPENGL_VEC3 OPENGL_2D_BASE::ProjecteVec3(OPENGL_VEC3 vec)
{
	GLdouble x,y,z;
	GLint viewport[4];
	GLdouble modelview_mat[16] = {0.0, };
	GLdouble projection_mat[16] = {0.0, };
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_mat);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_mat);

	gluProject(vec.GetX(), vec.GetY(), vec.GetZ(), modelview_mat, projection_mat, viewport,  &x, &y, &z);

	return OPENGL_VEC3(x,y,z);
}

//////////////////////////////////////////////////////////////////////////
// OPENGL_2D_LINE_TEXT
OPENGL_2D_LINE_TEXT::OPENGL_2D_LINE_TEXT(FONT_NAME font_name
							, int font_size
							, FONT_STYLE font_style
							, OPENGL_FONT_TYPE opengl_font_type)
	: OPENGL_2D_BASE(font_name, font_size, font_style, opengl_font_type)
{
}

OPENGL_2D_LINE_TEXT::~OPENGL_2D_LINE_TEXT()
{
}

void OPENGL_2D_LINE_TEXT::PrintLineText(int pos_x, int pos_y, int num_line, OPENGL_COLOR text_color, const char *fmt, ...)
{
	PushScreenCoordinateMatrix();

	char str[256] = {0,};
	va_list ap;
	if (fmt == NULL)
		*str=0;
	else 
	{
		va_start(ap, fmt);
		vsprintf(str, fmt, ap);
		va_end(ap);
	}

	glColor3f(text_color.GetRed(), text_color.GetGreen(), text_color.GetBlue());
	glRasterPos2i(pos_x, (pos_y - GetFont()->Descender()) - (GetFont()->LineHeight() * num_line));
	GetFont()->Render(str);

	PopScreenCoordinateMatrix();
}

void OPENGL_2D_LINE_TEXT::PrintLineText(int pos_x, int pos_y, int num_line, OPENGL_COLOR text_color, bool is_bg, OPENGL_COLOR bg_color, const char *fmt, ...)
{
	PushScreenCoordinateMatrix();

	char str[256] = {0,};
	va_list ap;
	if (fmt == NULL)
		*str=0;
	else 
	{
		va_start(ap, fmt);
		vsprintf(str, fmt, ap);
		va_end(ap);
	}

	static const int box_offset_x = 2;
	static const int box_offset_y = 2;

	if(is_bg)
	{
		glColor3f(bg_color.GetRed(), bg_color.GetGreen(), bg_color.GetBlue());

		FTBBox bb = ft_font_->BBox(str);
		float x_up = bb.Upper().Xf() + box_offset_x;
		float x_lo = bb.Lower().Xf() - box_offset_x;
		float y_up = bb.Upper().Yf() + box_offset_y; 
		float y_lo = bb.Lower().Yf() - box_offset_y; 

		glPushMatrix();
		glTranslatef(pos_x, (pos_y - GetFont()->Descender()) - (GetFont()->LineHeight() * num_line), 0);
		glBegin(GL_QUADS);
		glVertex2f(x_lo, y_lo);
		glVertex2f(x_up, y_lo);
		glVertex2f(x_up, y_up);
		glVertex2f(x_lo, y_up);
		glEnd();
		glPopMatrix();
	}

	glColor3f(text_color.GetRed(), text_color.GetGreen(), text_color.GetBlue());
	glRasterPos2i(pos_x, (pos_y - GetFont()->Descender()) - (GetFont()->LineHeight() * num_line));
	GetFont()->Render(str);

	PopScreenCoordinateMatrix();
}

//////////////////////////////////////////////////////////////////////////
// OPENGL_2D_BLOCK_TEXT
#define LEFT_RIGHT_LENGTH 130
OPENGL_2D_BLOCK_TEXT::OPENGL_2D_BLOCK_TEXT(FONT_NAME font_name
									, int font_size
									, FONT_STYLE font_style
									, OPENGL_FONT_TYPE opengl_font_type
									, TEXT_BLOCK_LAYOUT block_layout)
	: OPENGL_2D_BASE(font_name, font_size, font_style, opengl_font_type)
	, font_layout_(0)
	, pos_(0,0)
	, margin_(10,10)
	, block_layout_(block_layout)
	, left_right_length_(LEFT_RIGHT_LENGTH)
{
	font_layout_ = new FTSimpleLayout;
	font_layout_->SetFont(GetFont());
	font_layout_->SetLineLength(left_right_length_);
}

OPENGL_2D_BLOCK_TEXT::~OPENGL_2D_BLOCK_TEXT()
{
	if(font_layout_)
		delete font_layout_;
}

void OPENGL_2D_BLOCK_TEXT::SetLeftRightLength(int lrLength)
{
	left_right_length_ = lrLength;
	ChangeLayout();
}

void OPENGL_2D_BLOCK_TEXT::Init()
{
	ChangeLayout();
}

void OPENGL_2D_BLOCK_TEXT::ChangeLayout()
{
	OPENGL_SIZE screen_size = GetScreenSize();

	switch(block_layout_)
	{
	case TEXT_BLOCK_TOP: 
		{
			pos_.SetX(margin_.Width());
			pos_.SetY(screen_size.Height() - margin_.Height() - GetFont()->LineHeight());
			font_layout_->SetLineLength(screen_size.Width() - (margin_.Width()*2));
			font_layout_->SetAlignment(FTGL::ALIGN_LEFT);
		}
		break;
	case TEXT_BLOCK_LEFT: 
		{
			pos_.SetX(margin_.Width());
			pos_.SetY(screen_size.Height() - margin_.Height() - GetFont()->LineHeight());
			font_layout_->SetLineLength(left_right_length_);
			font_layout_->SetAlignment(FTGL::ALIGN_LEFT);
		}
		break;
	case TEXT_BLOCK_RIGHT:
		{
			pos_.SetX(screen_size.Width() - left_right_length_ - margin_.Width());
			pos_.SetY(screen_size.Height() - margin_.Height() - GetFont()->LineHeight());
			font_layout_->SetLineLength(left_right_length_);
			font_layout_->SetAlignment(FTGL::ALIGN_LEFT);
		}
		break;
	case TEXT_BLOCK_BOTTOM:
		{
			pos_.SetX(margin_.Width());
			pos_.SetY(margin_.Height());
			font_layout_->SetLineLength(screen_size.Width() - (margin_.Width()*2));
			font_layout_->SetAlignment(FTGL::ALIGN_LEFT);
		}
		break;
	}
}

void OPENGL_2D_BLOCK_TEXT::PrintLineText(int num_line, OPENGL_COLOR text_color, const char *fmt, ...)
{
	ChangeLayout();

	PushScreenCoordinateMatrix();

	char str[256] = {0,};
	va_list ap;
	if (fmt == NULL)
		*str=0;
	else 
	{
		va_start(ap, fmt);
		vsprintf(str, fmt, ap);
		va_end(ap);
	}

	glColor3f(text_color.GetRed(), text_color.GetGreen(), text_color.GetBlue());
	switch(block_layout_)
	{
	case TEXT_BLOCK_TOP:		
		glRasterPos2i(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line));
		font_layout_->Render(str, -1, true);
		break;
	case TEXT_BLOCK_LEFT: 
		glRasterPos2i(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line));
		font_layout_->Render(str, -1, true);
		break;
	case TEXT_BLOCK_RIGHT:
		glRasterPos2i(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line));
		font_layout_->Render(str, -1, true);
		break;
	case TEXT_BLOCK_BOTTOM:
		glRasterPos2i(pos_.X(), (pos_.Y()- GetFont()->Descender()) + (GetFont()->LineHeight() * num_line));
		font_layout_->Render(str, -1, true);
		break;
	}

	PopScreenCoordinateMatrix();
}

void OPENGL_2D_BLOCK_TEXT::PrintLineText(int num_line, OPENGL_COLOR text_color, bool is_bg, OPENGL_COLOR bg_color, const char *fmt, ...)
{
	ChangeLayout();

	PushScreenCoordinateMatrix();

	char str[256] = {0,};
	va_list ap;
	if (fmt == NULL)
		*str=0;
	else 
	{
		va_start(ap, fmt);
		vsprintf(str, fmt, ap);
		va_end(ap);
	}

	static int box_offset_x = 2;
	static int box_offset_y = 2;

	if(is_bg)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glColor4f(bg_color.GetRed(), bg_color.GetGreen(), bg_color.GetBlue(), 0.5f);

		FTBBox bb = ft_font_->BBox(str);
		float x_up = bb.Upper().Xf() + box_offset_x;
		float x_lo = bb.Lower().Xf() - box_offset_x;
		float y_up = bb.Upper().Yf() + box_offset_y; 
		float y_lo = bb.Lower().Yf() - box_offset_y; 

		glPushMatrix();
		switch(block_layout_)
		{
		case TEXT_BLOCK_TOP:		
			glTranslatef(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line), 0);
			break;
		case TEXT_BLOCK_LEFT: 
			glTranslatef(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line), 0);
			break;
		case TEXT_BLOCK_RIGHT:
			glTranslatef(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line), 0);
			break;
		case TEXT_BLOCK_BOTTOM:
			glTranslatef(pos_.X(), (pos_.Y()- GetFont()->Descender()) + (GetFont()->LineHeight() * num_line), 0);
			break;
		}

		glBegin(GL_QUADS);
		glVertex2f(x_lo, y_lo);
		glVertex2f(x_up, y_lo);
		glVertex2f(x_up, y_up);
		glVertex2f(x_lo, y_up);
		glEnd();
		glPopMatrix();

		glDisable(GL_BLEND);
	}

	glColor3f(text_color.GetRed(), text_color.GetGreen(), text_color.GetBlue());
	switch(block_layout_)
	{
	case TEXT_BLOCK_TOP:		
		glRasterPos2i(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line));
		font_layout_->Render(str, -1, true);
		break;
	case TEXT_BLOCK_LEFT: 
		glRasterPos2i(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line));
		font_layout_->Render(str, -1, true);
		break;
	case TEXT_BLOCK_RIGHT:
		glRasterPos2i(pos_.X(), (pos_.Y()- GetFont()->Descender()) - (GetFont()->LineHeight() * num_line));
		font_layout_->Render(str, -1, true);
		break;
	case TEXT_BLOCK_BOTTOM:
		glRasterPos2i(pos_.X(), (pos_.Y()- GetFont()->Descender()) + (GetFont()->LineHeight() * num_line));
		font_layout_->Render(str, -1, true);
		break;
	}

	PopScreenCoordinateMatrix();
}

void OPENGL_2D_BLOCK_TEXT::PrintBlockText(OPENGL_COLOR text_color, const char *fmt, ...)
{
	char str[256] = {0,};
	va_list ap;
	if (fmt == NULL)
		*str=0;
	else 
	{
		va_start(ap, fmt);
		vsprintf(str, fmt, ap);
		va_end(ap);
	}

	PrintLineText(0, text_color, str);
}

void OPENGL_2D_BLOCK_TEXT::PrintBlockText(OPENGL_COLOR text_color, bool is_bg, OPENGL_COLOR bg_color, const char *fmt, ...)
{
	char str[256] = {0,};
	va_list ap;
	if (fmt == NULL)
		*str=0;
	else 
	{
		va_start(ap, fmt);
		vsprintf(str, fmt, ap);
		va_end(ap);
	}

	static const int offset_x = 5;
	static const int offset_y = 5;
	if(is_bg)
	{
		ChangeLayout();

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		PushScreenCoordinateMatrix();

		glColor4f(bg_color.GetRed(), bg_color.GetGreen(), bg_color.GetBlue(), 0.5f);
		float line_width = font_layout_->GetLineLength();
		
		FTBBox bb = font_layout_->BBox(str);
		float x_up = line_width + offset_x;
		float x_lo = bb.Lower().Xf() - offset_x;
		float y_up = bb.Upper().Yf() + offset_y;
		float y_lo = bb.Lower().Yf() - offset_y;

		if(font_layout_->GetAlignment() == FTGL::ALIGN_RIGHT)
		{
			x_up = bb.Upper().Xf() + offset_x;
			x_lo = bb.Upper().Xf() - line_width - offset_y;
		}

		glPushMatrix();
		glTranslatef(pos_.X(), pos_.Y()- GetFont()->Descender(), 0);

		glBegin(GL_QUADS);
		glVertex2f(x_lo, y_lo);
		glVertex2f(x_up, y_lo);
		glVertex2f(x_up, y_up);
		glVertex2f(x_lo, y_up);
		glEnd();

		glColor3f(text_color.GetRed(), text_color.GetGreen(), text_color.GetBlue());
		glBegin(GL_LINE_LOOP);
		glVertex2f(x_lo, y_lo);
		glVertex2f(x_up, y_lo);
		glVertex2f(x_up, y_up);
		glVertex2f(x_lo, y_up);
		glEnd();

		glPopMatrix();
		PopScreenCoordinateMatrix();

		glDisable(GL_BLEND);
	}

	PrintLineText(0, text_color, str);
}

//////////////////////////////////////////////////////////////////////////
// OPENGL_2D_MENU_TEXT
#define MENU_LENGTH 250
OPENGL_2D_MENU_TEXT::OPENGL_2D_MENU_TEXT(FONT_NAME font_name
			, int font_size
			, FONT_STYLE font_style
			, OPENGL_FONT_TYPE opengl_font_type)
	: OPENGL_2D_BASE(font_name, font_size, font_style, opengl_font_type)
	, font_layout_(0)
	, pos_(0,0)
	, margin_(10,10)
	, offset_(5,5)
{
	font_layout_ = new FTSimpleLayout;
	font_layout_->SetFont(GetFont());
	font_layout_->SetLineLength(MENU_LENGTH);
}

OPENGL_2D_MENU_TEXT::~OPENGL_2D_MENU_TEXT()
{
	if(font_layout_)
		delete font_layout_;
}

void OPENGL_2D_MENU_TEXT::Init()
{
	ResetPosition();
}

void OPENGL_2D_MENU_TEXT::SetLineLength(int line_length)
{
	font_layout_->SetLineLength(line_length);
}

void OPENGL_2D_MENU_TEXT::ResetPosition()
{
	OPENGL_SIZE screen_size = GetScreenSize();

	pos_.SetX(margin_.Width());
	pos_.SetY(screen_size.Height() - margin_.Height() - GetFont()->Ascender());
}

OPENGL_POINT OPENGL_2D_MENU_TEXT::GetSubMenuPos(int selected_line)
{
	ResetPosition();

	float x_up = font_layout_->GetLineLength() + margin_.Width();
	float y_up = -(selected_line  * ft_font_->LineHeight());

	return OPENGL_POINT(x_up + pos_.X(), y_up + pos_.Y());
}

void OPENGL_2D_MENU_TEXT::PrintMenu(const char* str, OPENGL_COLOR text_color, int hightlight_line, OPENGL_COLOR highlight_color, bool is_bg, OPENGL_COLOR bg_color)
{
	ResetPosition();
	PrintSubMenu(pos_, str, text_color, hightlight_line, highlight_color, is_bg, bg_color);
}

void OPENGL_2D_MENU_TEXT::PrintSubMenu(OPENGL_POINT pos, const char* str, OPENGL_COLOR text_color, int hightlight_line, OPENGL_COLOR highlight_color, bool is_bg, OPENGL_COLOR bg_color)
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// BG
	if(is_bg)
	{
		PushScreenCoordinateMatrix();

		glColor4f(bg_color.GetRed(), bg_color.GetGreen(), bg_color.GetBlue(), 0.5f);
		float line_width = font_layout_->GetLineLength();

		FTBBox bb = font_layout_->BBox(str);
		float x_up = line_width + offset_.Width();
		float x_lo = -offset_.Width();
		float y_up = bb.Upper().Yf() + offset_.Height();
		float y_lo = bb.Lower().Yf() - offset_.Height();

		glPushMatrix();
		glTranslatef(pos.X(), pos.Y(), 0);

		glBegin(GL_QUADS);
		glVertex2f(x_lo, y_lo);
		glVertex2f(x_up, y_lo);
		glVertex2f(x_up, y_up);
		glVertex2f(x_lo, y_up);
		glEnd();

		glColor4f(text_color.GetRed(), text_color.GetGreen(), text_color.GetBlue(), 0.5f);
		glBegin(GL_LINE_LOOP);
		glVertex2f(x_lo, y_lo);
		glVertex2f(x_up, y_lo);
		glVertex2f(x_up, y_up);
		glVertex2f(x_lo, y_up);
		glEnd();

		glPopMatrix();
		PopScreenCoordinateMatrix();
	}

	// Highlight menu
	PushScreenCoordinateMatrix();

	FTBBox bbox = font_layout_->BBox(str);
	float h_line = ft_font_->LineHeight();
	float w_line = font_layout_->GetLineLength();
	float x_lo = -offset_.Width();
	float x_up = w_line + offset_.Width();;
	float y_lo = bbox.Lower().Yf(); 
	float y_up = bbox.Upper().Yf(); 
	y_lo = -offset_.Height() - ((hightlight_line + 1) * h_line);
	y_up = -offset_.Height() - (hightlight_line  * h_line);

	glPushMatrix();
	glTranslatef(pos.X(), pos.Y() + ft_font_->LineHeight(), 0);

	glColor4f(highlight_color.GetRed(), highlight_color.GetGreen(), highlight_color.GetBlue(), 0.5f);
	// Draw the front face
	glBegin(GL_QUADS);
	glVertex3f(x_lo, y_lo, 0.0f);
	glVertex3f(x_up, y_lo, 0.0f);
	glVertex3f(x_up, y_up, 0.0f);
	glVertex3f(x_lo, y_up, 0.0f);
	glEnd();
	glPopMatrix();

	PopScreenCoordinateMatrix();

	glDisable(GL_BLEND);

	// Print menu
	PushScreenCoordinateMatrix();

	glColor3f(text_color.GetRed(), text_color.GetGreen(), text_color.GetBlue());
	glRasterPos2i(pos.X(), pos.Y());
	font_layout_->Render(str, -1, true);

	PopScreenCoordinateMatrix();
}