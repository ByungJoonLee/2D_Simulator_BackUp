// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>

// TODO: reference additional headers your program requires here
#include <boost/asio.hpp>
#include <windows.h>

#include <GL/glew.h>
#include <GL/glut.h>

#include "MACROS.h"
#include "COMMON_DEFINITION.h"
#include "ARRAY.h"
#include "ARRAY_2D.h"
#include "DYNAMIC_ARRAY.h"
#include "GRID_STRUCTURE_2D.h"


// #define USE_MANIPULATOR_TEST

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

#ifdef _DEBUG
#define DEBUG_NEW   new( _NORMAL_BLOCK, __FILE__, __LINE__)
#else
#define DEBUG_NEW
#endif