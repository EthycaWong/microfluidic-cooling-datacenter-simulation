// raygui implementation and its dark style live in C (the style header is not valid C++).
#define RAYGUI_IMPLEMENTATION
#include "raylib.h"
#include "raygui.h"
#include "styles/dark/style_dark.h"

void ViewerLoadDarkStyle(void) { GuiLoadStyleDark(); }
