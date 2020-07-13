#include "Canvas.h"

using namespace rt;

// Canvas constructor. The pixels array size is width * height
Canvas::Canvas(const int width, const int height) : m_width{ width }, m_height{ height }{
    pixels = new Color[m_height * m_width];
    for (int i = 0; i < (m_height * m_width); ++i) {
        pixels[i] = Color();
    }
}

// Draws RGB color values into a single canvas pixel.
void Canvas::write_pixel(int x, int y, Color color) {
    int i = (x * (m_width) + y);
   this-> pixels[i] = color;
}

// Returns the specified canvas pixel RGB color value.
Color Canvas::get_pixel(int x, int y) {
    return pixels[(x * (m_width) + y)];
}