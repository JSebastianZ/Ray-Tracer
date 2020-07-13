#include "Light.h"

using namespace rt;


// Default constructor with light source at the origin and no illumination.
Light::Light() : m_source{ 0,0,0 }, m_intensity{ 0,0,0 }{}

// Light object constructor with specified position and color.
Light::Light(Point source, Color intensity) : m_source{ source }, m_intensity{ intensity }{}