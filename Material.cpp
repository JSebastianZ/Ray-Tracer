#include "Material.h"

using namespace rt;

// Default constructor, with color black, and all attributes = 0.
Material::Material() : m_color{ 1,1,1 }, m_ambient{ 0.1 }, m_diffuse{ 0.9 }, m_specular{ 0.9 },
m_shininess{ 200.0 }, m_reflective{ 0 } {}

// Constructor of material object with specified reflection values and color.
Material::Material(Color color, real ambient, real diffuse, real specular, real shininess, real reflective) :
	m_color{ color }, m_ambient{ ambient }, m_diffuse{ diffuse }, m_specular{ specular }, 
	m_shininess{ shininess }, m_reflective{ reflective }{}

// Constructor with pattern.

// Constructor of material object with specified reflection values and color.
Material::Material(Color color, real ambient, real diffuse, real specular, real shininess, Pattern pattern) :
	m_color{ color }, m_ambient{ ambient }, m_diffuse{ diffuse }, m_specular{ specular },
	m_shininess{ shininess }, m_pattern{ pattern } {  m_pattern.m_pattern_exist = true; }