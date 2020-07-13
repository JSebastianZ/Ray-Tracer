#include "Tuple.h"

using namespace rt;

	// Holds the value along the x axis.
	real m_x = 0;
	
	// Holds the value along the y axis.
	real m_y = 0;

	// Holds the value along the z axis.
	real m_z = 0;

	// Indicates primitive type. 0 = vector / 1 = point.
	int m_w = 0;

	// Default constructor.
	Tuple::Tuple() : m_x(0), m_y(0), m_z(0), m_w(0) {}

	// Constructs a tuple according to parameters passed.
	Tuple::Tuple(real x, real y, real z, int w)
		: m_x{ x }, m_y{ y }, m_z{ z }, m_w{ w }{}

	// Destructor.
	Tuple::~Tuple() {}
