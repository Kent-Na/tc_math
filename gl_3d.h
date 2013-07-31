#include "math_2d.h"
#include "math_3d.h"

namespace cd{
	inline void vertex(const vector3 &v){
		glVertex3dv(v.a);
	}
	inline void vertex(const vector2 &v){
		glVertex2dv(v.a);
	}
}
