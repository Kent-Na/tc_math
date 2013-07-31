#ifndef	CD_MATH_3D
#define CD_MATH_3D

namespace cd{
	struct vector3{
		union{
			struct{
				double v0;
				double v1;
				double v2;
			};
			struct{
				double x;
				double y;
				double z;
			};
			double a[3];
		};

		vector3():v0(0.0),v1(0.0),v2(0.0){ }
		
		vector3(const double* val){
			for (size_t i = 0; i<3; i++)
				a[i] = val[i];
		}

		vector3(double val0, double val1, double val2){
			v0 = val0;
			v1 = val1;
			v2 = val2;
		}

		void print() const{
			printf("{%e,%e,%e}\n", v0, v1, v2);
		}

		double square_magnitude() const{
			return v0*v0+v1*v1+v2*v2;
		}

		double magnitude() const{
			return sqrt(square_magnitude());
		}

		vector3 normalize() const{
			double s = magnitude();
			return vector3(v0/s, v1/s, v2/s);
		}

		vector3 set_magnitude(double val) const{
			double s = val/magnitude();
			return vector3(v0*s, v1*s, v2*s);
		}

	};

	inline vector3 operator - (const vector3 &l){
		return vector3(-l.v0, -l.v1, -l.v2);
	}

	inline vector3 operator + (const vector3 &l){
		return vector3(+l.v0, +l.v1, +l.v2);
	}

	inline vector3 operator * (const vector3 &l, const double &r){
		return vector3(l.v0*r, l.v1*r, l.v2*r);
	}

	inline vector3 operator / (const vector3 &l, const double &r){
		return vector3(l.v0/r, l.v1/r, l.v2/r);
	}


	
	struct matrix4{
		double a[16];
		
		// 0 4 8 c 
		// 1 5 9 d
		// 2 6 a e
		// 3 7 b f

		static matrix4 unit(){
			matrix4 out;
			for (size_t i=0; i<4; i++)
				out.a[i*4+i]=1.0;
			return out;
		}

		static matrix4 scale(double val){
			matrix4 out;
			for (size_t i=0; i<3; i++)
				out.a[i*4+i]=val;
			out.a[3*4+3]=1.0;
			return out;
		}

		static matrix4 translate(vector3 offset){
			matrix4 out = unit();
			for (size_t i = 0; i<3; i++){
				out.a[4*3+i] = offset.a[i];
			}
			return out;
		}
		matrix4(){
			for (size_t i = 0; i<16; i++)
				a[i] = 0.0;
		}
	};

	inline matrix4 operator * (const matrix4 &l, const matrix4 &r){
		matrix4 out;
		for (size_t i = 0; i<4; i++){
			for (size_t j=0; j<4; j++){
				for(size_t k=0; k<4; k++){
					out.a[i*4+j] += l.a[k*4+j]*r.a[i*4+k];
				}
			}
		}
		return out;
	}

	inline vector3 operator * (const matrix4 &l, const vector3 &r){
		//assume r4 = 1
		vector3 out;
		for (size_t i = 0; i<3; i++){
			for (size_t j = 0; j<3; j++){
				out.a[i] += l.a[j*4+i]*r.a[j];
			}
			out.a[i] += l.a[3*4+i];
		}
		return out;
	}

	struct quatanion{
		double n;
		vector3 v;
	
		quatanion():n(0.0) { }

		quatanion(double val_re, const vector3& val_im):
			n(val_re), v(val_im) { }

		quatanion(double val_re, double val0, double val1, double val2){
			n = val_re;
			v.v0 = val0;
			v.v1 = val1;
			v.v2 = val2;
		}

		static quatanion rotation(const vector3 &axis, const double angle){
			return quatanion(cos(angle*0.5),
					axis.set_magnitude(sin(angle*0.5)));
		}

		vector3 rotate(const vector3 &l){
			const double x = v.x;	
			const double y = v.y;	
			const double z = v.z;	

			const double nn=n*n;
			const double xx=x*x;
			const double yy=y*y;
			const double zz=z*z;
			return vector3(
				2*(x*y*l.y +x*z*l.z +n*y*l.z -n*z*l.y)
				+l.x*(nn +xx -yy -zz) ,
				2*(y*z*l.z +y*x*l.x +n*z*l.x -n*x*l.z)
				+l.y*(nn +yy -zz -xx) ,
				2*(z*x*l.x +z*y*l.y +n*x*l.y -n*y*l.x)
				+l.z*(nn +zz -xx -yy)
			);
		}

		vector3 rotate(double val0, double val1, double val2){
			return rotate(vector3(val0, val1, val2));
		}

		matrix4 rotate(const matrix4 &l){
			matrix4 out;
			for (size_t i = 0; i<4; i++){
				const vector3 tmp = rotate(vector3(l.a+4*i));
				for (int j = 0; j<3; j++)
					out.a[4*i+j] = tmp.a[j];
			}
			out.a[15] = 1.0;
			return out;
		}
	};

}

#endif
