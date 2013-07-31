#ifndef	CD_MATH_2D
#define CD_MATH_2D

namespace cd{
	struct vector2{
		union{
			struct{
				double v0;
				double v1;
			};
			struct{
				double x;
				double y;
			};
			double a[2];
		};

		vector2():v0(0.0),v1(0.0){ }
		
		vector2(const double* val){
			for (size_t i = 0; i<2; i++)
				a[i] = val[i];
		}

		vector2(double val0, double val1){
			v0 = val0;
			v1 = val1;
		}

		void print() const{
			printf("{%e,%e}\n", v0, v1);
		}

		double square_magnitude() const{
			return v0*v0+v1*v1;
		}

		double magnitude() const{
			return sqrt(square_magnitude());
		}

		vector2 normalize() const{
			double s = magnitude();
			return vector2(v0/s, v1/s);
		}

		vector2 set_magnitude(double val) const{
			double s = val/magnitude();
			return vector2(v0*s, v1*s);
		}

	};

	inline vector2 operator - (const vector2 &l){
		return vector2(-l.v0, -l.v1);
	}

	inline vector2 operator + (const vector2 &l){
		return vector2(+l.v0, +l.v1);
	}

	inline vector2 operator + (const vector2 &l, const vector2 &r){
		return vector2(l.v0+r.v0, l.v1+r.v1);
	}

	inline vector2 operator - (const vector2 &l, const vector2 &r){
		return vector2(l.v0-r.v0, l.v1-r.v1);
	}

	inline vector2 operator * (const vector2 &l, const double &r){
		return vector2(l.v0*r, l.v1*r);
	}

	inline vector2 operator * (const double &l, const vector2 &r){
		return r*l;
	}

	inline vector2 operator / (const vector2 &l, const double &r){
		return vector2(l.v0/r, l.v1/r);
	}

    inline vector2 operator += (vector2 &l,const vector2 &r){
        l.x+=r.x;
        l.y+=r.y;
        return l;
    }

    inline vector2 abs_elem(const vector2 &v){
        return vector2(abs(v.x), abs(v.y));
    }


    ///////////
    //matrix
    //

    struct matrix2{
        //
        //  |0 2|
        //  |1 3|
        //
        
    public:
        double a[4];
        
        matrix2(){
            for (int i= 0; i<4; i++)
                a[i]= 0;
        }

        static matrix2 unit(){
            matrix2 m;
            for (int i= 0; i<2; i++)
                m.a[i*2+i] = 1.0;
            return m;
        }

        matrix2(double a0, double a1, double a2, double a3){
            a[0] = a0;
            a[1] = a1;
            a[2] = a2;
            a[3] = a3;
        }
        
        inline matrix2 transpose() const{
            return matrix2 (a[0],a[2],a[1],a[3]);
        }
        
        matrix2 inverse() const{
            double det = 1.0f/(a[0]*a[3]-a[1]*a[2]);
            return matrix2( a[3]*det, -a[1]*det,
                           -a[2]*det,  a[0]*det);
        }
    };

    inline vector2 operator * (const matrix2 &l,const vector2 &r){
        return vector2(r.x*l.a[0]+r.y*l.a[2],
                       r.x*l.a[1]+r.y*l.a[3]);
    }

    inline bool DcCrossPosition(
            const vector2 &x0, const vector2 &x1,
            const vector2 &y0, const vector2 &y1,
            vector2* result){
        
        const vector2 vx = x1-x0;
        const vector2 vy = y1-y0;

        const matrix2 m(
            vx.x, vy.x,
            vx.y, vy.y);

        const matrix2 m_inv = m.inverse();

        const vector2 x0d = m_inv*x0;
        const vector2 x1d = m_inv*x1;
        const vector2 y0d = m_inv*y0;
        const vector2 y1d = m_inv*y1;

        const double x = y0d.x;
        const double y = x0d.y;

        const bool condx = x0d.x<x && x<x1d.x;
        const bool condy = y0d.y<y && y<y1d.y;

        *result = m*vector2(x,y);
        return condx && condy;
    }
}

#endif
