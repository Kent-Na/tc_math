#ifndef	CD_MATH_ND
#define CD_MATH_ND

namespace cd{
	template <size_t n>
	struct vector{
		double a[n];

		vector(){
			for (int i=0; i<n; i++)
				a[i]=0.0;
		}
		vector(const double* val){
			for (int i=0; i<n; i++)
				a[i]=val[i];
		}

		/*
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
		*/
	};
	
	/*
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
	*/

	template<size_t n>
	struct matrix{
		double a[n*n];
	
		//m(i,j)

		// +-> j
		// |
		// i

		static matrix<n> unit(){
			matrix<n> out;
			for (size_t i=0; i<n; i++)
				out(i,i)=1.0;
			return out;
		}

		matrix(){
			for (size_t i = 0; i<n*n; i++)
				a[i] = 0.0;
		}
		matrix(double* val){
			for (size_t i = 0; i<n*n; i++)
				a[i] = val[i];
		}

		double& operator ()(size_t i, size_t j){
			return a[j*n+i];
		}
		const double& operator ()(size_t i, size_t j)const{
			return a[j*n+i];
		}
		matrix<n> inv(){
			return inverse();
		}
		matrix<n> inverse(){
			matrix<n> m = *this;	
			matrix<n> out = unit();
			
			for (int i = 0; i<n; i++){
				const double s = m(i,i);
				for (int j = 0; j<n; j++){
					m(i,j)  /= s;
					out(i,j)/= s;
				}

				for (int j = 0; j<n; j++){
					if (j==i) continue;
					const double s = m(j,i);
					for (int k = i; k<n; k++){
						m(j,k) -= s*m(i,k);
					}
					for (int k = 0; k<n; k++){
						out(j,k)-= s*out(i,k);
					}
				}
			}
			return out;
		}
	};

	template<size_t n>
	inline matrix<n> operator * (const matrix<n> &l, const matrix<n> &r){
		matrix<n> out;
		for (size_t i = 0; i<n; i++){
			for (size_t j=0; j<n; j++){
				for(size_t k=0; k<n; k++){
					out(i,j) += l(i,k)*r(k,j);
				}
			}
		}
		return out;
	}

	template<size_t n>
	inline vector<n> operator * (const matrix<n> &l, const vector<n> &r){
		vector<n> out;
		for (size_t i = 0; i<n; i++){
			for (size_t j = 0; j<n; j++){
				out.a[i] += l.a[j*n+i]*r.a[j];
			}
		}
		return out;
	}

}

#endif
