
namespace cd{
	struct imaginary{
		double val;

		imaginary(){ }
		explicit imaginary(double value){
			val = value;	
		}
	};

	inline imaginary operator * (const imaginary &l, const double&r){
		return imaginary(l.val*r);
	}

	inline imaginary operator * (const double &l, const imaginary &r){
		return imaginary(r*l);
	}

	struct complex{
		double re;
		double im;

		complex(){ }
		complex(double re):re(re),im(0){}
		complex(double re, double im):re(re),im(im){}

		complex conjugate() const{
			return complex(re, -im);
		}
		void print() const{
			printf("%e+%ei\n", re, im);
		}
	};

	inline complex operator + (const double &l, const imaginary &r){
		return complex(l,r.val);
	}

	inline complex operator + (const imaginary &l, const double &r){
		return r+l;
	}

	inline complex operator - (const complex &l){
		return complex(-l.re, -l.im);
	}

	inline complex operator + (const complex &l, const complex &r){
		return complex(l.re+r.re, l.im+r.im);
	}

	inline complex& operator += (complex &l, const complex &r){
		return l = l+r;
	}

	inline complex operator - (const complex &l, const complex &r){
		return complex(l.re-r.re, l.im-r.im);
	}

	inline complex operator * (const complex &l, const complex &r){
		return complex(l.re*r.re-l.im*r.im, l.re*r.im+r.re*l.im);
	}

	inline complex operator * (const complex &l, const double&r){
		return complex(l.re*r, r*l.im);
	}

	inline complex operator * (const double &l, const complex &r){
		return r*l;
	}

	inline complex operator * (const complex &l, const imaginary &r){
		return complex(-l.im*r.val, l.re*r.val);
	}
	inline complex operator * (const imaginary &l, const complex &r){
		return r*l;
	}
}
