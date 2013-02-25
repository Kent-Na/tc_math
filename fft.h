namespace cd{
	template <size_t n, bool inverse = false>
	struct complex_fft_n2{
		complex *table;
		complex_fft_n2(){
			const size_t size = 1<<n;
			table = new complex[size];
			for (size_t i=0; i<size; i++){
				const double theta = i*2*M_PI/size;
				if (inverse)
					table[i] = complex(cos(theta),-sin(theta));
				else
					table[i] = complex(cos(theta),sin(theta));
			}
		}
		~complex_fft_n2(){
			delete[] table;
		}
		
		void execute(complex* data){
			const size_t size = 1<<n;

			for (size_t k=0; k<n; k++){
				const size_t offset = size>>(k+1);
				const size_t sub_size = 1<<k;
				for (size_t j=0; j<sub_size; j++){
					const size_t j_offset = j*offset*2;

					//i=0
					{
						const size_t i_dush = j_offset;
						const size_t i_offset = i_dush+offset;
						const complex tmp = data[i_dush]+data[i_offset];
						data[i_offset] = data[i_dush]-data[i_offset];
						data[i_dush] = tmp;
					}

					for (size_t i=1; i<offset; i++){
						const size_t i_dush = i+j_offset;
						const size_t i_offset = i_dush+offset;
						const complex tmp = data[i_dush]+data[i_offset];
						data[i_offset] = 
							(data[i_dush]-data[i_offset])*table[i<<k];
						data[i_dush] = tmp;
					}
				}
			}

			{
				size_t i = 0;
				for (size_t j=1; j<size-1; j++){
					for (size_t k=size>>1; k>(i^=k); k>>=1);
					if (j<i){
						const complex tmp = data[j];
						data[j] = data[i];
						data[i] = tmp;
					}
				}
			}
		}
	};

	template <size_t n, bool inverse = false>
	struct real_fft_n2{
		complex_fft_n2<n-1, inverse> base;
		
		//Result will be overwriten in *data 
		void execute(double* data){
			complex* c_data = (complex*)data;
			if (!inverse)
				base.execute(c_data);
			const size_t size = 1<<n;
			const size_t half_size = size/2;
			const size_t quarter_size = size/4;

			auto w = [size](double x)->complex{
				const double theta = x*2*M_PI/size;
				if (inverse)
					return complex(cos(theta),-sin(theta));
				else
					return complex(cos(theta),sin(theta));
			};

			//i == 0
			{
				const complex c0 = c_data[0];
				c_data[0].re = c0.re+c0.im;
				if (inverse)
					c_data[0].im = -c0.re+c0.im;
				else
					c_data[0].im = c0.re-c0.im;
			}

			for (size_t i=1; i<quarter_size; i++){
				complex c0 = c_data[i];
				complex c1 = c_data[half_size-i];
				c_data[i] = 
					c0-(0.5+imaginary(0.5)*w(i))*(c0-c1.conjugate());
				c_data[half_size-i] = 
					c1-(0.5+imaginary(0.5)*w(half_size-i))*
					(c1-c0.conjugate());
			}

			if (inverse){
				base.execute(c_data);
				for (size_t i=0; i<half_size; i++)
					c_data[i] = c_data[i].conjugate();
			}
		}
	};
}
