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
	template <size_t n, bool inverse = false>
	struct complex_fft{
		complex_fft_n2<n, inverse> fft;
		complex_fft_n2<n, !inverse> ifft;
		complex work_buffer[1<<n];
		complex signal_buffer[1<<n];

		void execute(complex* data, size_t size){
			const size_t base_size = 1<<n;

			if (size>base_size/2){
				printf("faili\n");
				return;
			}

			auto w = [size](double x)->complex{
				const double theta = x*2*M_PI/size;
				if (inverse)
					return complex(cos(theta),-sin(theta));
				else
					return complex(cos(theta),sin(theta));
			};
			for (int i= 0; i<size; i++)
				signal_buffer[i]=data[i]*w(i*i/2.0);
			for (int i= size; i<base_size; i++)
				signal_buffer[i] = 0;

			for (int i= 0; i<base_size; i++)
				work_buffer[i] = 0;
			for (int i= 0; i<size; i++){
				work_buffer[i] = w(-i*i/2.0);
				if (i!=0)
				work_buffer[base_size-i]=w(-i*i/2.0);
			}

			fft.execute(work_buffer);
			fft.execute(signal_buffer);

			for (int i= 0; i<base_size; i++)
				work_buffer[i] = signal_buffer[i]*work_buffer[i];

			ifft.execute(work_buffer);

			for (int i= 0; i<base_size; i++)
				work_buffer[i] = work_buffer[i]*w(i*i/2.0);
			for (int i= 0; i<size; i++)
				data[i]=work_buffer[i];
		}
	};

	template <bool inverse = false>
	struct slow_fft{

		void execute(complex* data_in, complex* data_out, size_t size){
			auto w = [size](double x)->complex{
				const double theta = x*2*M_PI/size;
				if (inverse)
					return complex(cos(theta),-sin(theta));
				else
					return complex(cos(theta),sin(theta));
			};
			for (int i=0; i<size; i++){
				complex tmp = 0;
				for (int j=0; j<size; j++){
					tmp = tmp+data_in[j]*w(i*j);
				}
				data_out[i] = tmp;
			}
		}
		void execute(double* data_in, complex* data_out, size_t size){
			complex* c_data_in = new complex[size];
			for (int i=0; i<size; i++){
				c_data_in[i] = data_in[i];
			}
			execute(c_data_in, data_out, size);
			delete[] c_data_in;
		}
	};
	/*
	template <size_t n>
	struct complex_fft{
		complex_fft_n2<> base_fft;
		complex chirp[n];
		complex_fft(size_t size){
			for ()

		}

		void execute(complex* data){

		}

	};
	*/
}
