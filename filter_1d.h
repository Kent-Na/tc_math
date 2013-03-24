namespace cd{

	//overwrite result to data1
	template<typename T>
	void convolution(T *data1, T *data2, size_t size){
		fft(data1, size);
		fft(data2, size);
		for (size_t i = 0; i<size; i++)
			data1[i]*=data2[i];
		ifft(data1, size);
	}

	template<typename T>
	void convolution(T *data1, 
			std::unary_function<size_t, T> f2, size_t size){
		T *data2 = new T[size];
		for (size_t i = 0; i<size; i++)
			data2[i] = f2(i);
		convolution<>(data1, data2, size);
		delete[] data2;
	}

	void real_part(complex* data, size_t size){
		double *r_data = (double*)data;
		for (size_t i = 0; i<size; i++){
			r_data[i]=data[i].re;
		}
	}

	void magnitude(complex* data, size_t size){
		double *r_data = (double*)data;
		for (size_t i = 0; i<size; i++){
			r_data[i]=sqrt(data[i].re*data[i].re*
					data[i].im*data[i].im);
		}
	}

	double gaussian(double x, double a){
		return sqrt(a/M_PI)*pow(M_E, -a*x*x);
	}

	template <typename T>
	void block_lpf(T *data, size_t size, size_t k_size){

		size_t i = 0;
		//k_2_m = k_2_p if k_size is even
		const size_t k_2_m = k_size/2;
		const size_t k_2_p = k_size-k_2_m;

		T buffer[k_2_p];

		auto sum = [](T* begin, T* end)->T{
			T sum = 0;
			for (T* itr = begin; itr<end; itr++)
				sum += *itr;
			return sum;
		};

		for (NULL; i<k_2_p; i++){
			buffer[i] = 
				sum(data, data+i+k_2_p)/(i+k_2_p);
		}

		for (NULL; i+k_2_m<size; i++){
			T tmp = buffer[i%k_2_p];
			buffer[i%k_2_p] = 
				sum(data+i-k_2_m, 
					data+i+k_2_p)/(k_size);
			data[i-k_2_p] = tmp;
		}

		for (NULL; i<size; i++){
			T tmp = buffer[i%k_2_p];
			buffer[i%k_2_p] = 
				sum(data+i-k_2_m, 
					data+size)/(size-i+k_2_m);
			data[i-k_2_p] = tmp;
		}

		for (NULL; i<size+k_2_p; i++){
			data[i-k_2_p] = buffer[i%k_2_p];
		}
	}

	void gaussian_lpf(double *data, size_t size, double a){
		double *out = new double[size];
		for (int i=0; i<(int)size; i++){
			double tmp = 0;
			for (int j=0; j<(int)size; j++){
				tmp += data[j]*gaussian(j-i,a);
			}
			out[i] = tmp;
		}
		memcpy(data, out, size*sizeof(double));
		delete[] out;
		return;
	}
}
