#include "includes_lib.hpp"
#include "Random.hpp"

Calc::Calc() {
	pick_out_mem();
}

Calc::~Calc() {
	free_mem();
}

void Calc::pick_out_mem() {
	mass_x = new double *[max_conf];
	for (auto i = 0; i < max_conf; i++) {
		mass_x[i] = new double [count_elem];
	}

	mass_M = new double [max_conf];
	mass_M_M = new double [max_conf];

	for (auto i = 0; i < max_conf; i++) {
		mass_M[i] = 0.0;
		mass_M_M[i] = 0.0;
	}
}

void Calc::free_mem() {
	for (int i = 0; i < max_conf; i++) {
		delete [] mass_x[i];
	}
	delete [] mass_x;
	delete [] mass_M;
	delete [] mass_M_M;

	cout << "Memory free\n";	
}

void Calc::init_mass( auto x, int k, int m ) {
	mass_x[k][m] = x;
}

void Calc::calc_M_D_A_E() {
	for (auto i = 0; i < max_conf; i++) {
		for (auto j = 0; j < count_elem; j++)
		{
			mass_M[i] += mass_x[i][j];
			mass_M_M[i] += mass_x[i][j] * mass_x[i][j];		
		}
		mass_M[i] /= count_elem;
		mass_M_M[i] /= count_elem;
	}

	for (auto i = 0; i < max_conf; i++) {
		m += mass_M[i] / max_conf;
		m_m += mass_M_M[i] / max_conf;
	} 

	d = m_m - m * m;
	g = sqrt(d);	

	for (auto i = 0; i < max_conf; i++) {
		for (auto j = 0; j < count_elem; j++) {
			a += mass_x[i][j] * mass_x[i][j] * mass_x[i][j] - 2 * mass_x[i][j] * mass_x[i][j] * m
			+ mass_x[i][j] * m * m - mass_x[i][j] * mass_x[i][j] * m
			+ 2 * mass_x[i][j] * m * m - m * m * m;
		}
		
	}
	a /= max_conf * count_elem;
	a /= g * g * g;

	for (auto i = 0; i < max_conf; i++) {
		for (auto j = 0; j < count_elem; j++) {
			e += ( mass_x[i][j] * mass_x[i][j] + m * m ) * ( mass_x[i][j] * mass_x[i][j] + m * m )
						- 4 * ( mass_x[i][j] * mass_x[i][j] *  mass_x[i][j] * m - mass_x[i][j] * mass_x[i][j] * m * m
						+ mass_x[i][j] * m * m * m );
		}
	}
	e /= max_conf * count_elem;
	e /= g * g * g * g;
	e -= 3;
}

void Calc::get_M_D_A_E(ofstream &f, string str) {
	f << str << "\n"; 
	f << "M[x]: " << m << "\n" << "D[x]: " << d << "\n"
	<< "As[x]: " << a << "\n" << "Ex[x]: " << e << "\n\n";
}

void Calc::clear_all() {
	m = 0.0; 
	m_m = 0.0; 
	d = 0.0; 
	g = 0.0; 
	e = 0.0; 
	a = 0.0;

	for (auto i = 0; i < max_conf; i++) {
		mass_M[i] = 0.0;
		mass_M_M[i] = 0.0;
		for (auto j = 0; j < count_elem; j++) {
			mass_x[i][j] = 0.0;
		}
	}
}

void Calc::generate_value(Generator_Type &generator_seed, Distribution_Type_0 distr) {
	
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generator_seed, Distribution_Type_1 distr) {
	
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generator_seed, Distribution_Type_2 distr) {
	
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generator_seed, Distribution_Type_3 distr) {
	
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generator_seed, Distribution_Type_4 distr) {
	
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generator_seed, Distribution_Type_5 distr) {
	
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generate_seed, Distribution_Type_6 distr) {
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generate_seed, Distribution_Type_7 distr) {
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generate_seed, Distribution_Type_8 distr) {
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}

void Calc::generate_value(Generator_Type &generate_seed, Distribution_Type_9 distr) {
	for (auto i = 0; i < max_conf; i++) {
		uniform_int_distribution<uint64_t> distr_seed(0, numeric_limits<uint64_t>::max());
		Generator_Type generator( distr_seed( generator_seed ) );
		for (int j = 0; j < count_elem; j++) {
			auto value = distr( generator );
			init_mass(value, i, j); 
		}
	}
	calc_M_D_A_E();
}