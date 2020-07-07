#ifndef Random_H
#define Random_H

using namespace std;

using Generator_Type = knuth_b;
//using Generator_Type = mt19937;
//using Generator_Type = default_random_engine;
//using Generator_Type = minstd_rand;

//using Distribution_Type_0 = bernoulli_distribution;
using Distribution_Type_0 = normal_distribution<double>;
using Distribution_Type_1 = binomial_distribution<int>;
using Distribution_Type_2 = chi_squared_distribution<double>;
using Distribution_Type_3 = exponential_distribution<double>;
using Distribution_Type_4 = fisher_f_distribution<double>;
using Distribution_Type_5 = gamma_distribution<double>;
using Distribution_Type_6 = lognormal_distribution<double>;
using Distribution_Type_7 = poisson_distribution<int>;
using Distribution_Type_8 = student_t_distribution<double>;
using Distribution_Type_9 = geometric_distribution<int>;

constexpr auto max_conf = 100;
constexpr auto count_elem = 1000000;

constexpr double p = 0.5;
constexpr int t = 5;
constexpr double n = 3.0;
constexpr auto seed = 1;
constexpr double lambda = 3.5;
constexpr double m = 10.0, b = 9.0;
constexpr double alpha = 2.0, beta = 2.0;
constexpr double g = 0.0, h = 1.0;
constexpr double my = 4.1;
constexpr double s = 10.0;
constexpr double k = 0.3; 

constexpr uint64_t count_threads = 2;
constexpr uint64_t max_distribution = 10;

Generator_Type generator_seed( seed );

class Calc {

public:
	Calc();
	~Calc();
	void get_M_D_A_E(ofstream &f, string str); //Mean, Dispersion, Asymmetry, Excess
	void clear_all();
	void generate_value(Generator_Type &generate_seed, Distribution_Type_0 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_1 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_2 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_3 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_4 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_5 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_6 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_7 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_8 distr);
	void generate_value(Generator_Type &generate_seed, Distribution_Type_9 distr);
private:
	void pick_out_mem();
	void free_mem();
	void init_mass( auto x, int k, int m );
	void calc_M_D_A_E();
	double **mass_x;
	double *mass_M;
	double *mass_M_M;
	double m = 0.0, m_m = 0.0, d = 0.0, g = 0.0, e = 0.0, a = 0.0;
};

#include "Random_impl.cpp"

#endif