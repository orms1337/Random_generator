#include "includes_lib.hpp"
#include "Random.hpp"


int main()
{
	//unsigned int start_time =  clock();
	//Generator_Type generator_seed( seed );
	Calc obj;

	ofstream file;
	file.open("data_knuth_b.txt");

	if (file.is_open()) {
		cout << "File is open\n";
		cout << "============\n";
	}

	#pragma omp parallel num_threads(count_threads)  
	{
		#pragma omp for private(obj)
		for (auto count = 0; count < max_distribution; count++ ) {
			switch ( count ) {
				case 0: {
					Distribution_Type_0 distr( 5.0, 2.0 );
					cout << "Normall\n";
					obj.generate_value(generator_seed, distr);
					obj.get_M_D_A_E(file, "Normall");
					obj.clear_all();
					break;
				}
				case 1: {
					Distribution_Type_1 distr( t, p );
					cout << "Binomial\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Binomial");
					obj.clear_all();
					break;
				}
				case 2: {
					Distribution_Type_2 distr( n );
					cout << "Chi squared distribution\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Chi squared distribution");
					obj.clear_all();
					break;
				}
				case 3: {
					Distribution_Type_3 distr( lambda );
					cout << "Exponential distribution\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Exponential distribution");
					obj.clear_all();
					break;
				}
				case 4: {
					Distribution_Type_4 distr( m, b );
					cout << "Fisher f distribution\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Fisher f distribution");
					obj.clear_all();
					break;
				}
				case 5: {
					Distribution_Type_5 distr( alpha, beta );
					cout << "Gamma distribution\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Gamma distribution");
					obj.clear_all();
					break;
				}
				case 6: {
					Distribution_Type_6 distr( g, h );
					cout << "Lognormal distribution\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Lognormal distribution");
					obj.clear_all();
					break;
				}
				case 7: {
					Distribution_Type_7 distr( my );
					cout << "Poisson distribution\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Poisson distribution");
					obj.clear_all();
					break;
				}
				case 8: {
					Distribution_Type_8 distr( s );
					cout << "Student t distribution\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Student t distribution");
					obj.clear_all();
					break;
				}
				case 9: {
					Distribution_Type_9 distr( k );
					cout << "Geometric distribution\n";
					obj.generate_value(generator_seed, distr);	
					obj.get_M_D_A_E(file, "Geometric distribution");
					obj.clear_all();
					break;
				}
			}
		}
	}
	int clocks = clock() / 1000;
	file << "Time (sec): " << setprecision(3) << clocks << "\n";
	file.close();
	//unsigned int end_time = clock();
	//unsigned int search_time = end_time - start_time;
	//cout << "Time: " << search_time / 1000;
	//cout << clock() / 1000 << "\n";
	//system("swapoff -a && sudo swapon -a");
}

