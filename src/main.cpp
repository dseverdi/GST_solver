#include "GST.h"
#include <iostream>
#include <string>

using namespace GraphX;



enum sol_type_t { OPT, APPROX, HEURISTIC };
enum model_type_t { NR, GR, ALL };


int main(int argc, char** argv)
{

	std::cout << "-----------------------------------------------------------------------------" << std::endl;

	std::cout << "\n\n\r         Approximating Group Steiner Tree Problem                      " << std::endl;
	std::cout << "\n\r\t version: 0.9														   " << std::endl;
	std::cout << "\n\r\t author: S. Jelic, D. Severdija 									   " << std::endl;
	std::cout << "\n\n\r-----------------------------------------------------------------------------" << std::endl;



	std::string stp_file;
	sol_type_t solver_type;
	model_type_t model_type = ALL;
	short num_args = 0;


	if (argc>1)
	{
		for (int i = 0; i<argc; i++)
		{

			if (strcmp(argv[i], "-OPT") == 0) { solver_type = OPT; num_args++; }
			else if (strcmp(argv[i], "-APPROX") == 0) { solver_type = APPROX; num_args++; }
			else if (strcmp(argv[i], "-HEURISTIC") == 0) { solver_type = HEURISTIC; num_args++; }
			if (strcmp(argv[i], "-NR") == 0) { model_type = NR; num_args++; std::cout << "model: NR" << std::endl; }
			else if (strcmp(argv[i], "-GR") == 0) { model_type = GR; num_args++; std::cout << "model: GR" << std::endl; }
			if (strcmp(argv[i], "-f") == 0) { stp_file = std::string(argv[i + 1]); num_args++; }

		}
		if (num_args <2) { std::cout << "missing argument: \n  usage: GST_solver -OPT|APPROX|HEURISTIC -f [stp_input]" << std::endl; exit(1); }
	}


	else
	{
		// test without arguments from console
		solver_type = HEURISTIC;
		stp_file = "instances\\01_rand_n=60_m=462_k=11.stp";
		// stp_file = "instances\\rand_12_18_5.stp";
	}


	GST_instance test1(stp_file);
	GST_instance test2(stp_file);


	GST_solver sol_NR(test1);
	GST_solver sol_GR(test2);

	switch (solver_type)
	{
	case OPT:
		if (model_type == NR || model_type == ALL) sol_NR.solve_opt(type1);
		if (model_type == GR || model_type == ALL) sol_GR.solve_opt(type2);
		break;

	case APPROX:
		if (model_type == NR || model_type == ALL) sol_NR.solve_approx(type1);
		if (model_type == GR || model_type == ALL) sol_GR.solve_approx(type2);
		break;

	case HEURISTIC:
		GST_solver sol_H(test1);
		sol_H.solve_heuristic("ANTCOLONY");
		break;

	}

#if defined( _MSC_VER )
	system("PAUSE");
#endif 

	return EXIT_SUCCESS;
}


