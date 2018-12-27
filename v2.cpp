#include <iostream>
#include <omp.h>
#include "v2.hpp"
#include <Eigen/Eigenvalues>
using namespace std;
using namespace Eigen;
#include <fstream>


/*void Jacobi_roots(int n);

void  Jacobi_roots(int n) {
	Eigen::MatrixXd J=Eigen::MatrixXd::Zero(n,n);
	int i;
	for(i=1;i<=n-1;++i) {
		J(i-1,i)=i/sqrt(4*i^2-1);
		J(i,i-1)=i/sqrt(4*i^2-1);
	}
	
EigenSolver<MatrixXd> es;

es.compute(J, false);
cout << "The eigenvalues of A are: " << es.eigenvalues().col(0) << endl;


        /*SelfAdjointEigenSolver<MatrixXd> eigensolver(J);
	if (eigensolver.info() != Success) abort();
	Eigen::VectorXd Leg_roots=eigensolver.eigenvalues();
       	cout << "The eigenvalues of A are:\n" << Leg_roots << endl;*/
//}

int main(int argc, char* argv[]) {
	int nLevels		=	atoi(argv[1]);
	int N		=	atoi(argv[2]);
	double epsilon = 10^-17;
	double start, end;
	double c	=	5.0;
    	//double p		=	(log(1/epsilon))/log(c);
	int p		=	25;

	//Jacobi_roots(3);

	//Eigen::VectorXd y	=	Eigen::VectorXd::Zero(N);
	Eigen::VectorXd x	=	Eigen::VectorXd::Zero(N);
	//x.resize(N);
	//x	=	Eigen::VectorXd::Zero(N);
	/*if (N == 64)
		ifstream myfile( "Leg_Nodes64.txt" );
	else if (N == 128)
		ifstream myfile( "Leg_Nodes128.txt" );
	else if (N == 256)
		ifstream myfile( "Leg_Nodes256.txt" );
	else if (N == 512)
		ifstream myfile( "Leg_Nodes512.txt" );
	else if (N == 1024)
		ifstream myfile( "Leg_Nodes1024.txt" );
	else if (N == 2048)
		ifstream myfile( "Leg_Nodes2048.txt" );
	else if (N == 4096)
		ifstream myfile( "Leg_Nodes4096.txt" );*/
	
	ifstream myfile( "Leg_Nodes512.txt" );
        int i=0;
	double a;
        while (myfile >> a)
        {
        	x(i)	=	a;
		i++;
		//std::cout << std::endl << a << std::endl;
        }
	////////////////////////////////////////////////////////////////	
	Eigen::VectorXd y(N);
	/*Eigen::VectorXd x	=	(-1-1.0/N)*Eigen::VectorXd::Ones(N);
	for (int k=1; k<=N; k++){
		x(k-1)	=	x(k-1) + 2.0*k/N;
	}*/
	for (int k=1; k<=N; ++k) {
		y(k-1)=(-cos((k-0.5)/N*PI));
		//std::cout << std::endl << y(k-1) << std::endl;
	}
	//y = x + 0.2/N*Eigen::VectorXd::Random(N);
	////////////////////////////////////////////////////////////////

	

	Eigen::VectorXd alpha(N);
	
	for (int k=0; k<N; ++k) {
		alpha(k)	=	exp(-4*x(k)*x(k));
	}
	//alpha		=	0.5*(Eigen::VectorXd::Ones(N) + Eigen::VectorXd::Random(N));
	
	int k,j;
	Eigen::VectorXd s	=	Eigen::VectorXd::Zero(N);
	Eigen::VectorXd g(N);

	start	=	omp_get_wtime();
	int nneg = 1;
	for(j=0;j<N;j++) {
		nneg	=	1;
		for(k=0;k<N;k++) {
			if(k == j)
				continue;
			if( x(j)-x(k) < 0) {
 				nneg = nneg * -1;
			} 
			s(j) = s(j) + log(fabs(x(j)-x(k)));
		}
		s(j)	=	nneg * exp(-s(j));
		g(j)	=	alpha(j) * s(j);
	}
	//std::cout << std::endl << "x: " << x << std::endl;
	
	//std::cout << std::endl << "N: " << N << std::endl;
	//std::cout << std::endl << "&y: " << &y << std::endl;
	//std::cout << std::endl << "y[0]: " << y[0] << std::endl;
	/*std::cout << std::endl << "y[0]: " << y[0] << "&y: " << &y[0] << "y[0] using pointer: " << *(&y[0]) << std::endl;
	std::cout << std::endl << "y[1]: " << y[1] << "&y: " << &y[1] << "y[1] using pointer: " << *(&y[1]) << std::endl;
	std::cout << std::endl << "y[2]: " << y[2] << "&y: " << &y[2] << "y[2] using pointer: " << *(&y[2]) << std::endl;*/
	//std::cout << std::endl << "y[1] using pointer: " << *(&y[1]) << std::endl;
	//std::cout << std::endl << "y_orig: " << y << std::endl;
	FMM1DTree* A	=	new FMM1DTree(nLevels, p, N, x, y, g);
 	
	

	A->set_Standard_Cheb_Nodes();
	A->createTree();

	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;

	start	=	omp_get_wtime();

	A->assign_Center_Location();
	A->assign_Leaf_Charges();

	end		=	omp_get_wtime();
	double timeAssignCharges=	(end-start);
	std::cout << std::endl << "Time taken to assemble the charges is: " << timeAssignCharges << std::endl;

	start	=	omp_get_wtime();

	A->create_operators(); 

	end		=	omp_get_wtime();
	double timeCreateOperators =	(end-start);
	std::cout << std::endl << "Time taken to create operators is: " << timeCreateOperators << std::endl;
		
	start	=	omp_get_wtime();

	A->evaluate_multipoles(); 

	end		=	omp_get_wtime();
	double timeEvaluateMultipoles =	(end-start);
	std::cout << std::endl << "Time taken to evaluate multipoles is: " << timeEvaluateMultipoles << std::endl;

	start	=	omp_get_wtime();

	A->evaluate_locals();
	
	end		=	omp_get_wtime();
	double timeEvaluateLocals =	(end-start);
	std::cout << std::endl << "Time taken to evaluate locals is: " << timeEvaluateLocals << std::endl;

	start	=	omp_get_wtime();

	A->evaluate_far_field();

	end		=	omp_get_wtime();
	double timeEvaluate_far_field =	(end-start);
	std::cout << std::endl << "Time taken to evaluate Evaluate far field is: " << timeEvaluate_far_field << std::endl;

	start	=	omp_get_wtime();

	A->evaluate_near_field();

	end		=	omp_get_wtime();
	double timeEvaluate_near_field =	(end-start);
	std::cout << std::endl << "Time taken to evaluate Evaluate near field is: " << timeEvaluate_near_field << std::endl;

	start	=	omp_get_wtime();

	A->evaluate_total_field(); 
	A->poly_interpolation();

	end		=	omp_get_wtime();
	double timeEvaluate_total_field =	(end-start);
	std::cout << std::endl << "Time taken to evaluate Evaluate near field is: " << timeEvaluate_total_field << std::endl;

	double totalTime	=	timeCreateTree + timeAssignCharges + timeCreateOperators + timeEvaluateMultipoles + timeEvaluateLocals + timeEvaluate_far_field + timeEvaluate_near_field + timeEvaluate_total_field;
	std::cout << std::endl << "Total Time taken is: " << totalTime << std::endl;
	
	double totalTimeInit 	=	timeCreateTree + timeAssignCharges + timeCreateOperators;
	double totalTimeEval 	=	timeEvaluateMultipoles + timeEvaluateLocals + timeEvaluate_far_field + timeEvaluate_near_field;
	std::cout << std::endl << "Total Time taken for algorithm Initialisation is: " << totalTimeInit << std::endl;
	std::cout << std::endl << "Total Time taken for algorithm Evaluation is: " << totalTimeEval << std::endl;

	start	=	omp_get_wtime();

	A->direct_evaluation1();
	end		=	omp_get_wtime();
	double timeDirect_Evaluation1 =	(end-start);
	std::cout << std::endl << "Time taken for Direct_Evaluation 1 is: " << timeDirect_Evaluation1 << std::endl;


	start	=	omp_get_wtime();

	A->direct_evaluation2();
	end		=	omp_get_wtime();
	double timeDirect_Evaluation2 =	(end-start);
	std::cout << std::endl << "Time taken for Direct_Evaluation 2 is: " << timeDirect_Evaluation2 << std::endl;


	// error evaluation
	A->error_evaluation();
	
	//double Error_inf_norm	=	A->error_infinity_norm_evaluation();
	double Error_L2_norm1	=	A->error_L2_norm_evaluation1();
	//std::cout << std::endl << "E_inf: " << Error_inf_norm << std::endl;
	

	double Error_L2_norm2	=	A->error_L2_norm_evaluation2();
	//std::cout << std::endl << "E_inf: " << Error_inf_norm << std::endl;
	std::cout << std::endl << "error_L2_norm_evaluation2: " << A->error_L2_norm_evaluation2() << std::endl;
	std::cout << std::endl << "error_infinity_norm_evaluation2: " << A->error_infinity_norm_evaluation2() << std::endl;
}


