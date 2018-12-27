#include <vector>
#include <Eigen/Dense>
#include <cmath>

const double PI	=	3.1415926535897932384;
Eigen::MatrixXd M_L(25,25);
Eigen::MatrixXd M_R(25,25);
Eigen::MatrixXd S_L(25,25);
Eigen::MatrixXd S_R(25,25);
Eigen::MatrixXd T1(25,25);
Eigen::MatrixXd T2(25,25);
Eigen::MatrixXd T3(25,25);
Eigen::MatrixXd T4(25,25);

class FMM1DLine {
public:
	int lineNumber;
	int parentNumber;
	int childrenNumbers[2];
	int neighborNumbers[2];
	int innerNumbers[2];
	int outerNumbers[2];

	FMM1DLine () {
		lineNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<1; ++l) {
			childrenNumbers[l]	=	-1;
                        neighborNumbers[l]	=	-1;
                        innerNumbers[l]		=	-1;
                        outerNumbers[l]		=	-1;
		}
	}
        
	double center;
	Eigen::VectorXd multipoles;
	Eigen::VectorXd locals;
	std::vector<int> local_x_index;
	//	following are defined only at the leaf
	//Eigen::VectorXd alpha;
	//Eigen::VectorXd location_of_charge_x;
};


class FMM1DTree {
public:
	int nLevels;			//	Number of levels in the tree.
        int nChebNodes;		//=p	//	Number of Chebyshev nodes along one direction.
	int N;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	int s;					//	no. of points per level at finest level
	Eigen::VectorXd y;				//	locations where f needs to be evaluated (size=N)
	Eigen::VectorXd x;				//	locations where xk 's are located (size=N)
	Eigen::VectorXd alpha;	
	std::vector<int> nLinesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> lineLength;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	
	std::vector<std::vector<FMM1DLine> > tree;	//	The tree storing all the information.

	//	Chebyshev nodes
	std::vector<double> standardChebNodes1D;

	//Eigen::MatrixXd cd(25,25);

	//	Different Operators
	Eigen::MatrixXd selfInteraction;		//	Needed only at the leaf level.
	Eigen::MatrixXd neighborInteraction[8];		//	Neighbor interaction only needed at the leaf level.
	//Eigen::MatrixXd M_L(25,25);					//	Transfer from multipoles of 2 children to multipoles of parent.
 	/*Eigen::MatrixXd M_R(25,25);		
	Eigen::MatrixXd S_L(25,25);		
	Eigen::MatrixXd S_R(25,25);		
	Eigen::MatrixXd T1(25,25);	
	Eigen::MatrixXd T2(25,25);	
	Eigen::MatrixXd T3(25,25);	
	Eigen::MatrixXd T4(25,25);				//	Transfer from locals of parent to locals of 2 children.

	/*M_L	=	Eigen::MatrixXd(25,25);					//	Transfer from multipoles of 2 children to multipoles of parent.
	M_R	=	Eigen::MatrixXd(25,25);
	S_L	=	Eigen::MatrixXd(25,25);
	S_R	=	Eigen::MatrixXd(25,25);
	T1	=	Eigen::MatrixXd(25,25);
	T2	=	Eigen::MatrixXd(25,25);
	T3	=	Eigen::MatrixXd(25,25);
	T4	=	Eigen::MatrixXd(25,25);*/

	Eigen::VectorXd y_in_Line;
	Eigen::VectorXd x_in_Line;		
	Eigen::VectorXd Far_Field;		//Far field at N yk points
	Eigen::VectorXd Near_Field;		//Net field at N yk points	
	Eigen::VectorXd Total_Field;		//Far field at N yk points
	Eigen::VectorXd Direct_Field1;
	Eigen::VectorXd Direct_Field2;
	Eigen::VectorXd Error1;
	Eigen::VectorXd Error2;
	
	FMM1DTree(int nLevels, int nChebNodes, int N, const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& alpha) {
		this->s				=	s;
		this->y				=	y;
		this->x				=	x;
		this->alpha			=	alpha;
                //std::cout << std::endl << "y: " << y << std::endl;
		this->nLevels			=	nLevels;
		//std::cout << std::endl << "this->nLevels " << this->nLevels << std::endl;
		this->nChebNodes		=	nChebNodes;//=p
		L				=	1;
	        nLinesPerLevel.push_back(1);//level 0
		lineLength.push_back(L);//semi length
		for (int k=1; k<=nLevels; ++k) {
			nLinesPerLevel.push_back(2*nLinesPerLevel[k-1]);
			lineLength.push_back(0.5*lineLength[k-1]);
		}
		//this->smallestBoxSize			=	lineLength[nLevels];
		//this->a					=	smallestBoxSize;
		this->N					=	N;
	}

        //	set_Standard_Cheb_Nodes
	void set_Standard_Cheb_Nodes() {
		for (int k=1; k<=nChebNodes; ++k) {
			standardChebNodes1D.push_back(-cos((k-0.5)/nChebNodes*PI));
		}
        }

	void createTree() {
		//	First create root and add to tree
		FMM1DLine root;
		root.lineNumber		=	0;
		root.parentNumber	=	-1;
		
		for (int l=0; l<2; ++l) {
			root.childrenNumbers[l]	=	l;
                        root.neighborNumbers[l]	=	-1;
           		root.innerNumbers[l]	=	-1;
			root.outerNumbers[l]	=	-1;
		}
		
		std::vector<FMM1DLine> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<FMM1DLine> level;
			for (int k=0; k<nLinesPerLevel[j]; ++k) {
				FMM1DLine line;
				line.lineNumber		=	k;
				line.parentNumber	=	k/2;
				for (int l=0; l<2; ++l) {
					line.childrenNumbers[l]	=	2*k+l;
				}
				level.push_back(line);
			}
			tree.push_back(level);
		}
	}

	void assign_Center_Location() {
		int J, K;
		tree[0][0].center	=	0.0;		
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*lineLength[j];
			for (int k=0; k<nLinesPerLevel[j]; ++k) {
				K	=	2*k;
				tree[J][K].center	=	tree[j][k].center-shift;
				tree[J][K+1].center	=	tree[j][k].center+shift;
				
			}
		}
	}


	double lag_poly(double eval_at_x,int j,int degree){//actual degree is degree-1//in this program degree is p
	    int k;	
	    
	    // method 1
		/////////////////////////////////////////////////////////////
	    /*int n_neg =1;
	    double u=0.0;
	    for (k=1;k<=degree;k++){
	        if((k-1)!=j){
		   if(eval_at_x-standardChebNodes1D[k-1] < 0)
			n_neg = n_neg*-1;
		   if(standardChebNodes1D[j]-standardChebNodes1D[k-1] < 0)
			n_neg = n_neg*-1;
	           u=u+log(fabs(eval_at_x-standardChebNodes1D[k-1]))-log(fabs(standardChebNodes1D[j]-standardChebNodes1D[k-1]));
		}
	    }
	    u=n_neg*exp(u);*/
		/////////////////////////////////////////////////////////////

            //method 2

	    double u=1.0;
	    
	    for (k=1;k<=degree;k++){
	        if((k-1)!=j)
	           u=u*(eval_at_x-standardChebNodes1D[k-1])/(standardChebNodes1D[j]-standardChebNodes1D[k-1]);
	    }
	/////////////////////////////////////////////////////////////

	    return (u);
	}


	void create_operators() {
		int i,j;
                //resize all the operators to size=nChebNodes
		M_L.resize(nChebNodes,nChebNodes);
		M_R.resize(nChebNodes,nChebNodes);
		S_L.resize(nChebNodes,nChebNodes);
		S_R.resize(nChebNodes,nChebNodes);
		T1.resize(nChebNodes,nChebNodes);
		T2.resize(nChebNodes,nChebNodes);
		T3.resize(nChebNodes,nChebNodes);
		T4.resize(nChebNodes,nChebNodes);

       		for(i=0; i<nChebNodes; i++){
		       for(j=0; j<nChebNodes; j++){
		       	   M_L(i, j)	=	lag_poly(3*standardChebNodes1D[i]/(6+standardChebNodes1D[i]),j,nChebNodes);
			   M_R(i, j)	=	lag_poly(3*standardChebNodes1D[i]/(6-standardChebNodes1D[i]),j,nChebNodes);

			   /*M_L(i, j)	=	lag_poly(standardChebNodes1D[i]/(6+standardChebNodes1D[i]),j,nChebNodes);
			   M_R(i, j)	=	lag_poly(standardChebNodes1D[i]/(6-standardChebNodes1D[i]),j,nChebNodes);*/

			   /*M_L(i, j)	=	lag_poly(standardChebNodes1D[i]/(2+standardChebNodes1D[i]),j,nChebNodes);
			   M_R(i, j)	=	lag_poly(standardChebNodes1D[i]/(2-standardChebNodes1D[i]),j,nChebNodes);*/

      			   S_L(i, j)	=	lag_poly((standardChebNodes1D[i]-1)/2,j,nChebNodes);
			   S_R(i, j)	=	lag_poly((standardChebNodes1D[i]+1)/2,j,nChebNodes);
			   T1(i, j)	=	lag_poly(3/(standardChebNodes1D[i]-6),j,nChebNodes);
			   T2(i, j)	=	lag_poly(3/(standardChebNodes1D[i]-4),j,nChebNodes);
			   T3(i, j)	=	lag_poly(3/(standardChebNodes1D[i]+4),j,nChebNodes);
			   T4(i, j)	=	lag_poly(3/(standardChebNodes1D[i]+6),j,nChebNodes);
		        }
		}  
 	 }



	//	Binary Search
	//	to find the line on which yk lies
	int binary_search(int l0, int ln, double yk) {
		int lm = (l0+ln)/2;
		//std::cout << std::endl << "lm" << lm << std::endl;
		if (l0 > ln) {
			return(-1);
			
		}
		else if (l0 == ln) {
			return(lm);
		}
		else
				if (yk <= tree[nLevels][lm].center + lineLength[nLevels]) {
					return(binary_search(0, lm, yk));
					
				}
				else
					return(binary_search(lm+1, ln, yk));
		}


//creating tree,operators
//numbering of children,neighbours,interaction list
//multipoles,local,near,near+far,error,time,main


	void assign_Leaf_Charges() {
		//#pragma omp parallel for
		int k,l;
		x_in_Line.resize(N);
		//std::cout << std::endl << "lineLength[nLevels]" << lineLength[nLevels] << std::endl;
		//x_in_Line(0) = binary_search(0, nLinesPerLevel[nLevels]-1, x(0));
		
		for (k=0; k<N; k++) {
			x_in_Line(k) = binary_search(0, nLinesPerLevel[nLevels]-1, x(k));
			//std::cout << std::endl << "x_in_Line( " << k <<  ")= " << x_in_Line(k) << std::endl;
			for (l=0; l < nLinesPerLevel[nLevels]; l++) {
				if(x_in_Line(k) == l)
					tree[nLevels][l].local_x_index.push_back(k);
			}
		}
		/*for (l=0; l < nLinesPerLevel[nLevels]; l++) {
			std::cout << std::endl << "line: " << l << std::endl;
				for (k=0; k < tree[nLevels][l].local_x_index.size(); k++) 
					std::cout << std::endl << tree[nLevels][l].local_x_index[k] << std::endl;
		}*/
		/*for (int k=0; k<nLinesPerLevel[nLevels]; ++k) {
			tree[nLevels][k].alpha		=	0.5*(Eigen::VectorXd::Ones(s) + Eigen::VectorXd::Random(s));
			tree[nLevels][k].location_of_charge_x		=	tree[nLevels][k].center*Eigen::VectorXd::Ones(s) + lineLength[nLevels]*Eigen::VectorXd::Random(s);
		}*/
	}


	//multipoles 
	void evaluate_multipoles() {
		int i,j,k,l,g;
		//	multipoles at finest level
		for (int k=0; k<nLinesPerLevel[nLevels]; ++k) {//lines
			tree[nLevels][k].multipoles	=	Eigen::VectorXd::Zero(nChebNodes);
			for(i=0; i<nChebNodes; i++){//p term multipole
				
				for(j=0; j<(tree[nLevels][k].local_x_index.size()); j++ ) {//xk 's in line l
					tree[nLevels][k].multipoles(i)	=	tree[nLevels][k].multipoles(i) + alpha(tree[nLevels][k].local_x_index[j])*standardChebNodes1D[i]/(3*lineLength[nLevels]-standardChebNodes1D[i]*(x(tree[nLevels][k].local_x_index[j])-tree[nLevels][k].center));
				}
			}
		}
		//	M2M = multipoles at coarser levels	
		for (l=nLevels-1; l>=1; l--) {
			for (int k=0; k<nLinesPerLevel[l]; ++k) {
				//std::cout << std::endl << "l: " << l << "   k: " << k << std::endl;
				tree[l][k].multipoles = M_L*tree[l+1][2*k].multipoles + M_R*tree[l+1][2*k+1].multipoles;
			}
		}
	}
			
		

	//locals 
	void evaluate_locals() {
		int k, l;
		tree[1][0].locals = Eigen::VectorXd::Zero(nChebNodes);
		tree[1][1].locals = Eigen::VectorXd::Zero(nChebNodes);
		for (l=1; l<nLevels; l++) { //	evaluating for all children of level l
			for (int k=0; k<nLinesPerLevel[l]; ++k) {
				tree[l+1][2*k].locals = S_L*tree[l][k].locals;
				tree[l+1][2*k+1].locals = S_R*tree[l][k].locals;
				
				/*if(2*k-2 >= 0 && 2*k-2 < nLinesPerLevel[l+1]){
					tree[l+1][2*k].locals = tree[l+1][2*k].locals + T1*tree[l+1][2*k-2].multipoles ;
					tree[l+1][2*k+1].locals = tree[l+1][2*k+1].locals + T1*tree[l+1][2*k-2].multipoles;
				}
				if(2*k-1 >= 0 && 2*k-1 < nLinesPerLevel[l+1]){
					tree[l+1][2*k+1].locals = tree[l+1][2*k+1].locals + T2*tree[l+1][2*k-1].multipoles;
				}
				if(2*k+2 >= 0 && 2*k+2 < nLinesPerLevel[l+1]){
					tree[l+1][2*k].locals = tree[l+1][2*k].locals + T3*tree[l+1][2*k+2].multipoles;
				}
				if(2*k+3 >= 0 && 2*k+3 < nLinesPerLevel[l+1]){
					tree[l+1][2*k].locals = tree[l+1][2*k].locals + T4*tree[l+1][2*k+3].multipoles ;
					tree[l+1][2*k+1].locals = tree[l+1][2*k+1].locals + T4*tree[l+1][2*k+3].multipoles;
				}*/

				if(2*k-2 >= 0 && 2*k-2 < nLinesPerLevel[l+1]){
					tree[l+1][2*k].locals = tree[l+1][2*k].locals + T3*tree[l+1][2*k-2].multipoles ;
					tree[l+1][2*k+1].locals = tree[l+1][2*k+1].locals + T4*tree[l+1][2*k-2].multipoles;
				}
				if(2*k-1 >= 0 && 2*k-1 < nLinesPerLevel[l+1]){
					tree[l+1][2*k+1].locals = tree[l+1][2*k+1].locals + T3*tree[l+1][2*k-1].multipoles;
				}
				if(2*k+2 >= 0 && 2*k+2 < nLinesPerLevel[l+1]){
					tree[l+1][2*k].locals = tree[l+1][2*k].locals + T2*tree[l+1][2*k+2].multipoles;
				}
				if(2*k+3 >= 0 && 2*k+3 < nLinesPerLevel[l+1]){
					tree[l+1][2*k].locals = tree[l+1][2*k].locals + T1*tree[l+1][2*k+3].multipoles ;
					tree[l+1][2*k+1].locals = tree[l+1][2*k+1].locals + T2*tree[l+1][2*k+3].multipoles;
				}
				//std::cout << std::endl <<  "   multipoles: " << tree[l+1][2*k-2].multipoles  << std::endl;
				//tree[l+1][2*k].locals = S_L*tree[l][k].locals + flag1*T3*tree[l+1][2*k-2].multipoles + flag3*T2*tree[l+1][2*k+2].multipoles + flag4*T3*tree[l+1][2*k+3].multipoles ;

				//tree[l+1][2*k+1].locals = S_R*tree[l][k].locals + flag1*T4*tree[l+1][2*k-2].multipoles + flag2*T3*tree[l+1][2*k-1].multipoles + flag4*T2*tree[l+1][2*k+3].multipoles ;
				}
		}
	}



	//	far Field Evaluation
	void evaluate_far_field() {
		int j,k;
		Far_Field = Eigen::VectorXd::Zero(N);
		y_in_Line.resize(N);
		//std::cout << std::endl << "y "  << y << std::endl;
		for (k=0; k<N; k++) {
			y_in_Line(k) = binary_search(0, nLinesPerLevel[nLevels]-1, y(k));
			//std::cout << std::endl << "y_in_Line(: " << k <<  "   " << y_in_Line(k) << std::endl;
			for (j=0; j<nChebNodes; j++) {
				Far_Field(k)	=	Far_Field(k) + tree[nLevels][y_in_Line(k)].locals(j)*lag_poly((y(k)-tree[nLevels][y_in_Line(k)].center)/lineLength[nLevels], j, nChebNodes);
			}
		}
	}	



	//	near Field Evaluation
	void evaluate_near_field() {
		int j,k;
		Near_Field = Eigen::VectorXd::Zero(N);
		for (k=0; k<N; k++) {
			if(y_in_Line(k) >= 1) {
				for(j=0; j<(tree[nLevels][y_in_Line(k)-1].local_x_index.size()); j++ ){ //s no. of charges in each interval at finest level
					Near_Field(k) = Near_Field(k) + alpha((tree[nLevels][y_in_Line(k)-1].local_x_index[j]))/(y(k)-x((tree[nLevels][y_in_Line(k)-1].local_x_index[j])));
				}
			}
			if(y_in_Line(k)+1 < nLinesPerLevel[nLevels]) {
				for(j=0; j<(tree[nLevels][y_in_Line(k)+1].local_x_index.size()); j++ ){ //s no. of charges in each interval at finest level
					Near_Field(k) = Near_Field(k) + alpha((tree[nLevels][y_in_Line(k)+1].local_x_index[j]))/(y(k)-x((tree[nLevels][y_in_Line(k)+1].local_x_index[j])));
				}
			}
			for(j=0; j<(tree[nLevels][y_in_Line(k)].local_x_index.size()); j++ ){ //s no. of charges in each interval at finest level
				Near_Field(k) = Near_Field(k) + alpha((tree[nLevels][y_in_Line(k)].local_x_index[j]))/(y(k)-x((tree[nLevels][y_in_Line(k)].local_x_index[j])));
			}
		}
	}


	//	near Field Evaluation
	void evaluate_total_field() {
		Total_Field = Near_Field + Far_Field;
	}


	
	//interpolation
	void poly_interpolation() {
		int k,j;
		Eigen::VectorXd r	=	Eigen::VectorXd::Zero(N);
		int nneg = 1;
		for(j=0; j<N; j++) {
			nneg = 1;
			for(k=0; k<N; k++) {
				if( y(j)-x(k) < 0) {
					nneg = nneg*-1;
				}
				r(j) = r(j) + log(fabs(y(j)-x(k)));
			}
			r(j)	=	nneg * exp(r(j));
			Total_Field(j)	=	Total_Field(j) * r(j);
		}
		
	}
	

	//	Direct Evaluation 1
	void direct_evaluation1() {
		int j,k,l;
		Eigen::VectorXd g = Eigen::VectorXd::Ones(N);
		Eigen::VectorXd p = Eigen::VectorXd::Ones(N);
		Eigen::VectorXd r = Eigen::VectorXd::Ones(N);
		Direct_Field1 = Eigen::VectorXd::Zero(N);
		for (j=0; j<N; j++) {
			for (k=0; k<N; k++) {
				if(k==j) continue;
				g(j) = g(j)/(x(j)-x(k));
			}
		}


		for (l=0; l<N; l++) {
			for (k=0; k<N; k++) {
				r(l) = r(l)*(y(l)-x(k));
				p(k) = exp(-4*x(k)*x(k))/(y(l)-x(k));
				Direct_Field1(l) = Direct_Field1(l) + p(k)*g(k);
			}
			Direct_Field1(l) = Direct_Field1(l)*r(l);
		}
	}



//	Direct Evaluation 2
	void direct_evaluation2() {
		int j,k,l;
		int nneg1 = 1;
		int nneg2 = 1;
		Eigen::VectorXd g = Eigen::VectorXd::Zero(N);
		Eigen::VectorXd p = Eigen::VectorXd::Zero(N);
		Eigen::VectorXd r = Eigen::VectorXd::Zero(N);
		Direct_Field2 = Eigen::VectorXd::Zero(N);
		for (j=0; j<N; j++) {
			nneg2	=	1;
			for (k=0; k<N; k++) {
				if(k==j) continue;
				if( x(j)-x(k) < 0) {
					nneg2 = nneg2 * -1;
				}
				g(j) = g(j)+log(fabs(x(j)-x(k)));
			}
			g(j)	=	nneg2 * exp(-g(j));
		}



		for (l=0; l<N; l++) {
			nneg1	=	1;
			for (k=0; k<N; k++) {
				if( y(l)-x(k) < 0) {
					nneg1 = nneg1*-1;
				}
				r(l) = r(l) + log(fabs(y(l)-x(k)));
				p(k) = exp(-4*x(k)*x(k))/(y(l)-x(k));
				Direct_Field2(l) = Direct_Field2(l) + p(k)*g(k);
			}
			Direct_Field2(l) = Direct_Field2(l)*exp(r(l))*nneg1;
		}
	}




	//	Error_infinity norm Evaluation
	void error_evaluation() {
		int k, j,i;
		//std::cout << std::endl << "Direct_Field.size()"  << Direct_Field.size() << std::endl;
		//std::cout << std::endl << "Total_Field.size()" << Total_Field.size() << std::endl;
		Error1		=	Direct_Field1 - Total_Field;
		Error2		=	Direct_Field2 - Total_Field;
		//std::cout << std::endl << "alpha" << alpha << std::endl;
		//std::cout << std::endl << "y" << y << std::endl;
		//std::cout << std::endl << "x" << x << std::endl;
		
		//std::cout << std::endl << "Far_Field" << Far_Field << std::endl;
		//std::cout << std::endl << "Near_Field" << Near_Field << std::endl;
		//std::cout << std::endl << "Algo_Field" << Total_Field << std::endl;

		//std::cout << std::endl << "Direct_Field1" << Direct_Field1 << std::endl;
		//for(i=0;i<500;i++)
		//	std::cout << std::endl << "Direct_Field2" << Direct_Field2 ;//<< std::endl;

                /*Eigen::MatrixXd comp_fields(N,10);
		comp_fields.col(0) = x;
		comp_fields.col(1) = y;
		comp_fields.col(2) = x-y;
		comp_fields.col(3) = Total_Field;
		comp_fields.col(4) = Direct_Field2;
		comp_fields.col(5) = Direct_Field1;
		comp_fields.col(6) = Far_Field;
		comp_fields.col(7) = Near_Field;
		comp_fields.col(8) = Error1;
		comp_fields.col(9) = Error2;
		std::cout << std::endl << "comp_fields" << comp_fields;*/

		//std::cout << std::endl << "Error1" << Error1 << std::endl;
		//std::cout << std::endl << "Error2" << Error2 << std::endl;
		//std::cout << std::endl << "Far_Field-Direct_Field" << Far_Field-Direct_Field << std::endl;
		
		//std::cout << std::endl << "standardChebNodes1D" <<  std::endl;
		//for (k=0;k<nChebNodes;k++)
		//	std::cout << std::endl << standardChebNodes1D[k] << std::endl;
		//double lag_pol_check	=	lag_poly(standardChebNodes1D[1],2,nChebNodes);
		//std::cout << std::endl << "lag_pol_check" << lag_pol_check << std::endl;
		//std::cout << std::endl << "M_L: " << std::endl;
		//std::cout << std::endl << M_L << std::endl;
		
		/*Eigen::MatrixXd all_multipoles_leaf=Eigen::MatrixXd::Ones(nChebNodes,nLinesPerLevel[nLevels]);
		
		for(k=0;k<nLinesPerLevel[nLevels];k++)
			all_multipoles_leaf.col(k) = tree[nLevels][k].multipoles;
		std::cout << std::endl << "all_multipoles_leaf: " << std::endl;
		std::cout << std::endl << all_multipoles_leaf << std::endl;*/

		/*Eigen::MatrixXd all_multipoles_leaf_1=Eigen::MatrixXd::Ones(nChebNodes,nLinesPerLevel[nLevels-2]);
		
		for(k=0;k<nLinesPerLevel[nLevels-2];k++)
			all_multipoles_leaf_1.col(k) = tree[nLevels-2][k].multipoles;
		std::cout << std::endl << "all_multipoles_leaf_1: " << std::endl;
		std::cout << std::endl << all_multipoles_leaf_1 << std::endl;*/
		
		/*Eigen::MatrixXd locals_leaf_1=Eigen::MatrixXd::Ones(nChebNodes,nLinesPerLevel[nLevels]);
		
		for(k=0;k<nLinesPerLevel[nLevels];k++)
			locals_leaf_1.col(k) = tree[nLevels][k].locals;
		std::cout << std::endl << "locals_leaf_1: " << std::endl;
		std::cout << std::endl << locals_leaf_1 << std::endl;*/
		//std::cout << std::endl << 0.5*(Eigen::VectorXd::Ones(5) + Eigen::VectorXd::Random(5)) << std::endl;
		
		/*Eigen::MatrixXd xy = Eigen::MatrixXd::Ones(N,2);
		xy.col(0) = x;
		xy.col(1) = y;
		std::cout << std::endl << "xy" << xy << std::endl;*/

		//std::cout << std::endl << "nChebNodes: " << nChebNodes << std::endl;
	}



	//	Error_infinity norm Evaluation
	/*double error_infinity_norm_evaluation() {
		double Error_inf_norm;
		Error_inf_norm	=	Error.lpNorm<Infinity>();
	}*/

		

//	Error_L2 norm Evaluation
	double error_L2_norm_evaluation1() {
		double Error_L2_norm1;
		Error_L2_norm1	=	Error1.norm();
		//Error_L2_norm1	=	Far_Field.norm();
		return (Error_L2_norm1);
	}


	double error_L2_norm_evaluation2() {
		double Error_L2_norm2;
		Error_L2_norm2	=	Error2.norm();
		return (Error_L2_norm2);
	}

	double error_infinity_norm_evaluation1() {
		Eigen::VectorXd error_li(N);
		for (int k=0; k<N; ++k) {
			error_li(k)	=	fabs(Error1(k));
		}
		//Error_L2_norm1	=	Far_Field.norm();
		return (error_li.maxCoeff());
	}

	double error_infinity_norm_evaluation2() {
		Eigen::VectorXd error_li(N);
		for (int k=0; k<N; ++k) {
			error_li(k)	=	fabs(Error2(k));
		}
		//Error_L2_norm1	=	Far_Field.norm();
		return (error_li.maxCoeff());
	}

};		
