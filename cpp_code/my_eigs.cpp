//This is all the c++ you need to get spectral data for the operator M used in MBO for modularity.
// Only call this function from the directory in which it resides, or you can get segmentation faults. It was taylored specifically to interact with the matlab scripts in this folder and is not safe under other conditions.
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>


#include "./from_chris/SCC_vector.h"
#include "./from_chris/SparseOp/SparseOp.h"
#include "./from_chris/RandOpNd/RandVectorOp.h"
#include "./from_chris/RayleighChebyshev/RayleighChebyshev.h"

typedef SparseOp<SCC::vector> matrix;
typedef SCC::vector Vector;
typedef const int length;

using namespace std;

class M_op //Parent class of M and L ops
{
	public:
		M_op()
		{
			N=0;
		}
		M_op(const M_op& M_)//It is really important that the case where M_ is null be handled correctly, or else segmentation faults occur.
		{
			if(M_.N>0)
			{
				W = M_.W;
				k_ori = M_.k_ori;
				k = M_.k;
				gamma = M_.gamma;
				twom_ori = M_.twom_ori;
				N = M_.N;
				eig_type = M_.eig_type;
			}
			else
			{
				N=0;
			}
		}
		M_op(matrix W_, Vector k_ori_, double twom_ori_, double gamma_, string eig_type_)
		{
			W=W_;
			k_ori = k_ori_;
			gamma = gamma_;
			twom_ori = twom_ori_;

			N = W.getRowDimension();
			Vector k(N);
			for(int i=0; i<N; i++)
				k(i)=1;
			W.apply(k);
			this->k = k;
			eig_type = eig_type_;
			if(eig_type.compare("M") == 0)
			{
			}
			else if(eig_type.compare("L")==0)
			{
			}
			else
			{
				cout << "Bad type in M_op constructor" << endl;
			}
		}
		void apply(Vector& u)
		{
			Vector u_(u);
			W.apply(u_);
			//u = k*u;
			if(eig_type.compare("M") == 0)
			{
				u = k*u - u_ + 2*gamma/twom_ori*k_ori.dot(u)*k_ori;
			}
			else if(eig_type.compare("L")==0)
			{
				u=k*u - u_;
			}
			else
			{
				cout << "Bad type in M_op apply" << endl;
			}
		}
		long getRowDimension()
		{
			return N;
		}
	protected:
		matrix W;
		Vector k_ori;
		Vector k;
		double gamma;
		double twom_ori;
		long N;
		string eig_type;
};

//This encapsulates the interface with the RC code.
void get_eigs(M_op& M, long eigCount, vector<Vector>& V, vector<double>& D, int seed)
{
	RayleighChebyshev<Vector, M_op, SCC::RandVectorOp<Vector> > RC_op;
	//RC_op.setVerboseFlag();
	//RC_op.setEigDiagnosticsFlag();
	//RC_op.setVerboseSubspaceFlag();
	double iterationTol = 0.01;
	double minEigValue = 0;
	double maxEigValue = 0;
	SCC::RandVectorOp<Vector> randOp;
	//randOp.resetWithRandomSeed();
	randOp.resetSeed(seed);
	long N = M.getRowDimension();
	Vector u(N);
	randOp.randomize(u);
	RC_op.getInitialRCspectralEstimates(iterationTol, u, M, randOp, minEigValue, maxEigValue);

	RC_op.setMaxInnerLoopCount(100);

	double minEigValueEst = minEigValue;
	double maxEigValueBound = maxEigValue;
	double subspaceTol = 0.1;
	long subspaceIncrementSize = 8;
	long bufferSize = 8;
	randOp.randomize(u);
	Vector vStart = u;
	RC_op.getMinEigenSystem(eigCount, minEigValueEst, maxEigValueBound, subspaceTol, subspaceIncrementSize, bufferSize, vStart, M, randOp, D, V);
}

//Parser
M_op parse_W(int& eigCount, int thread, int& seed)
{
	std::ostringstream stringStream;
	stringStream << "params_";
	stringStream << thread;
	stringStream << ".dat";
	string copyOfStr = stringStream.str();
	ifstream params(copyOfStr);

	if(!params.is_open())
		cout << "Params did not open!" << endl;

	string line; 
	getline(params,line);
	istringstream iss( line );
	string result;
	getline(iss, result, ' ');
	long n_nnz = atoi( result.c_str() ); 
	getline(iss, result, ' ');
	long N = atoi( result.c_str() ); 
	getline(iss, result, ' ');
	double gamma = atof( result.c_str() ); 
	getline(iss, result, ' ');
	double twom = atof( result.c_str() ); 
	getline(iss, result, ' ');
	eigCount = atoi( result.c_str() ); 
	getline(iss, result, ' ');
	string eig_type = result;

	if(n_nnz == 0)//This shouldn't happen
	{
		cout << "Warning: zero matrix passed to my_eigs. This should not happen." << endl;
		eigCount = -1*eigCount;
	}

	std::ostringstream stringStream1;
	stringStream1 << "W_";
	stringStream1 << thread;
	stringStream1 << ".dat";
	copyOfStr = stringStream1.str();
	ifstream file(copyOfStr);

	if(!file.is_open())
		cout << "W File did not open!" << endl;

	matrix W(N,N);

	long i=0;long j=0;double val=0;
	for(int row=0; row<n_nnz; row++)
	{
		getline(file,line);
		istringstream iss( line );

		getline(iss, result, ' ');
		i = atoi( result.c_str() ); 
		getline(iss, result, ' ');
		j = atoi( result.c_str() ); 
		getline(iss, result, ' ');
		val = atof( result.c_str() ); 
		W.setOperatorData(i-1,j-1,val);
	}
	W.compact();
	W.sortColumnIndices();

	std::ostringstream stringStream2;
	stringStream2 << "k_ori_";
	stringStream2 << thread;
	stringStream2 << ".dat";
	copyOfStr = stringStream2.str();
	ifstream k_file(copyOfStr);

	if(!k_file.is_open())
		cout << "k_file did not open!" << endl;

	Vector k_ori(N);
	for(int i=0; i<N; i++)
	{
		getline(k_file,line);
		istringstream iss( line );

		getline(iss, result, ' ');
		k_ori(i) = atof( result.c_str() ); 
	}

	std::ostringstream stringStream3;
	stringStream3 << "seed_";
	stringStream3 << thread;
	stringStream3 << ".dat";
	copyOfStr = stringStream3.str();
	ifstream seed_file(copyOfStr);

	if(!seed_file.is_open())
		cout << "seed_file did not open!" << endl;

	getline(seed_file,line);
	istringstream isss( line );
	getline(isss, result, ' ');
	seed = atoi( result.c_str() ); 

	M_op M(W,k_ori,twom,gamma, eig_type);

	return M;
}

void write_file(vector<Vector> eigVectors, vector<double> eigValues, int thread)
{
	std::ostringstream stringStream2;
	stringStream2 << "V_";
	stringStream2 << thread;
	stringStream2 << ".txt";
	string copyOfStr = stringStream2.str();
	ofstream arrayData(copyOfStr);
	arrayData << scientific;

	const int N = eigVectors[0].getSize();
	const int eigCount = eigValues.size();

	for(int j=0;j<eigCount;j++){
		for(int i=0;i<N;i++){
			arrayData<<eigVectors[j](i)<<endl;
		}}

	std::ostringstream stringStream1;
	stringStream1 << "D_";
	stringStream1 << thread;
	stringStream1 << ".txt";
	copyOfStr = stringStream1.str();
	ofstream D_file(copyOfStr);
	D_file << scientific;

	for(int i=0; i<eigCount; i++)
	{
		D_file<<eigValues[i]<<endl;
	}
}

void test(vector<Vector> eigVectors, vector<double> eigValues, M_op M)
{
	length eigCount = eigValues.size();
	cout << "Printing percent errors" << endl;
	for(int i=0; i<eigCount; i++)
	{
		Vector temp = eigVectors[i] ;
		M.apply(temp);
		cout << (eigValues[i]*eigVectors[i] - temp).norm2() / (temp).norm2() << endl;
	}
}

int main(int argc, char * argv[])
{
	try
	{
		int thread = 1;
		if(argc >1)
		{
			thread = atoi(argv[1]); //This is a suffix to be used when reading a writing files. It makes it possible to run more than one instance of this program at a time.
		}

		int eigCount;//This is assigned in the next line
		int seed;
		M_op M = parse_W(eigCount,thread,seed);

		length N = M.getRowDimension();

		Vector u(N,1);

		vector<Vector> eigVectors;
		vector<double> eigValues;
		if(eigCount<0)//a zero matrix was passed
		{
			for(int i=0; i<(-1*eigCount); i++)
			{
				eigValues.push_back(0);
				eigVectors.push_back(Vector(N,0));
			}
		}
		else
		{
			get_eigs(M, eigCount, eigVectors, eigValues, seed);
		}

		write_file(eigVectors, eigValues, thread);

		//test(eigVectors,eigValues,M);
	}
	catch(const std::bad_alloc& e)
	{
		cout << e.what() << endl;
	}
	catch(...)
	{
		cout << "unknown exception" << endl;
	}
}
