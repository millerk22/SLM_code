#include<algorithm>
#include<iostream>
#include<fstream>
#include<sstream>
#include<tuple>
#include<vector>

using namespace std;

struct sparse{
  vector<unsigned int> i;
  vector<tuple<unsigned int,double>> j;
};

vector<vector<double>> parse_mat(string name, int thread, unsigned int N)
{
  ifstream mystream(name+"_"+to_string(thread)+".dat");

  if(!mystream.is_open())
    cout << name << " did not open!" << endl;

  vector<vector<double>> A(N, vector<double>());
  string str;
  unsigned int i=0;
  while(getline(mystream, str))
  {
    stringstream ss(str);
    while(!ss.eof()){
      double temp;
      ss >> temp;
      A[i].push_back(temp);
    }
    A[i].pop_back(); // Somehow the last element is being read twice
    i++;
  }

  return A;
}
sparse parse_spr(string name, int thread)
{
  ifstream mystream(name+"_"+to_string(thread)+".dat");

  if(!mystream.is_open())
    cout << name << " did not open!" << endl;

  sparse A;
  while(!mystream.eof()){
    int tempi, tempj;
    double tempv;
    mystream >> tempi >> tempj >> tempv;
    A.i.push_back(tempi);
    A.j.push_back(make_tuple(tempj,tempv));
  }
  A.i.pop_back(); // Somehow the last element is being read twice
  A.j.pop_back();

  return A;
}

vector<int> parse_vec(string name, int thread)
{
  ifstream mystream(name+"_"+to_string(thread)+".dat");

  if(!mystream.is_open())
    cout << name << " did not open!" << endl;

  vector<int> myvec;
  while(!mystream.eof()){
    int temp;
    mystream >> temp;
    myvec.push_back(temp);
  }
  myvec.pop_back(); // Somehow the last element is being read twice

  return myvec;
}

double parse_dbl(string name, int thread)
{
  ifstream mystream(name+"_"+to_string(thread)+".dat");

  if(!mystream.is_open())
    cout << name << " did not open!" << endl;

  double mydbl;
  mystream >> mydbl;

  return mydbl;
}

int main(int argc, char * argv[]){

  //Take thread as user input
  int thread = 1;
  if(argc >1)
    thread = atoi(argv[1]);

  //Parse user inputs
  vector<int> g = parse_vec("g",thread);
  sparse      A = parse_spr("A",thread);
  vector<int> k = parse_vec("k",thread);
  double   twom = parse_dbl("twom",thread);
  double  gamma = parse_dbl("gamma",thread);
  vector<int> diagB = parse_vec("diagB",thread);

  //Print inputs
  cout << "g" << endl;
  for(unsigned int i=0; i<g.size(); i++){
    cout << g[i] << endl;
  }
  cout << "A" << endl;
  for(unsigned int i=0; i<A.i.size(); i++){
    //cout << A.i[i] << " " << A.j[i] << " " << A.v[i] << endl;
    cout << A.i[i] << " " << get<0>(A.j[i]) << " " << get<1>(A.j[i]) << endl;
  }
  cout << "k" << endl;
  for(unsigned int i=0; i<k.size(); i++){
    cout << k[i] << endl;
  }
  cout << "twom " << twom << endl;
  cout << "gamma " << gamma << endl;

  //Determine constants
  unsigned int N = g.size();
  unsigned int nhat = *max_element(begin(g), end(g));

//  vector<vector<unsigned int>> Nbrs(A.i.size(),vector<unsigned int>());
//  for(unsigned int x=0; x<A.i.size(); x++){
//    Nbrs[x].push_back(A.j[x]);
//  }
//
  vector<vector<double>> Cut = parse_mat("Cut", thread, nhat);
  vector<int> vol = parse_vec("vol", thread);

  cout << "Cut" << endl;
  for(unsigned int i=0; i<N; i++){
    for(unsigned int j=0; j<nhat; j++){
      cout << Cut[i][j] << " ";
    }
    cout << endl;
  }

  double tol = 10^(-4);
  double pass_improvement = std::numeric_limits<double>::max();
  while(2*pass_improvement/twom > tol){
    pass_improvement=0;
    for(unsigned int i=0; i<N; i++){
      unsigned int a_old=g[i];

      vector<double> I(nhat);
      for(unsigned int a=0; a<nhat; a++){
	I[a] = Cut[i][a] - gamma/twom * k[i]*vol[a];
      }

      unsigned int a = *max_element(begin(I), end(I));
      double maxvalue = I[a];

      if(a != a_old){
	g[i]=a;

	pass_improvement += maxvalue + diagB[i] - Cut[i][a_old] + gamma/twom * k[i] * vol[a_old];

	for(unsigned int j=0; j<A[i].size(); j++){
	  Cut[A[i][j],a_old] -= A[Nbrs[i][j],i];
	}

      }

    }
  }

}
