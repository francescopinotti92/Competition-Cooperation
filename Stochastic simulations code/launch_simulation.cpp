#include <iostream>                                              
#include <fstream>          
#include <string>           
#include <sstream>        
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <sys/stat.h>
#include <random>

using namespace std;
using std::vector;

//--DECLARE FUNCTIONS--//
vector <int> prevalence_compartment(const vector<int>& v, const vector<int>& S);
int which_community(int i, int N, int Nc);

int main(int argc, char *argv[])
{
  //----DECLARE VARIABLES----//
  
  //---EXTERNAL PARAMETER DESCRIPTION---//
  /*
    [1] (int)     network size
    [2] (int)     maximum time allowed for simulations
    [3] (double)  beta1
    [4] (double)  beta2
    [5] (double)  alpha
    [6] (double)  C1
    [7] (double)  C2
    [8] (double)  mu
    [9] (double)  average network degree
    [10] (int)    number of communities
    [11] (double) community structure strength
    [12] (int)    B1's starting community
    [13] (int)    B2's starting community
    [14] (int)    A's starting community 
  */
  
  //---PARAMETERS---//
  int N = stoi(argv[1]); //lattice side size
  int Tmax = stoi(argv[2]);

  double B1= stod(argv[3]), B2= stod(argv[4]), A= stod(argv[5]); //beta1, beta2 and alfa
  double C1= stod(argv[6]), C2= stod(argv[7]); //cooperation factors
  double mu= stod(argv[8]); //recovery probability
  double ave_degree = stod(argv[9]);
  int Nc = stoi(argv[10]);
  double Pin = stod(argv[11]);

  int where_B1_starts = stoi(argv[12]);
  int where_B2_starts = stoi(argv[13]);
  int where_A_starts = stoi(argv[14]);
  
  //---VARIABLES---//
  double rn;
  int v,w;
  int rn_int;
  int checkA;
  int encounterA;
  
  vector<int> ICvec;
  vector<int> S; //state vector: 0 = Susceptible, 1 = B1, 2 = B2, 3 = A, 4 = AB1, 5 = AB2 
  vector<int> Smemo; //buffer vector needed to update S at the end of each timestep
  vector<int> memo; //buffer vector needed to update S at the end of each timestep
  vector<int> infected; //contains the infected nodes
  vector<int> checked; //tells if a node has already been checked during the infection step
  vector<int> encounter;
  vector<int> infected_result;
  vector<vector<int> > contact_list;
  vector<double> bbeta;
  vector<double> Cs;
 
  // Network diagnostics
  int edges;

  // Community network creation
  int nstubs;
  vector<vector<int> > stubs_in;
  vector<vector<int> > stubs_out;
  vector<int> weights;
  vector<int> indexes;
  vector<int> check_size;
  vector<double> rweights;
  int first_choice, second_choice;

  bool connected;
  vector<vector<int> > patch_contact;
  vector<int> patch_status;

  // get seed from present time
  unsigned int random_seed = time(NULL);
  srand(random_seed + 30000);

  // define distributions used
  mt19937 mt(random_seed);
  uniform_real_distribution<double> rand1(0,1);
  uniform_int_distribution<int> randint1(0, N - 1);
  uniform_int_distribution<int> randint2(0, 1);
  uniform_int_distribution<int> randint_within(0,N/Nc - 1);
  uniform_int_distribution<int> randint_Nc(0,Nc - 1);
  poisson_distribution<int> randpoisson(ave_degree);
  exponential_distribution<double> randexp1(1);

  bbeta = {B1, B2};
  Cs = {C1, C2};
  ICvec = {100,100,100};

  //---INITIAL CONDITIONS---//
  
  memo.resize(0);
  infected.resize(0);
  checked.resize(0);
  S.resize(0);
  S.resize(N,0);
  
  // initial conditions
  
  //B1
  for(int i = 0; i < ICvec[0]; i++)
    {
      rn_int = where_B1_starts*(N/Nc) + randint_within(mt); 
      if (S.at(rn_int) == 0)
	{
	  S.at(rn_int) = 1;
	  infected.push_back(rn_int);
	}
      else i--;
    }
  //B2
  for(int i = 0; i < ICvec[1]; i++)
    {
      rn_int = where_B2_starts*(N/Nc) + randint_within(mt);
      if (S[rn_int] == 0)
	{
	  S[rn_int] = 2;
	  infected.push_back(rn_int);
	}
      else i--;
    }
      //A
  for(int i = 0; i < ICvec[2]; i++)
    {
      rn_int = where_A_starts*(N/Nc) + randint_within(mt);
      if (S[rn_int] == 0)
	{
	  S[rn_int] = 3;
	  infected.push_back(rn_int);
	}
      else i--;
    }
  
  //CREATE NETWORK
  
  connected = false;
  while (connected == false)
    {
      edges = 0;
      contact_list.resize(0);
      contact_list.resize(N, vector<int>());
      stubs_in.resize(0);
      stubs_in.resize(Nc, vector<int> {});
      stubs_out.resize(0);
      stubs_out.resize(Nc, vector<int> {});
      
      patch_contact.resize(0);
      patch_contact.resize(Nc, vector<int> {});
      
      //generate degrees
      for(int i = 0; i < N; i++)
	{
	  nstubs = randpoisson(mt);
	  for (int j = 0; j < nstubs; j++)
	    {
	      if (rand1(mt) < Pin) {stubs_in.at(which_community(i,N,Nc)).push_back(i);}
	      else {stubs_out[which_community(i,N,Nc)].push_back(i);}
	    }
	}
      
      // within community configuraion model
      random_shuffle (stubs_in.begin(), stubs_in.end());
      for (int kc = 0; kc < Nc; kc++)
	{
	  if (stubs_in[kc].size() % 2 != 0)
	    {
	      randint2.param(uniform_int_distribution<int>::param_type(0, stubs_in[kc].size()-1));
	      stubs_in[kc].push_back(stubs_in[kc][randint2(mt)]);
	    }
	  random_shuffle (stubs_in[kc].begin(), stubs_in[kc].end());
	  for (int k = 0; k < stubs_in[kc].size();)
	    {
	      if ((stubs_in[kc][k] != stubs_in[kc][k+1]) and (find(contact_list.at(stubs_in[kc][k]).begin(),contact_list.at(stubs_in[kc][k]).end(),stubs_in[kc][k+1])==contact_list.at(stubs_in[kc][k]).end()))
		{
		  contact_list[stubs_in[kc][k]].push_back(stubs_in[kc][k+1]);
		  contact_list[stubs_in[kc][k+1]].push_back(stubs_in[kc][k]);
		}
	      
	      k += 2;
	    }
	}
      
      //extra community configuration model                             
      
      weights.resize(0);
      rweights.resize(0);
      rweights.resize(Nc,0);
      indexes.resize(0);
      indexes.resize(Nc,0);
      check_size.resize(0);
      check_size.resize(Nc,0);
      
      //compute stubs numbers                                                                                                                                                                                       
      for (int kc = 0; kc < Nc; kc++)
	{
	  weights.push_back(stubs_out[kc].size());
	  if (stubs_out[kc].size() == 0) check_size[kc] = 1;
	  random_shuffle (stubs_out[kc].begin(), stubs_out[kc].end());
	}
      
      //pick edges until only one community has stubs                                                                                                                                                                
      while (accumulate(check_size.begin(), check_size.end(), 0) < Nc-1 )
	{
	  first_choice = -1;
	  second_choice = -1;
	  for (int k = 0; k < Nc; k++)
	    {
	      if (weights[k] > 0)
		{
		  rweights[k] = (double)randexp1(mt)/weights[k];
		  if (first_choice == -1 and check_size[k] == 0) first_choice = k;
		  if (rweights[k] < rweights[first_choice] and check_size[k] == 0) first_choice = k;
		}
	    }
	  for (int k = 0; k < Nc; k++)
	    {
	      if (second_choice == -1 and check_size[k] == 0 and first_choice != k) second_choice = k;
	      if (rweights[k] < rweights[second_choice] and check_size[k] == 0 and first_choice != k) second_choice = k;
	    }
	  weights[first_choice] -= 1;
	  if (weights[first_choice] == 0) check_size[first_choice] = 1;
	  weights[second_choice] -= 1;
	  if (weights[second_choice] == 0) check_size[second_choice] = 1;
	  v = stubs_out[first_choice].at(indexes[first_choice]);
	  w = stubs_out[second_choice].at(indexes[second_choice]);
	  indexes[first_choice] += 1;
	  indexes[second_choice] += 1;
	  
	  if (find(contact_list.at(v).begin(),contact_list.at(v).end(),w)==contact_list.at(v).end()) 
	    {
	      contact_list.at(v).push_back(w);                                                                                                                                                                      
	      contact_list.at(w).push_back(v);
	      if ((find(patch_contact.at(which_community(v,N,Nc)).begin(),patch_contact.at(which_community(v,N,Nc)).end(),which_community(w,N,Nc))==patch_contact.at(which_community(v,N,Nc)).end()))
		{
		  patch_contact.at(which_community(v,N,Nc)).push_back(which_community(w,N,Nc));
		  patch_contact.at(which_community(w,N,Nc)).push_back(which_community(v,N,Nc));
		}
	    }	  
	}
      
      //verify if the community network is connected with a snowball (SI) algorithm
      //pathological situations may still appear
      
      if (Nc > 1)
	{
	  patch_status.resize(0);
	  patch_status.resize(Nc,0);
	  if (patch_contact.at(0).size() > 0)
	    {
	      patch_status.at(0) = 1;
	      for (int ti = 0; ti < Nc+1; ti++)
		{
		  for (int ni = 0; ni < Nc; ni++) 
		    {
		      if (patch_status.at(ni) == 1) 
			{
			  for (auto itt = patch_contact.at(ni).begin(); itt != patch_contact.at(ni).end(); ++itt)
			    {
			      patch_status.at(*itt) = 1;
			    }
			}
		    }
		}
	    }
	  if (accumulate(patch_status.begin(), patch_status.end(), 0) == Nc) connected = true;
	}
      else
	{
	  connected = true;
	}
    }

  //---START SIMULATION---//
  for(int t = 0; t < Tmax; ++t)
    {
      memo.resize(0);
      checked.resize(0);
      checked.resize(N,0);
			      
      Smemo.resize(0);
      Smemo.resize(N,0);
      
      //---OUTPUT---//
      infected_result = prevalence_compartment(infected,S);
      
      //---Detailed output---/
      cout << t << " ";
      for (auto itt = infected_result.begin(); itt != infected_result.end(); ++itt) cout << " " << *itt;
      cout << endl;
      
      //---INFECTION---//
      //Loop through single infected
      for(vector<int>::iterator it = infected.begin(); it != infected.end(); ++it)
	{
	  for(vector<int>::iterator it1 = contact_list[*it].begin(); it1 != contact_list[*it].end(); ++it1)
	    {
	      if (S[*it1] < 4 and checked[*it1] == 0) //if *it1 is susceptible,B1,B2 or A, unchecked node is found try to infect it
		{
		  encounter.resize(0);
		  encounterA = 0;
		  for (auto it2 = contact_list[*it1].begin(); it2 != contact_list[*it1].end(); ++it2) 
		    {
		      if (S[*it2] >= 3) encounterA += 1;
		      if ((S[*it2] == 1) or (S[*it2] == 4)) encounter.push_back(1);
		      if ((S[*it2] == 2) or (S[*it2] == 5)) encounter.push_back(2);			      
		    }
		  if (S[*it1] == 0) //if *it1 is susceptible
		    {
		      checkA = 0;
		      rn = rand1(mt);
		      if(rn < 1-pow(1-A,encounterA)) checkA = 1;
		      random_shuffle (encounter.begin(), encounter.end());
		      for(auto enc = encounter.begin(); enc != encounter.end(); ++enc)
			{
			  if (rand1(mt) < bbeta[*enc - 1]) 
			    {
			      Smemo[*it1] += *enc;
			      break;
			    }
			}
					      
		      Smemo[*it1] += 3*checkA;
		      if (Smemo[*it1] > 0) memo.push_back(*it1);
		    }
		  else if (S[*it1] == 1) //if *it1 is in B1 compartment
		    {
		      if (rand1(mt) < 1-pow(1-Cs[0]*A,encounterA))
			{
			  Smemo[*it1] = 3;
			  memo.push_back(*it1);
			}
		    }
		  else if (S[*it1] == 2) //if *it1 is in B2 compartment
		    {
		      if (rand1(mt) < 1-pow(1-Cs[1]*A,encounterA))
			{
			  Smemo[*it1] = 3;
			  memo.push_back(*it1);
			}
		    }
		  else if (S[*it1] == 3) //if *it1 is in A compartment
		    {
		      random_shuffle (encounter.begin(), encounter.end());
		      for(auto enc = encounter.begin(); enc != encounter.end(); ++enc)
			{
			  if (rand1(mt) < Cs[*enc - 1] * bbeta[*enc - 1])
			    {
			      Smemo[*it1] = *enc;
			      memo.push_back(*it1);
			      break;
			    }
			}
		    }
		  checked[*it1] = 1;
		}
	    }
	}
      
      //---RECOVERY---//
      
      if (t >= 0) //MUST BE GREATER OR EQUAL THAN ZERO
	{
	  for(auto it = infected.begin(); it < infected.end();)
	    {
	      if (Smemo[*it] == 0)
		{
		  if (S[*it] < 4)
		    {
		      if (rand1(mt) < mu)
			{
			  S[*it] = 0;
			  swap(*it, infected.back()); 
			  infected.pop_back();
			}
		      else ++it;
		    }
		  else
		    {
		      rn = rand1(mt);
		      if (rn > (1-mu)*(1-mu))
			{
			  if(rn - (1-mu)*(1-mu) < mu*mu)
			    {
			      S[*it] = 0;
			      swap(*it, infected.back()); 
			      infected.pop_back();
			    }
			  else
			    {
			      if(rn - (1-mu)*(1-mu) - mu*mu < mu*(1-mu)) S.at(*it) -= 3;
			      else S[*it] = 3;
			      ++it;
			    }
			}
		      else ++it;
		    }
		}
	      else
		{
		  ++it;
		}
	    } 
	}
      
      //---UPDATE---//        
      for(vector<int>::iterator it = memo.begin(); it < memo.end(); ++it)
	{
	  if (S[*it] == 0)
	    {
	      S[*it] += Smemo[*it];
	      Smemo[*it] = 0;                                                                                              
	      infected.push_back(*it);
	    }
	  else
	    {
	      S[*it] += Smemo[*it];
	      Smemo[*it] = 0;                                                                                              
	    }
	}
      memo.resize(0);
    }
}

//count infected by compartment
vector <int> prevalence_compartment(const vector<int>& v, const vector<int>& S)
{
  vector <int> prevalence {0,0,0,0,0};
  for (auto it = v.begin(); it != v.end(); ++it)
    {
      //      cout << *it << " " << S.at(*it) << endl;
      prevalence.at(S[*it]-1) += 1;
    }
  return prevalence;
}

int which_community(int i, int N, int Nc) {return i/(N/Nc);}

