
/*********************************************************************************
This program uses the C API for Gurobi to solve the zero-forcing problem using
the integer program that comes from the infection perspective of zero forcing

Author: Caleb Fast
Date: 2017

Compilation: g++ -std=c++11 -o find_forcing_set_infectionIP find_forcing_set_infectionIP.cpp -L<path to gurobi libraries> -l<gurobi###> -lgurobi_c++
Example: g++ -std=c++11 -o find_forcing_set_infectionIP find_forcing_set_infectionIP.cpp -L/opt/gurobi651/linux64/lib -lgurobi65 -lgurobi_c++

Usage: Change max_timesteps on line 105 to be a valid bound on maximum number of timesteps.
To run: ./find_forcing_set_infectionIP <path to .edg file>

**********************************************************************************/


#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<set>
#include<map>
#include<queue>
#include<cfloat>
#include<string>
#include<cstring>
#include<ctime>
#include<cmath>
#include<climits>
#include<stdlib.h>
#include<algorithm>
#include<chrono>
#include "../forcing_lib/forcing_lib.hpp"


/**INCLUDE PATH TO gurobi_c++.h**/
#include "/opt1/opt.LOCAL/gurobi752/linux64/include/gurobi_c++.h"
//Examples
//For my Laptop
//#include "/opt/gurobi651/linux64/include/gurobi_c++.h"

#define FLOAT_TOL 0.00001
#define STOPPING_TIME 7200
#define PROBLEM_TYPE 1 //1 for IP, 0 for LP

using namespace std;
using namespace forcing;



/*callback function for seperating cuts*/
int __stdcall cut_callback(GRBmodel *model, void *cbdata, int where, void *usrdata){

  if(where != GRB_CB_POLLING){
    struct callback_data* mydata = (struct callback_data *) usrdata;
    //clock_t current_time = clock();
    auto current = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> running_time = current - mydata->chrono_start_time;
    //double running_time = (current_time - mydata->start_time)/static_cast<double>(CLOCKS_PER_SEC) ;
    /*int error;
      error = GRBcbget(cbdata, where, GRB_CB_RUNTIME, &running_time);

      if(error){
      cerr << "ERROR: " << error << " getting runtime." << endl;
      }*/
    if(running_time > std::chrono::seconds(STOPPING_TIME)){
      GRBterminate(model);
    }
  }

  if(where == GRB_CB_MIPNODE){
    struct callback_data* mydata = (struct callback_data *) usrdata;
    int error=0;

    //cout << "found new incumbent" << endl;
    const int num_nodes = mydata->num_nodes;
    int fort_size = mydata->forts.size();
    //cout << "forts.size() is " << fort_size << endl;

    double* incumbent_solution = (double*) calloc(mydata->num_variables, sizeof(double));
    error = GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, incumbent_solution);

    std::vector<double> weights(mydata->num_variables, 0);
    for(int i=0; i<mydata->num_variables; ++i){
      weights[i] = incumbent_solution[i];
    }


    //Find the violated fort constraints



    auto time_start = std::chrono::high_resolution_clock::now();
    std::set<std::set<int> > forts;

    std::set<int> forbidden_nodes;


    mydata->our_graph.find_minimal_fort_LP(mydata->env, forts, forbidden_nodes, weights, 1.0);
    //Get forbidden nodes as intersection of previous forts

    //std::vector<double> one_weights(num_nodes, 1.0);
    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, forbidden_nodes, weights, 1.0, mydata->num_forts_variable_in); //Current best

    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0); //Current best

    //mydata->our_graph.find_max_intersection_fort_LP(mydata->env, forts, mydata->num_forts_variable_in, fort_size, weights, 1.0);
    //mydata->our_graph.find_maximal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    //mydata->our_graph.find_minimal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    //mydata->our_graph.find_forts_LP(mydata->env, forts, mydata->nodes_in_cover, fort_number_weights, 1.0);

    //Add the fort constraints for start vertices
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      //int v = mydata->our_graph.determine_facet_AB_LP(mydata->env, *it);

      int v = -1;
      if(v >= 0){
	std::vector<int> cind;
	std::vector<double> cval;
        cind.push_back(v);
        cval.push_back(1.0);
        for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
          cind.push_back(*jt);
          cval.push_back(1.0);
        }
        error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 2);
        if(error){
          cout << "ERROR at: add callback gomory constraint " << error << endl;
          return 1;
        }
      }
      else{
	//mydata->forts.push_back(*it);
	//std::set<int> new_border;
	bool is_in_border;
	std::vector<int> cind;
	//std::vector<int> cind2;
	std::vector<double> cval;
	std::set<int> a_variables;
	//cout << "adding constraint with: ";
	int min_degree = INT_MAX;
	for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
	  const int degree = mydata->our_graph.nodes[*jt].get_degree();
	  if(degree < min_degree){
	    min_degree = degree;
	  }
	}


	for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
	  //mydata->forts_node_in[*jt].push_back(mydata->forts.size()-1);
	  is_in_border = false;
	  for(int i=0; i<mydata->our_graph.nodes[*jt].get_degree(); ++i){
	    if(it->count(mydata->our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
	      //cout << "node is on border " << endl;
	      //new_border.insert(*jt);
	      is_in_border = true;
	      //mydata->num_forts_variable_in[*jt]++;
	      //cind.push_back(get_x_v_S_number(num_nodes,*jt));
	      //cind2.push_back(mydata->end_variable_start + *jt);
	      //cval.push_back(-1.0);
	      break;
	    }
	  }
	  if(is_in_border){
	    cind.push_back(*jt);
	    cval.push_back(-1.0);
	  }
	  else{
	    cind.push_back(*jt);
	    cval.push_back(-1.0/min_degree);
	  }
	}


        //cind2.push_back(mydata->end_variable_start + *jt);


	//Find out the minimum "degree" of the nodes adjacent to each node in the fort
	//fort_degree is the number of neighbors inside the fort
	//min_out_degree is the number of neighbors outside the fort
	/*
	  const int node_degree = mydata->our_graph.nodes[*jt].get_degree();
	  int fort_degree = 0;
	  int min_out_degree = num_nodes;
	  for(int i=0; i<node_degree; ++i){
	    const int neighbor = mydata->our_graph.nodes[*jt].get_adj_list_member(i);
	      //If not in fort, check neighbors degree of nodes in fort
	        if(it->count(neighbor) == 0){
		    const int neighbor_degree = mydata->our_graph.nodes[neighbor].get_degree();
		        int degree = 0;
			    for(int j=0; j<neighbor_degree; ++j){
			          if(it->count(mydata->our_graph.nodes[neighbor].get_adj_list_member(j))==1){
          degree++;
	        }
		    }
		        if(degree < min_out_degree){
			      min_out_degree = degree;
			          }
				    }
				      else{
				          fort_degree++;
					    }
					    }
					    if(fort_degree == 0){
					      fort_degree++;
					      }
					      if(min_out_degree < 2){
					      std::cerr << "ERROR: Fort is not actually a fort!" << endl;
					      }
        for(int i=0; i<it->size(); ++i){
         // cval.push_back((fort_degree > (min_out_degree-1)) ? -1.0/((double)min_out_degree-1) : -1.0/(double)fort_degree);
        }
	*/
	//cout << endl;
	//cout << "cind.size() is " << cind.size() << endl;

	//mydata->fort_borders.push_back(new_border);
	mydata->num_forts++;
	error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0);
	if(error){
	  cout << "ERROR at: add callback fort constraint " << error << endl;
	  return 1;
	}
      }
    }

    //clock_t time_end = clock();
    auto time_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> running_time = time_end - time_start;
    mydata->time_in_LP += running_time.count();
  }

  return 0;
};
//*/



int main(int argc, char* argv[]){


    if(argc != 2){
        cerr << "Error: usage is find_forcing_set <path to .edg file>" << endl;
        return 1;
    }

    GRBenv *env = NULL;
    GRBmodel *zero_forcing = NULL;
    GRBmodel *zero_forcing_linear = NULL;
    GRBmodel *col_gen = NULL;
    std::vector<int> end1;
    std::vector<int> end2;
    std::vector<double> rhs;

    int error = 0;
    std::vector<int> blank_vector;
    forcing::graph our_graph;
    forcing::edge_list edges;


  /******File Reader*****************************************************************/
    fstream file;
    string filename = argv[1];

    file.open(argv[1]);

    int num_edges;
    int num_nodes;
    int max_timesteps;
    int min_degree;
    //edge_list e_edges;
    //edge_list e2_edges;
    if(file.is_open()){

    file >> num_nodes;
    file >> num_edges;
    max_timesteps = 8;//num_nodes*3/4;  //THIS BOUND IS ONLY VALID FOR CUBIC GRAPHS
    }
    else{
      cerr << "ERROR: could not open file: " << filename << endl;
    }
    forcing::node blank_node;
    std::vector<forcing::node> nodes(num_nodes, blank_node);

    for(int k=0; k<num_edges;++k){
      int i,j;

      file >> i;
      file >> j;
      end1.push_back(i);
      end2.push_back(j);
      nodes[i].add_to_adj_list(j);
      nodes[j].add_to_adj_list(i);

      //e_edges.add_edge(i,j);
      edges.add_edge(i,j);
      edges.add_edge(j,i);
      our_graph.edges.add_edge(i,j);
      our_graph.edges.add_edge(j,i);
    }
    file.close();

    our_graph.nodes = nodes;
    num_edges = edges.get_size();

    //find min degree
    min_degree = num_nodes;
    for(int i=0; i<num_nodes; ++i){
      if(nodes[i].get_degree() < min_degree){
	min_degree = nodes[i].get_degree();
      }
    }


    auto setup_start_time = std::chrono::high_resolution_clock::now();
	//Gurobi model********************************************************************************

    //Create the Gurobi environment
    error = GRBloadenv(&env, NULL);
    if(error){
      cout << "ERROR at: load environment " << error << endl;
      return 1;
    }
    //error = GRBsetintparam(GRBgetenv(zero_forcing), "UpdateMode", 1);
    //const int num_variables = num_nodes*3+num_edges;


    /*Create the model*******************************************************/
    const int num_variables = 2*num_nodes+ num_edges;
    const int s_start = 0;
    const int x_start = s_start + num_nodes;
    const int y_start = x_start + num_nodes;

    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<x_start; ++ i){
      if(PROBLEM_TYPE == 1){
        types[i] = GRB_BINARY;
      }
      if(PROBLEM_TYPE == 0){
        types[i] = GRB_CONTINUOUS;
      }

      objective_coeffs[i] = 1.0;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1;
    }

    for(int i=x_start; i<y_start; ++ i){
      if(PROBLEM_TYPE == 1){
        types[i] = GRB_INTEGER;
      }
      if(PROBLEM_TYPE == 0){
        types[i] = GRB_CONTINUOUS;
      }
      objective_coeffs[i] = 0;
      lower_bounds[i] = 0;
      upper_bounds[i] = max_timesteps;
    }

    for(int i=y_start; i<num_variables; ++ i){
      if(PROBLEM_TYPE == 1){
        types[i] = GRB_BINARY;
      }
      if(PROBLEM_TYPE == 0){
        types[i] = GRB_CONTINUOUS;
      }
      objective_coeffs[i] = 0.0;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1.0;
    }

    error = GRBnewmodel(env, &zero_forcing, "forcing", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
      return 1;
    }
    free(objective_coeffs);
    free(lower_bounds);
    free(upper_bounds);
    free(types);

    error = GRBupdatemodel(zero_forcing);
    if(error){
      cout << "ERROR at: first model update" << endl;
      return 1;
    }

    //Now add the constraints

  /*Add a maximal set of disjoint forts*
  {
  bool stopper = false;
  std::set<std::set<int> > forts;
  std::set<int> forbidden_nodes;
  int counter = 0;
  while(!stopper){
    std::vector<double> weights(num_nodes,1.0);
    int size_change = -1*forts.size();
    cout << "forbidden nodes is " << endl;
    for(std::set<int>::iterator jt=forbidden_nodes.begin(); jt!=forbidden_nodes.end(); ++jt){
        cout << *jt << ",";
      }
      cout << endl;
    our_graph.find_forts_LP(env, forts, forbidden_nodes, weights, (double) num_nodes);
    size_change += forts.size();
    if(size_change < 1){
      stopper = true;
    }
  }
    cout << "number of disjoint forts is " << forts.size() << endl;

    //Add fort constraints for start variables
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      //mydata.forts.push_back(*it);
      //std::set<int> new_border;
      //bool is_in_border;
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;

      for(std::set<int>::iterator jt=it->begin(); jt!=it->end(); ++jt){
        cind.push_back(*jt);
        cval.push_back(-1.0);

      }

      //mydata.fort_borders.push_back(new_border);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0, NULL);
      if(error){
        std::cerr << "ERROR: " << error << " disjoint forts constraint" << endl;
        return 1;
      }

    }

  } //*/


    //Constraint 1: each vertex has some method of forcing
    for(int i=0; i<num_nodes; ++i){
      std::vector<int> cind;
      std::vector<double> cval;
      cind.push_back(s_start + i);
      cval.push_back(1.0);
      const int node_degree = our_graph.nodes[i].get_degree();
      for(int j = 0; j < node_degree; ++j){
        cind.push_back(y_start + our_graph.edges.get_edge_number(our_graph.nodes[i].get_adj_list_member(j),i));
        cval.push_back(1.0);
      }


      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 1.0, NULL);
      if(error){
        cout << "ERROR: " << error << " at constraint 1." << endl;
        return 1;
      }
    }

    //Constraint 2: time constraints for forcing
    for(int i=0; i<num_edges; ++i){
      std::vector<int> cind;
      std::vector<double> cval;
      cind.push_back(y_start + i);
      cval.push_back((double) max_timesteps + 1.0);
      //cind.push_back(y_start +our_graph.edges.get_edge_number(our_graph.edges.get_end2(i), our_graph.edges.get_end1(i)));
      //cval.push_back((double) max_timesteps - 1.0);

      const int edge_end1 = our_graph.edges.get_end1(i);
      const int edge_end2 = our_graph.edges.get_end2(i);

      cind.push_back(x_start + edge_end1);
      cval.push_back(1.0);
      cind.push_back(x_start + edge_end2);
      cval.push_back(-1.0);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, (double) max_timesteps, NULL);

      if(error){
        cout << "ERROR at: add loop prevention constraints " << error << endl;
        return 1;
      }
      cind.clear();
      cval.clear();
      cind.push_back(y_start + i);
      cval.push_back((double) max_timesteps + 1.0);
      cind.push_back(x_start + edge_end1);
      cval.push_back(1.0);
      cind.push_back(x_start + edge_end2);
      cval.push_back(-1.0);
      //Now add for 2nd neighbors
      const int neighbor_degree = our_graph.nodes[edge_end1].get_degree();
      for(int j=0; j<neighbor_degree; ++j){
        const int second_neighbor = our_graph.nodes[edge_end1].get_adj_list_member(j);
        if(second_neighbor != edge_end2){
          cind[1] = x_start + second_neighbor;
          error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, (double) max_timesteps, NULL);
          if(error){
            cout << "ERROR at: add loop prevention 2nd neighbor constraints " << error << endl;
            return 1;
          }
        }
      }
    }

  error = GRBupdatemodel(zero_forcing);
    //callback_data* blank_data;
    //Set Callback function for adding cuts
  callback_data mydata;
  mydata.previous_size = 0;
  mydata.counter = 0;
  mydata.min_degree = min_degree;
  //mydata.edges = edges;
  //mydata.nodes = nodes;
  mydata.num_forts_variable_in = std::vector<int>(num_variables, 0);
  mydata.num_nodes = num_nodes;
  mydata.num_forts = 0;
  mydata.num_variables = num_variables;
  //mydata.our_graph = our_graph;
  mydata.env = env;
  mydata.our_graph.nodes = our_graph.nodes;
  const int edge_list_size = our_graph.edges.get_size();
  for(int i=0; i<edge_list_size; ++i){
    mydata.our_graph.edges.add_edge(our_graph.edges.get_end1(i), our_graph.edges.get_end2(i));
  }
  mydata.forts_node_in.resize(num_nodes);

  error = GRBsetcallbackfunc(zero_forcing, cut_callback, (void *) &mydata);
  if(error){
	  cout << "ERROR at: GRBsetcallbackfunc" << error << endl;
	  return 1;
	}

  error = GRBsetintparam(GRBgetenv(zero_forcing), "Threads", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "Symmetry", 2);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "MIPFocus", 1);
  error = GRBupdatemodel(zero_forcing);
    auto setup_end_time = std::chrono::high_resolution_clock::now();
    auto start_time = std::chrono::high_resolution_clock::now();
    //clock_t end_time = clock();
    clock_t time_start = clock();
    auto chrono_time_start = std::chrono::high_resolution_clock::now();
    mydata.start_time = time_start;
    mydata.chrono_start_time = chrono_time_start;

  error = GRBoptimize(zero_forcing);
  //clock_t time_end = clock();
      auto end_time = std::chrono::high_resolution_clock::now();
  //double total_time = (time_end - time_start)/static_cast<double>(CLOCKS_PER_SEC);
  std::chrono::duration<double> total_time = end_time - start_time;
std::chrono::duration<double> setup_time = setup_end_time - setup_start_time;
	double IPtime;
	error = GRBgetdblattr(zero_forcing, GRB_DBL_ATTR_RUNTIME, &IPtime);
  cout << "Running time is: " << IPtime << endl;

  double optimum;
  error = GRBgetdblattr(zero_forcing, GRB_DBL_ATTR_OBJVAL, &optimum);
  cout << "ZFS size was: " << optimum << endl;
  double LP_bound;
  error = GRBgetdblattr(zero_forcing, GRB_DBL_ATTR_OBJBOUND, &LP_bound);

    ofstream fout;
  string basename = filename;
	basename = basename.substr(basename.find_last_of('/')+1);

	fout.open("infectionIP_setup_results_chrono.txt", ios::app);
	fout.precision(10);
	fout << endl << endl;
	fout << "Results for: " << basename << endl;
  fout << "ZFS size was: " << optimum << endl;
  fout << "Lower bound was: " << LP_bound << endl;
  fout << "Gurobi reported time was: " << IPtime << endl;
  fout << "Self-calculated (chrono) time was: " << total_time.count() << endl;
fout << "Setup time was: " << setup_time.count() << endl;
  fout << endl;
	fout.close();


  return 0;
  }


