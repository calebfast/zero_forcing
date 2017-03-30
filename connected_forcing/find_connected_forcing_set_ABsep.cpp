
/*********************************************************************************
This program uses the C API for Gurobi to solve the connected zero-forcing problem using
the fort cover IP and A,B-separation connectivity constraints

Author: Caleb Fast
Date: 2017

Compilation: g++ -std=c++11 -o find_connected_forcing_set_ABsep find_connected_forcing_set_ABsep.cpp -L<path to gurobi libraries> -l<gurobi###> -lgurobi_c++
Example: g++ -std=c++11 -o find_connected_forcing_set_ABsep find_connected_forcing_set_ABsep.cpp -L/opt/gurobi651/linux64/lib -lgurobi65 -lgurobi_c++

Usage: ./find_connected_forcing_set_ABsep <path to .edg file>

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
#include<stdlib.h>
#include<algorithm>

#include "../forcing_lib/forcing_lib.hpp"
#include<chrono>

/**INCLUDE PATH TO gurobi_c++.h**/
#include "<path to gurobi_c++.h>"
//Examples
//For my Laptop
//#include "/opt/gurobi651/linux64/include/gurobi_c++.h"

#define NUM_FOR_CG_CUT 3
#define SEP_START_TIME 0
#define FLOAT_TOL 0.00001
#define MIN_CLOSURE_SIZE 1
#define CUTOFF 0.9
#define STOPPING_TIME 7200
#define WAVE_STOP 10  //Controls the maximum size of sets that we run wavefront on to get cuts

using namespace std;
using namespace forcing;

int get_x_v_S_number(const int num_nodes, const int index){

  return index;
};

//Variable_in_forts should have 2*num_nodes empty elements, as should dominating vertices
void get_dominating_vertices(const int num_nodes, std::vector<std::set<int> >& forts, std::vector<std::set<int> >& fort_borders, std::vector<std::set<int> >& variable_in_forts, std::vector<node_subset>& closures, std::vector<std::vector<int> >& dominating_vertices){

  //Find the forts that contain each vertex
  for(int i=0; i<num_nodes; ++i){
    int counter = 0;

    for(std::vector<std::set<int> >::iterator it = fort_borders.begin(); it != fort_borders.end(); ++it){
      if(it->count(i) == 1){
        //cout << "variable " << i << " is in fort " << counter << endl;
        variable_in_forts[i].insert(counter);
      }
      //cout << " counter is " << counter << endl;
      //cout << "forts.size() is " << forts.size() << endl;
      //cout << "fort_borders.size() is " << fort_borders.size() << endl;
      counter++;
    }
  }

  //Find the forts that contain each subset variable
  for(int i=0; i < num_nodes; ++i){
    int counter = 0;
    for(std::vector<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
        if(closures[i].is_in_closure(*jt)){
          //cout << "variable " << num_nodes+i << " is in fort " << counter << endl;
          variable_in_forts[num_nodes+i].insert(counter);
          break;
        }
      }
      counter++;
    }
  }


  //Find which nodes dominate other nodes
  for(int i=0; i<2*num_nodes; ++i){
    for(int j=0; j<num_nodes; ++j){
      if(i != j){
      bool dominated = true;
      for(std::set<int>::iterator it =variable_in_forts[i].begin(); it!= variable_in_forts[i].end(); ++it){
        if(variable_in_forts[j].count(*it) != 1){
          dominated = false;
          break;
        }
      }
      if(dominated){
        dominating_vertices[j].push_back(i);
      }
    }
    }
  }

  //Find which subset variables dominate other subset variables
  for(int i=0; i<num_nodes; ++i){
    for(int j=0; j<num_nodes; ++j){
    if(i != j){
      bool dominated = true;
      for(std::set<int>::iterator it =variable_in_forts[num_nodes+i].begin(); it!= variable_in_forts[num_nodes+i].end(); ++it){
        if(variable_in_forts[num_nodes+j].count(*it) != 1){
          dominated = false;
          break;
        }
      }
      if(dominated){
        dominating_vertices[num_nodes+j].push_back(num_nodes+i);
      }
    }
    }
  }

};

double get_min_set_cover_number(GRBenv* env, const int num_nodes, std::vector<std::set<int> >& forts, std::vector<std::set<int> >& fort_borders, std::vector<node_subset>& closures, std::set<int>& nodes_in_cover){

  //First set up set cover IP
  //GRBenv *env2;
  //GRBloadenv(&env2, NULL);
  GRBmodel *model;
  int error;
  const int num_variables = num_nodes*2;
  char* types = (char*) calloc(num_variables, sizeof(char));
  double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
  double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
  double* lower_bounds = (double*) calloc(num_variables, sizeof(double));
  for(int i=0; i<num_variables; ++i){
    types[i] = GRB_BINARY;
    objective_coeffs[i] = 1.0;
    if(i >= num_nodes){
      objective_coeffs[i] = closures[i-num_nodes].get_weight();
    }
    upper_bounds[i] = 1.0;
    lower_bounds[i] = 0.0;
  }

  //Get dominated vertices
  std::vector<int> blank_vector;
  std::set<int> blank_set;
  std::vector<std::set<int> > variable_in_forts(2*num_nodes, blank_set);
  std::vector<std::vector<int> > dominating_vertices(2*num_nodes, blank_vector);
  get_dominating_vertices(num_nodes, forts, fort_borders, variable_in_forts, closures, dominating_vertices);

  /*set dominated variables to 0
  int counter = 0;
  for(int i=0; i<num_variables; ++i){
    //cout << "variable " << i << " dominates " << dominating_vertices[i].size() << " other variables." << endl;
    for(std::vector<int>::iterator jt = dominating_vertices[i].begin(); jt != dominating_vertices[i].end(); ++jt){
      if(*jt > i){
        upper_bounds[i] = 0.0;
        counter++;
      }
    }
  }*/
  //cout << counter << " nodes were dominated" << endl;

  error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model for intersection" << error << endl;
      return 1;
    }
    free(objective_coeffs);
    free(lower_bounds);
    free(upper_bounds);
    free(types);

    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: first model update" << endl;
      return 1;
    }

  /*This constraint just for testing*/

  /*Constraint 3: a variables are 0 if they are part of another positive a variable */
  for(int i=0; i<num_nodes; ++i){
    std::vector<int> cind;
    std::vector<double> cval;
    for(int j=0; j<num_nodes; ++j){
      if(closures[i].is_in_closure(j)){
        cind.push_back(num_nodes+j);
        cval.push_back(1.0);
      }
    }
    error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 1, NULL);
      if(error){
        std::cerr << "ERROR: " << error << "at a variable disjointness constraint" << endl;
        return 1;
      }
  }

      /* start variables are only nonzero if they are outside of a variable closures*/
    for(int i=0; i<num_nodes; ++i){
      for(int j=0; j<num_nodes; ++j){
        if(closures[i].is_in_closure(j)){
          std::vector<int> cind;
          std::vector<double> cval;
          cind.push_back(j);
          cval.push_back(1.0);
          cind.push_back(num_nodes+i);
          cval.push_back(1.0);
          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 1, NULL);
          if(error){
            std::cerr << "ERROR: " << error << "no start in closures" << endl;
            return 1;
          }
        }
      }
    }


  {
    std::vector<int> cind;
    std::vector<double> cval;
    for(int i=0; i<num_nodes; ++i){
      cind.push_back(num_nodes+i);
      cval.push_back(1.0);
    }
    error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
        cout << "ERROR in intersection cuts at: add constraint 1: " << error << endl;
        return 1;
      }
  }
  //*/

  //Now add the constraints
    const int num_forts = forts.size();
    for(int i=0; i<num_forts; ++i){
      std::vector<int> cind;
      std::vector<double> cval;
      for(int j=0; j<num_variables; ++j){
        if(variable_in_forts[j].count(i)==1){
          cind.push_back(j);
          cval.push_back(1.0);
        }
      }
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
        cout << "ERROR in intersection cuts at: add constraint 1: " << error << endl;
        return 1;
      }
    }

    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: second model update" << endl;
      return 1;
    }

    //cout << "optimizing model " << endl;
    GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    GRBsetintparam(GRBgetenv(model), "Threads", 1);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: first model update" << endl;
      return 1;
    }
    GRBoptimize(model);
    if(error){
      cout << "ERROR at: optimize" << endl;
      return 1;
    }

    double cover_number = -1;
    GRBgetdblattr(model, "ObjVal", &cover_number);
    for(int i=0; i<num_variables; ++i){
      double sol;
      GRBgetdblattrelement(model, "X", i, &sol);
      if(i < num_nodes){
        if(sol >= 1-FLOAT_TOL){
          nodes_in_cover.insert(i);
        }
      }
      if(i >= num_nodes){
        if(sol >= 1-FLOAT_TOL){
          for(int j=0; j<num_nodes; ++j){
            if(closures[i-num_nodes].is_in_closure(j)){
              nodes_in_cover.insert(j);
            }
          }
        }
      }
    }
    GRBfreemodel(model);
    //GRBfreeenv(env2);
    return cover_number;

};


void graph_viz_solution(std::vector<double>& node_labels, std::vector<double>& node_labels2 , forcing::graph& our_graph){

  ofstream gv_file;

  gv_file.open("LP_graph_cutoff.gv");
    //Begin by writing preamble

    if(gv_file.is_open()){
    gv_file << "##Command to generate layout is: neato -Tpng thisfile > thisfile.png" << endl;

    gv_file << "digraph LP_graph" << " {" << endl;

    //Write the nodes
    for(int i=0; i < node_labels.size(); i++){
      gv_file << i << " [label=\"" << std::setprecision(2) << i << ":" << node_labels[i] << ":" << node_labels2[i] << "\"];" << endl;
    }

    //Write the edges
    for(int i=0; i < our_graph.edges.get_size(); i++){
      gv_file << our_graph.edges.get_end1(i) << "->" << our_graph.edges.get_end2(i) << endl;
    }

    gv_file << "overlap=false;" << endl;
    gv_file << "label=\"" << "\";" << endl;
    gv_file << "fontsize=10" << endl;
    gv_file << "}" << endl;

    gv_file.close();

    }
    else{
    cerr << "ERROR: Could not open " << endl;
    }
}

void get_union_of_sets(std::vector<std::set<int> >& sets, std::set<int>& union_of_sets){

  for(std::vector<std::set<int> >::iterator it = sets.begin(); it != sets.end(); ++it){
    for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
      union_of_sets.insert(*jt);
    }
  }

};



bool check_connectivity(std::set<int>& chosen_nodes, std::vector<int>& positive_edges, std::set<int>& min_cut, std::set<int>& component, forcing::graph& our_graph){

    forcing::graph subgraph;
    forcing::edge_list subgraph_cut;
    std::unordered_map<int, int> reverse_node_map;
    our_graph.get_subgraph(chosen_nodes, subgraph, subgraph_cut, reverse_node_map);
    const int num_nodes = subgraph.nodes.size();
    const int num_edges = subgraph.edges.get_size();
    std::set<int> candidates;
    for(int i=0; i < num_nodes; ++i){
      candidates.insert(i);
    }

    std::set<int> subgraph_component;
    subgraph.get_single_component(subgraph_component, candidates);

    for(std::set<int>::iterator it = subgraph_component.begin(); it != subgraph_component.end(); ++it){
      component.insert(reverse_node_map[*it]);
    }

    if(subgraph_component.size() == num_nodes){
      return true;
    }
    else{
      return false;
    }
    /*
    std::vector<double> adj_matrix(num_nodes*num_nodes, 0);

    for(int i=0; i<num_edges; ++i){
      const int end1 = subgraph.edges.get_end1(i);
      const int end2 = subgraph.edges.get_end2(i);
      adj_matrix[end1*num_nodes+end2] = 1;
    }

    std::set< std::set<int> > cuts_under_one;
    double min_cut_value;
    forcing::node_identification_min_cut(num_nodes,adj_matrix, min_cut, min_cut_value, cuts_under_one);

    if(min_cut_value < 1 - FLOAT_TOL){
      return false;
    }
    else{
      return true;
    }
    */
};

/*callback function for seperating cuts*/
int __stdcall cut_callback(GRBmodel *model, void *cbdata, int where, void *usrdata){

  if(where != GRB_CB_POLLING){
    struct callback_data* mydata = (struct callback_data *) usrdata;
    //clock_t current_time = clock();
    //double running_time = (current_time - mydata->start_time)/static_cast<double>(CLOCKS_PER_SEC);
        	auto current = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> running_time = current - mydata->chrono_start_time;
    /*if(running_time > STOPPING_TIME){
      GRBterminate(model);
    }*/
    		if(running_time > std::chrono::seconds(STOPPING_TIME)){
	      GRBterminate(model);
      }
  }

  if(where == GRB_CB_MIPSOL){
    struct callback_data* mydata = (struct callback_data *) usrdata;
    int error=0;

    //cout << "found new incumbent" << endl;
    const int num_nodes = mydata->num_nodes;
    int fort_size = mydata->forts.size();
    //cout << "forts.size() is " << fort_size << endl;

    double* incumbent_solution = (double*) calloc(mydata->num_variables, sizeof(double));
    error = GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, incumbent_solution);


    std::vector<int> num_forts_node_in(num_nodes, 0);
    std::set<int> forbidden_nodes;
    std::vector<double> weights(mydata->num_variables, 0);
    for(int i=0; i<mydata->num_variables; ++i){
      weights[i] = incumbent_solution[i];
      if(i < num_nodes && incumbent_solution[i] > 1-FLOAT_TOL){
        forbidden_nodes.insert(i);
      }
    }
    /*
    for(int i=0; i<num_nodes; ++i){
      for(int j=0; j<mydata->node_subsets.size(); ++j){
        if(mydata->node_subsets[j].is_in_closure(i)){
          weights[i] += incumbent_solution[mydata->a_start+j];
        }
      }
    }*/
    //Find the end nodes of the current incumbent


    //Find the violated fort constraints


    /*std::vector<double> fort_number_weights(mydata->num_variables, 0);
    for(int i=0; i<mydata->num_variables; ++i){
      fort_number_weights[i]  = weights[i] + static_cast<double>(mydata->num_forts_variable_in[i])/(100.0*fort_size);
    }*/
    //clock_t time_start = clock();
    auto time_start = std::chrono::high_resolution_clock::now();
       std::set<std::set<int> > forts;
    //mydata->our_graph.find_fort_closure(forts, weights);


    //mydata->our_graph.find_fort_closure(forts, end_weights);
    //mydata->our_graph.find_fort_greedy_style_2(forts, weights);
    //mydata->our_graph.find_fort_greedy_style_2(forts, end_weights);
    mydata->our_graph.find_minimal_fort_LP(mydata->env, forts, forbidden_nodes, weights, 1.0);
    //mydata->our_graph.find_multiple_minimal_border_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);

    //Get forbidden nodes as intersection of previous forts
    //std::set<int> temp1(*(mydata->forts.end()-1));
    //std::set<int> temp2(*(mydata->forts.end()-2));
    //std::set<int> fort3(*(mydata->forts.end()-3));
    //std::set<int> fort4(*(mydata->forts.end()-4));
    //std::set<int> fort5(*(mydata->forts.end()-5));
    //std::set<int> temp1;
    //std::set<int> temp2;
    //set_intersection(fort1.begin(),fort1.end() , fort2.begin(), fort2.end(), std::inserter(temp1, temp1.begin()));
    //set_intersection(fort3.begin(),fort3.end() , fort4.begin(), fort4.end(), std::inserter(temp2, temp2.begin()));
    //set_intersection(temp1.begin(),temp1.end() , temp2.begin(), temp2.end(), std::inserter(forbidden_nodes, forbidden_nodes.begin()));

    //std::vector<double> one_weights(num_nodes, 1.0);
    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, forbidden_nodes, weights, 1.0, mydata->num_forts_variable_in); //Current best
    //forbidden_nodes.clear();
    /*for(int i=0; i<num_nodes; ++i){
      if(weights[i] >1-FLOAT_TOL){
        forbidden_nodes.insert(i);
      }
    }*/
    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, forbidden_nodes, fort_number_weights, 1.0, mydata->num_forts_variable_in); //Current best
    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, forbidden_nodes, weights, 1.0, mydata->num_forts_variable_in);
    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, forbidden_nodes, end_weights, 1.0, mydata->num_forts_variable_in); //Current best

    /*if(mydata->time_in_LP > mydata->previous_time + 1){
      cout << "time finding forts is " << mydata->time_in_LP << endl;
      cout << "current number of forts is " << mydata->forts.size() << endl;
      mydata->previous_time = mydata->time_in_LP;
    }*/
    //mydata->nodes_in_cover.clear();
    //int temp_forts_size = forts.size();

    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0); //Current best

    //mydata->our_graph.find_max_intersection_fort_LP(mydata->env, forts, mydata->num_forts_variable_in, fort_size, weights, 1.0);
    //mydata->our_graph.find_maximal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    //mydata->our_graph.find_minimal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    /*mydata->our_graph.find_forts_LP(mydata->env, forts, mydata->nodes_in_cover, fort_number_weights, 1.0);
    if(forts.size() == temp_forts_size){
      cout << "running second IP" << endl;
      mydata->nodes_in_cover.clear();
      mydata->our_graph.find_forts_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
      if(forts.size() == temp_forts_size){
        cout << "second IP found nothing" << endl;
      }
    }*/
    //mydata->our_graph.find_forts_LP(mydata->env, forts, mydata->nodes_in_cover, fort_number_weights, 1.0);

    //Add the fort constraints for start vertices
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      //mydata->forts.push_back(*it);
      //std::set<int> new_border;
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      std::set<int> a_variables;
      //cout << "adding constraint with: ";
      for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
        bool on_border = false;
         for(int i=0; i<mydata->our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(mydata->our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            cind.push_back(get_x_v_S_number(num_nodes,*jt));
            //cind2.push_back(mydata->end_variable_start+*jt);
            cval.push_back(-1.0);
            on_border = true;
            break;
          }
        }
        if(!on_border){
          cind.push_back(get_x_v_S_number(num_nodes,*jt));
          cval.push_back(-1.0/mydata->min_degree);
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

      for(int j=0; j<mydata->node_subsets.size(); ++j){
        if(mydata->node_subsets[j].is_in_closure(*jt)){
          //cout << "node is in closure " << endl;
          a_variables.insert(j);
        }
       }


      }



      if(error){
	cout << "ERROR beffore: add callback fort constraint " << error << endl;
	return 1;
      }
      //mydata->fort_borders.push_back(new_border);
      mydata->num_forts++;
      error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0);
      if(error){
	cout << "ERROR at: add callback fort constraint " << error << endl;
	return 1;
      }

    }
    auto time_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> lp_time = time_end - time_start;
    mydata->time_in_LP += lp_time.count();

    if(forts.size() == 0){
      //clock_t connect_time_begin = clock();
      auto connect_time_start = std::chrono::high_resolution_clock::now();
      //Check for connectivity
    {
    std::set<int> chosen_nodes;
    std::set<int> separator;
    std::set<int> min_cut;
    std::vector<int> positive_edges;
    std::vector<int> cind;
    std::vector<double> cval;
    int counter = 0;
    for(int i=0; i<num_nodes; ++i){
      if(incumbent_solution[i] > 1-FLOAT_TOL){
        chosen_nodes.insert(i);
        counter++;
        //cind.push_back(i);
        //cval.push_back(-1.0);
      }
      else{
        separator.insert(i);
        //cind.push_back(i);
        //cval.push_back(1.0);
      }
    }
    for(int i=0; i<mydata->num_edges; ++i){
      const int end1 = mydata->our_graph.edges.get_end1(i);
      const int end2 = mydata->our_graph.edges.get_end2(i);
      if(incumbent_solution[end1] > 1-FLOAT_TOL && incumbent_solution[end2] > 1-FLOAT_TOL){
        positive_edges.push_back(i);
      }
    }
    bool is_connected = true;
    std::vector<std::set<int> > components;
    std::set<int> component;
    is_connected = check_connectivity(chosen_nodes, positive_edges, min_cut, component, mydata->our_graph);
    if(is_connected){
      cout << "SOLUTION IS CONNECTED" << endl;
    }
    if(!is_connected){
      cout << "NOT CONNECTED: ADDING CONSTRAINT" << endl;



    /*for(std::vector<int>::iterator it = component.begin(); it != component.end(); ++it){
      const int current_node = *it;
      const int node_degree = mydata->our_graph.nodes[current_node].get_degree();
      for(int i=0; i<node_degree; ++i){
        const int neighbor = mydata->our_graph.nodes[current_node].get_adj_list_member(i);
        if(component.count(neighbor) == 0){
          separator.insert(neighbor);
        }
      }
    }*/

    const int a = *component.begin();
    int b;
    for(std::set<int>::iterator it = chosen_nodes.begin(); it != chosen_nodes.end(); ++it){
      if(component.count(*it) < 1){
        b = *it;
        break;
      }
    }
    mydata->our_graph.get_minimal_separator(separator, components, a, b);

    cind.push_back(a);
    //cout << "a = " << a << endl;
    cval.push_back(-1.0);
    cind.push_back(b);
    //cout << "b = " << b << endl;
    cval.push_back(-1.0);
    //cout << "Separator is: ";
    for(std::set<int>::iterator it = separator.begin(); it != separator.end(); ++it){
      //cout << *it << ", ";
      cind.push_back(*it);
      cval.push_back(1.0);
    }
    //cout << endl;
      mydata->counter++;
      error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, -1.0);
      if(error){
        cerr << "ERROR " << error << " at connectivity constraint" << endl;
        return 1;
      }
    }

    /*if(!is_connected){
      error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, (1.0 - counter));
      if(error){
        cerr << "ERROR " << error << " at connectivity constraint" << endl;
        return 1;
      }
    }*/

    }

      //clock_t connect_time_end = clock();
      auto connect_time_end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> connect_time = connect_time_end - connect_time_start;
      mydata->time_in_connectivity += connect_time.count();
    }




    /*Find which vertices are dominated
    std::set<int> blank_set;
    std::vector<int> blank_vector;
    std::vector<std::set<int> > variable_in_forts(2*num_nodes, blank_set);
    std::vector<std::vector<int> > dominating_vertices(2*num_nodes, blank_vector);
    get_dominating_vertices(num_nodes, mydata->forts, mydata->fort_borders, variable_in_forts, mydata->node_subsets, dominating_vertices);


  for(int i=0; i<2*num_nodes; ++i){
    error = GRBsetdblattrelement(model, "UB", i, 1.0);
    if(error){
      cout << "ERROR at: setting bounds " << error << endl;
        return 1;
    }
  }


  for(int i=0; i<num_nodes; ++i){
    for(std::vector<int>::iterator jt = dominating_vertices[i].begin(); jt != dominating_vertices[i].end(); ++jt){
      if(*jt > i){
        error = GRBsetdblattrelement(model, "UB", *jt, 1.0);
        if(error){
      cout << "ERROR at: setting bounds " << error << endl;
        return 1;
    }
      }
    }
  }
  for(int i=num_nodes; i<2*num_nodes; ++i){
    for(std::vector<int>::iterator jt = dominating_vertices[i].begin(); jt != dominating_vertices[i].end(); ++jt){
      if(*jt > i){
        error = GRBsetdblattrelement(model, "UB", *jt, 1.0);
        if(error){
      cout << "ERROR at: setting bounds " << error << endl;
        return 1;
    }
      }
    }
  }
GRBupdatemodel(model);

  int num_not_dominated = 0;
  for(int i=0; i<2*num_nodes; ++i){
    double temp;
    error = GRBgetdblattrelement(model, "UB", i, &temp);
    if(error){
      cout << "ERROR at: setting bounds " << error << endl;
        return 1;
    }
    if(temp >= 1.0 -FLOAT_TOL){
      cout << i << " is not dominated " << endl;
      num_not_dominated++;
    }
  }
  cout << "num_forts is " << mydata->forts.size() << endl;
  cout << "num_not_dominated is " << num_not_dominated << endl;
*/

    /*Add the CG cuts for nodes in 5 forts.*
    int status = 0;
    //Find 5 forts that give a cut
    //while(status == 0){
    if(mydata->forts.size() - mydata->forbidden_fort_indices.size() > NUM_FOR_CG_CUT){
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      std::vector<int> fort_indices;
      status = mydata->our_graph.find_fort_intersection_cuts(mydata->env, num_nodes, NUM_FOR_CG_CUT, mydata->forts, mydata->forbidden_fort_indices, fort_indices, weights);
      if(status == 0){
        //std::vector<double> node_coefficients(num_nodes, 0.0);
        std::set<int> node_indices;

        //Fill cind and cval vectors
        std::vector<double> fort_coefficients(2*num_nodes, 0);
    for(int i=0; i<num_nodes; ++i){
      for(int j=0; j<fort_indices.size(); ++j){
        if(mydata->forts[fort_indices[j]].count(i) == 1){

          fort_coefficients[mydata->a_start+i] += 1;
          //node_coefficients[i] += 1.0;
          //std::pair<std::set<int>::iterator, bool> temp;
          node_indices.insert(i);
          for(int k=0; k<mydata->our_graph.nodes[i].get_degree(); ++k){
          if(mydata->forts[fort_indices[j]].count(mydata->our_graph.nodes[i].get_adj_list_member(k)) != 1){
            //cout << "node is on border " << endl;
            fort_coefficients[i] += 1;

            //cind.push_back(get_x_v_S_number(num_nodes,i));
            //cval.push_back(1.0);
            break;
          }
        }
        }
      }
    }




    for(int i=0; i<2*num_nodes; ++i){
      if(fort_coefficients[i] > FLOAT_TOL){
        cind.push_back(i);
        cval.push_back(ceil(fort_coefficients[i]/2.0));
      }
    }

    for(int i=0; i<num_nodes; ++i){
      if(node_indices.count(i) == 1){

        //cind.push_back(i);
        //cind2.push_back(mydata->end_variable_start+i);
        //cval.push_back(1.0);
        //cind.push_back(mydata->a_start + i);
        //cind2.push_back(mydata->ae_start + i);
        //cval.push_back(1.0);
      }
    }
    //Add the cut
    error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, ceil((double)NUM_FOR_CG_CUT/2.0));
    if(error){
      cout << "ERROR at: add callback fort constraint " << error << endl;
      return 1;
    }

    mydata->forbidden_fort_indices.push_back(fort_indices[0]);
    //
    //cout << "forts in constraint are ";
    //for(std::vector<int>::iterator it =fort_indices.begin(); it!=fort_indices.end(); ++it){
    //  cout << *it << " ";
    //}
    //cout << endl;
  }
  else{
    //cout << "ended with status " << status << endl;
  }
  }
 // }*/

        //clock_t time_end = clock();

  }

   if(false && where == GRB_CB_MIPNODE){
    struct callback_data* mydata = (struct callback_data *) usrdata;
    int error=0;

    int status;
    error = GRBcbget(cbdata, where, GRB_CB_MIPNODE_STATUS, &status);
    if(status == GRB_OPTIMAL){
    const int num_nodes = mydata->num_nodes;


    double* incumbent_solution = (double*) calloc(mydata->num_variables, sizeof(double));
    error = GRBcbget(cbdata, where, GRB_CB_MIPNODE_REL, incumbent_solution);
    if(error){
      cout << "ERROR at: get current relaxation solution " << error << endl;
      return 1;
    }

    double running_time;
    error = GRBcbget(cbdata, where, GRB_CB_RUNTIME, &running_time);
    if(running_time > STOPPING_TIME){
      std::vector<double> node_labels(incumbent_solution, incumbent_solution+num_nodes);
      std::vector<double> node_labels2(incumbent_solution+num_nodes, incumbent_solution+2*num_nodes);
      graph_viz_solution(node_labels, node_labels2, mydata->our_graph);
      GRBterminate(model);
    }

    std::vector<int> num_forts_node_in(num_nodes, 0);
    std::set<int> forbidden_nodes;
    std::vector<double> weights(mydata->num_variables, 0);
    for(int i=0; i<mydata->num_variables; ++i){
      weights[i] = incumbent_solution[i];
    }
    for(int i=0; i<num_nodes; ++i){
     for(int j=0; j<mydata->node_subsets.size(); ++j){
        if(mydata->node_subsets[j].is_in_closure(i)){
          weights[i] += incumbent_solution[mydata->a_start+j];
        }
      }
    }

    std::set<std::set<int> > forts;
    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0, num_forts_node_in);
    //mydata->our_graph.find_maximal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);

     //Add the fort constraints for start vertices
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      mydata->forts.push_back(*it);
      std::vector<int> cind;
      std::vector<int> cind2;
      std::vector<double> cval;
      std::set<int> a_variables;
      //cout << "adding constraint with: ";
      for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
        mydata->forts_node_in[*jt].push_back(mydata->forts.size()-1);
        for(int i=0; i<mydata->our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(mydata->our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            cind.push_back(get_x_v_S_number(num_nodes,*jt));
            cind2.push_back(mydata->end_variable_start+*jt);
            cval.push_back(-1.0);
            break;
          }
        }



        //cind2.push_back(mydata->end_variable_start + *jt);


	//Find out the minimum "degree" of the nodes adjacent to each node in the fort
	//fort_degree is the number of neighbors inside the fort
	//min_out_degree is the number of neighbors outside the fort
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

	//cout << endl;
	//cout << "cind.size() is " << cind.size() << endl;

      for(int j=0; j<mydata->node_subsets.size(); ++j){
        if(mydata->node_subsets[j].is_in_closure(*jt)){
          //cout << "node is in closure " << endl;
          a_variables.insert(j);
        }
       }


      }

      for(std::set<int>::iterator kt = a_variables.begin(); kt != a_variables.end(); ++kt){
        cind.push_back(mydata->a_start + *kt);
        cind2.push_back(mydata->ae_start + *kt);
        cval.push_back(-1.0);
      }


      if(error){
	cout << "ERROR beffore: add callback fort constraint " << error << endl;
	return 1;
      }
      error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0);
      if(error){
	cout << "ERROR at: add callback fort constraint " << error << endl;
	return 1;
      }
      //error = GRBcblazy(cbdata, cind2.size(), &cind2[0], &cval[0], GRB_LESS_EQUAL, -1.0);
      if(error){
	cout << "ERROR at: add callback fort constraint " << error << endl;
	return 1;
      }
    }


     /*Add the CG cuts for nodes in 5 forts.*/
    int status = 0;
    //Find 5 forts that give a cut
    //while(status == 0){
    if(mydata->forts.size() - mydata->forbidden_fort_indices.size() > NUM_FOR_CG_CUT){
      std::vector<int> cind;
      std::vector<int> cind2;
      std::vector<double> cval;
      std::vector<int> fort_indices;
      status = mydata->our_graph.find_fort_intersection_cuts(mydata->env, num_nodes, NUM_FOR_CG_CUT, mydata->forts, mydata->forbidden_fort_indices, fort_indices, weights);
      if(status == 0){
        //std::vector<double> node_coefficients(num_nodes, 0.0);
        std::set<int> node_indices;

        //Fill cind and cval vectors
        std::vector<double> fort_coefficients(2*num_nodes, 0);
    for(int i=0; i<num_nodes; ++i){
      for(int j=0; j<fort_indices.size(); ++j){
        if(mydata->forts[fort_indices[j]].count(i) == 1){

          fort_coefficients[mydata->a_start+i] += 1;
          //node_coefficients[i] += 1.0;
          //std::pair<std::set<int>::iterator, bool> temp;
          node_indices.insert(i);
          for(int k=0; k<mydata->our_graph.nodes[i].get_degree(); ++k){
          if(mydata->forts[fort_indices[j]].count(mydata->our_graph.nodes[i].get_adj_list_member(k)) != 1){
            //cout << "node is on border " << endl;
            fort_coefficients[i] += 1;

            //cind.push_back(get_x_v_S_number(num_nodes,i));
            //cval.push_back(1.0);
            break;
          }
        }
        }
      }
    }




    for(int i=0; i<2*num_nodes; ++i){
      if(fort_coefficients[i] > FLOAT_TOL){
        cind.push_back(i);
        cind2.push_back(mydata->end_variable_start+i);
        cval.push_back(ceil(fort_coefficients[i]/2.0));
      }
    }

    for(int i=0; i<num_nodes; ++i){
      if(node_indices.count(i) == 1){

        //cind.push_back(i);
        //cind2.push_back(mydata->end_variable_start+i);
        //cval.push_back(1.0);
        //cind.push_back(mydata->a_start + i);
        //cind2.push_back(mydata->ae_start + i);
        //cval.push_back(1.0);
      }
    }
    //Add the cut
    error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, ceil((double)NUM_FOR_CG_CUT/2.0));
    if(error){
      cout << "ERROR at: add callback fort constraint " << error << endl;
      return 1;
    }
    //error = GRBcblazy(cbdata, cind2.size(), &cind2[0], &cval[0], GRB_GREATER_EQUAL, ceil((double)NUM_FOR_CG_CUT/2.0));
    if(error){
      cout << "ERROR at: add callback fort constraint " << error << endl;
      return 1;
    }

    mydata->forbidden_fort_indices.push_back(fort_indices[0]);
    //
    //cout << "forts in constraint are ";
    //for(std::vector<int>::iterator it =fort_indices.begin(); it!=fort_indices.end(); ++it){
    //  cout << *it << " ";
    //
    //cout << endl;
  }
  }
 // }*/

  }

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
    cout << "num_nodes is " << num_nodes << endl;
    max_timesteps = num_nodes*3/4;  //For now, this is the best bound I have on the number of timesteps required to force the graph
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


  /*Find branch decomposition*****************************************************
  int branch_width;
  int bnodes, bedges;
  int* bend1, *bend2, *nu;
  int throwaway = branchdecomp(num_nodes,    //*number of nodes in input graph
			   end1.size(),    //*number of edges in input graph *
			   &end1[0],     //*list of beginning of edges in input graph *
			   &end2[0],     //*list of ends of edges in input graph *
			   branch_width,   //*branch-width will return in this*
			   bnodes,  //*number of nodes in branch decomposition*
			   bedges,  //*number of edges in branch decomposition*
			   bend1,   //*list of beginning of edges in branch decomposition*
			   bend2,   //*list of ends of edges in branch decomposition*
			   nu);  //*bijection from edges of input graph to leaves of branch decomposition*

  forcing::graph bd_tree;
  for(int i =0; i < bnodes; ++i){
    forcing::node new_node;
    bd_tree.nodes.push_back(new_node);
  }
  for(int i=0; i<bedges; ++i){
    //cout << "node degree is: " << bd_tree.nodes[bend1[i]].get_degree() << endl;
    bd_tree.nodes[bend1[i]].add_to_adj_list(bend2[i]);
    //cout << "node degree is: " << bd_tree.nodes[bend1[i]].get_degree() << endl;
    //cout << endl;
    bd_tree.nodes[bend2[i]].add_to_adj_list(bend1[i]);
  }

  bd_tree.nodes[0].daddy = 0;
  std::queue<int> postdfs_order;
  forcing::dfs_order(bd_tree.nodes, postdfs_order, bnodes, num_edges/2, 0);

  //Fill the subsets for leaves
  forcing::fill_leaves(bd_tree.nodes, num_edges/2, &end1[0], &end2[0], nu);
  //Fill rest of subsets
  while(!postdfs_order.empty()){
    const int index = postdfs_order.front();
    postdfs_order.pop();
    bd_tree.nodes[index].subtree.insert(bd_tree.nodes[bd_tree.nodes[index].child1].subtree.begin(), bd_tree.nodes[bd_tree.nodes[index].child1].subtree.end());
    bd_tree.nodes[index].subtree.insert(bd_tree.nodes[bd_tree.nodes[index].child2].subtree.begin(), bd_tree.nodes[bd_tree.nodes[index].child2].subtree.end());
  }
    */


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
    //Start with number of forcing chains = to heuristic branch_width value
    //Create the model with 3N*branch_width binary variables
    const int num_variables = num_nodes + end1.size();
    const int a_start = num_nodes;
    const int vertex_time_start = num_nodes*2;
    const int edge_variable_start = num_nodes;
    const int loop_preventer_start = num_nodes*4+num_edges;

    const int ae_start = num_nodes*3;
    const int y_start = num_nodes*4;

    std::vector<node_subset> node_subsets;
    for(int i=0; i<num_nodes; ++i){
      node_subset new_subset;
      new_subset.add_to_subset(i);
      const int num_neighbors = our_graph.nodes[i].get_degree();
      for(int j=0; j<num_neighbors-1; ++j){
        new_subset.add_to_subset(our_graph.nodes[i].get_adj_list_member(j));
      }
      new_subset.set_weight();
      new_subset.update_closure(our_graph.nodes);
      node_subsets.push_back(new_subset);
    }

  //Set Callback function for adding cuts
    callback_data mydata;
    mydata.time_in_LP = 0;
    mydata.not_added = true;
    mydata.max_timesteps = max_timesteps;
    mydata.edge_variable_start = edge_variable_start;
    mydata.vertex_time_start = vertex_time_start;
    mydata.time_in_connectivity = 0;
    mydata.previous_size = 0;
  mydata.counter = 0;
  mydata.min_degree = min_degree;
  mydata.y_start = y_start;
  mydata.a_start = a_start;
  mydata.ae_start = ae_start;
  //mydata.end_variable_start = end_variable_start;
  mydata.num_edges = num_edges;
  //mydata.nodes = nodes;
  mydata.num_forts = 0;
  mydata.num_forts_variable_in = std::vector<int>(num_variables, 0);
  mydata.num_nodes = num_nodes;
  mydata.num_variables = num_variables;
  mydata.node_subsets = node_subsets;
  //mydata.our_graph = our_graph;
  mydata.env = env;
    mydata.our_graph.nodes = our_graph.nodes;
    const int edge_list_size = our_graph.edges.get_size();
    for(int i=0; i<edge_list_size; ++i){
      mydata.our_graph.edges.add_edge(our_graph.edges.get_end1(i), our_graph.edges.get_end2(i));
    }
  mydata.forts_node_in.resize(num_nodes);
  cout << "mydata->node_subset.size() is " << mydata.node_subsets.size() << endl;
  cout << "our_graph nodes.size() is " << our_graph.nodes.size() << endl;
  cout << "our_graph edge_list.size() is " << our_graph.edges.get_size() << endl;
  cout << "mydata.our_graph nodes.size() is " << mydata.our_graph.nodes.size() << endl;
  cout << "mydata.our_graph edge_list.size() is " << mydata.our_graph.edges.get_size() << endl;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));
    int* branching_priority = (int*) calloc(num_variables, sizeof(int));
    //Variables for first time stage (x_v,0 variables)
    for(int i=0; i<num_nodes; ++ i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      //Objective sums all nodes at first time stage
      objective_coeffs[i] = 1;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1;
      branching_priority[i] = 4;
    }
    for(int i=num_nodes; i<num_variables; ++ i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      //Objective sums all nodes at first time stage
      objective_coeffs[i] = 0;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1;
      branching_priority[i] = 4;
    }

    /*for(int i=a_start; i<vertex_time_start; ++ i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      //Objective sums all nodes at first time stage
      objective_coeffs[i] = (double) node_subsets[i-a_start].get_weight();
      lower_bounds[i] = 0;
      upper_bounds[i] = 1;
      branching_priority[i] =5;
    }
    for(int i=vertex_time_start; i<edge_variable_start; ++ i){
      types[i] = GRB_INTEGER;
      //types[i] = GRB_CONTINUOUS;
      //Objective sums all nodes at first time stage
      objective_coeffs[i] = 0.0;
      lower_bounds[i] = 0;
      upper_bounds[i] = max_timesteps;
      branching_priority[i] = 2;
    }
    for(int i=edge_variable_start; i<num_variables; ++ i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      //Objective sums all nodes at first time stage
      objective_coeffs[i] = 0.0;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1.0;
      branching_priority[i] =0;
    }*/



    /***********************/
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
    //set branching priority
    //error = GRBsetintattrarray(zero_forcing, "BranchPriority", 0, num_variables, branching_priority);
    if(error){
      cout << "ERROR at: set BranchPriority " << error << endl;
      return 1;
    }
    free(branching_priority);
    error = GRBupdatemodel(zero_forcing);

    if(error){
      cout << "ERROR at: write model: " << error << endl;
      return 1;
    }


 /*Constraint 1: Add fort constraints that do not contain subset variables*
  {
  std::set<std::set<int> > forts;
  for(int i=0; i<node_subsets.size(); ++i){
      std::set<int> new_fort;
      for(int j=0; j<num_nodes; ++j){
        if(!node_subsets[i].is_in_closure(j)){
          new_fort.insert(j);
        }
      }
    forts.insert(new_fort);
  }

    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      mydata.num_forts++;
      //mydata.forts.push_back(*it);
      //std::set<int> new_border;
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      std::set<int> a_variables;
      for(std::set<int>::iterator jt=it->begin(); jt!=it->end(); ++jt){
        bool on_border = false;
        for(int i=0; i<our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            //new_border.insert(*jt);
            //mydata.num_forts_variable_in[*jt]++;
            cind.push_back(get_x_v_S_number(num_nodes,*jt));
            //cind2.push_back(end_variable_start+*jt);
            cval.push_back(-1.0);
            on_border = true;
            break;
          }
        }
        if(!on_border){
          cind.push_back(get_x_v_S_number(num_nodes,*jt));
          cval.push_back(-1.0/min_degree);
        }
        for(int j=0; j<mydata.node_subsets.size(); ++j){
          if(mydata.node_subsets[j].is_in_closure(*jt)){
            //cout << "node is in closure " << endl;
            a_variables.insert(j);
          }
        }
      }


      //for(std::set<int>::iterator kt = a_variables.begin(); kt != a_variables.end(); ++kt){
      //  cind.push_back(mydata.a_start + *kt);
      //  //cind2.push_back(ae_start + *kt);
      //  //mydata.num_forts_variable_in[mydata.a_start + *kt]++;
      //  //cind2.push_back(ae_start + *kt);
      //  cval.push_back(-1.0);
      //}

      //mydata.fort_borders.push_back(new_border);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0, NULL);
      if(error){
        std::cerr << "ERROR: " << error << " disjoint forts constraint" << endl;
        return 1;
      }

    }
  }
  //*/


  /*Add a maximal set of disjoint forts*/
  {
  bool stopper = false;
  std::set<std::set<int> > forts;
  std::set<int> forbidden_nodes;
  int counter = 0;
  while(!stopper){
    std::vector<double> weights(num_nodes,1.0);
    int size_change = -1*forts.size();
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
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      std::set<int> a_variables;
      for(std::set<int>::iterator jt=it->begin(); jt!=it->end(); ++jt){
        bool on_border = false;
        for(int i=0; i<our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            //new_border.insert(*jt);
            //mydata.num_forts_variable_in[*jt]++;
            cind.push_back(get_x_v_S_number(num_nodes,*jt));
            //cind2.push_back(end_variable_start+*jt);
            cval.push_back(-1.0);
            on_border = true;
            break;
          }
        }
        if(!on_border){
          cind.push_back(get_x_v_S_number(num_nodes,*jt));
          cval.push_back(-1.0/min_degree);
        }
        for(int j=0; j<mydata.node_subsets.size(); ++j){
          if(mydata.node_subsets[j].is_in_closure(*jt)){
            a_variables.insert(j);
          }
        }

      }
      //for(std::set<int>::iterator kt = a_variables.begin(); kt != a_variables.end(); ++kt){
      //  cind.push_back(mydata.a_start + *kt);
      //  cval.push_back(-1.0);
      //}
      //mydata.fort_borders.push_back(new_border);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0, NULL);
      if(error){
        std::cerr << "ERROR: " << error << " disjoint forts constraint" << endl;
        return 1;
      }
      mydata.num_forts++;

    }

  } //*/



  /*Chosen nodes must have a chosen neighbor*/
  for(int i=0; i<num_nodes; ++i){
    std::vector<int> cind;
    std::vector<double> cval;
    cind.push_back(i);
    cval.push_back(1.0);
    const int node_degree = our_graph.nodes[i].get_degree();
    for(int j=0; j<node_degree; ++j){
      cind.push_back(our_graph.nodes[i].get_adj_list_member(j));
      cval.push_back(-1.0);
    }

    error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 0.0, NULL);
    if(error){
      std::cerr << "ERROR: " << error << " at edge constraints." << endl;
      return 1;
    }

  }//*/

  /*Add the edge constraints*/
  for(int i=0; i<end1.size(); ++i){
    std::vector<int> cind;
    std::vector<double> cval;
    const int end_1 = end1[i];
    const int end_2 = end2[i];

    cind.push_back(num_nodes+i);
    cval.push_back(2.0);
    cind.push_back(end_1);
    cval.push_back(-1.0);
    cind.push_back(end_2);
    cval.push_back(-1.0);

    error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 0.0, NULL);
    if(error){
      std::cerr << "ERROR: " << error << " at edge constraints." << endl;
      return 1;
    }

  }//*/

  /*There must be a spanning tree*/
  {
    std::vector<int> cind;
    std::vector<double> cval;
    for(int i=0; i<end1.size(); ++i){
      cind.push_back(num_nodes+i);
      cval.push_back(-1.0);
    }
    for(int i=0; i<num_nodes; ++i){
      cind.push_back(i);
      cval.push_back(1.0);
    }
    error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 1.0, NULL);
    if(error){
      std::cerr << "ERROR: " << error << " at edge constraints." << endl;
      return 1;
    }
  }//*/


	//End Gurobi Model***************************************************************************
  //error = GRBwrite(zero_forcing, "zero_forcing.lp");
	//Solve the model
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "PreSolve", 0);
  error = GRBsetintparam(GRBgetenv(zero_forcing), "PreCrush", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "Method", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "MIPFocus", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "BranchDir", -1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "MIRCuts", 2);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "GomoryPasses", 5);
  error = GRBsetintparam(GRBgetenv(zero_forcing), "LazyConstraints", 1);
 error = GRBsetintparam(GRBgetenv(zero_forcing), "Threads", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "Symmetry", 2);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "OutputFlag", 0);


  //mydata.orbits = orbits;
  error = GRBsetcallbackfunc(zero_forcing, cut_callback, (void *) &mydata);
  if(error){
	  cout << "ERROR at: GRBsetcallbackfunc" << error << endl;
	  return 1;
	}

  error = GRBupdatemodel(zero_forcing);

  clock_t time_start = clock();
  mydata.start_time = time_start;
      auto chrono_time_start = std::chrono::high_resolution_clock::now();
    mydata.chrono_start_time = chrono_time_start;
  error = GRBoptimize(zero_forcing);
    auto chrono_time_end = std::chrono::high_resolution_clock::now();

	double IPtime;
	error = GRBgetdblattr(zero_forcing, GRB_DBL_ATTR_RUNTIME, &IPtime);
  cout << "Running time is: " << IPtime << endl;

  clock_t time_end = clock();

 std::chrono::duration<double> chrono_time = chrono_time_end - chrono_time_start;
  double total_time = (time_end- time_start)/static_cast<double>( CLOCKS_PER_SEC );

  double optimum;
  error = GRBgetdblattr(zero_forcing, GRB_DBL_ATTR_OBJVAL, &optimum);
  double LP_bound;
  error = GRBgetdblattr(zero_forcing, GRB_DBL_ATTR_OBJBOUND, &LP_bound);

  cout << "time enforcing connectivity was: " << mydata.time_in_connectivity << endl;


//  ofstream fout;
//  string basename = filename;
//	basename = basename.substr(basename.find_last_of('/')+1);
//
//	fout.open("AB_sep_connected_results_chrono_final.txt", ios::app);
//	fout.precision(10);
//	fout << endl << endl;
//	fout << "Results for: " << basename << endl;
//  fout << "ZFS size was: " << optimum << endl;
//  fout << "Lower bound was: " << LP_bound << endl;
//  fout << "Number of forts was: " << mydata.num_forts << endl;
//  fout << "Time finding forts was: " << mydata.time_in_LP << endl;
//  fout << "Number of AB constraints was: " << mydata.counter << endl;
//  fout << "Time enforcing connectivity was: " << mydata.time_in_connectivity << endl;
//  fout << "Gurobi reported time was: " << IPtime << endl;
//  fout << "Self-calculated time was: " << total_time << endl;
//  fout << "Chrono time was: " << chrono_time.count() << endl;
//	fout << endl;
//	fout.close();

  /***Output forcing set information*/

  std::vector<double> edge_value;
  std::vector<double> node_labels;
  std::vector<double> node_labels2;
  cout << "optimal variables are: " << endl;
  cout << "start node variables: " << endl;
  for(int i=0; i<num_nodes;++i){
    double sol_x;
    error = GRBgetdblattrelement(zero_forcing,GRB_DBL_ATTR_X,i,&sol_x);
    cout << "x_" << i  << " = " << sol_x << endl;
    node_labels.push_back(sol_x);
  }
  cout << "subset variables: " << endl;
  for(int i=0; i<num_nodes;++i){
    double sol_x;
    error = GRBgetdblattrelement(zero_forcing,GRB_DBL_ATTR_X,i + num_nodes,&sol_x);
    for(int j=0; j<num_nodes; ++j){
      if(node_subsets[i].is_in_closure(j)){
      node_labels[j] += sol_x;
      }
    }
    node_labels2.push_back(sol_x);
    cout << "x_" << i  << " = " << sol_x << endl;
  }
  for(int i=0; i<num_edges; ++i){
    double sol_x;
    error = GRBgetdblattrelement(zero_forcing,GRB_DBL_ATTR_X,i + edge_variable_start,&sol_x);
    edge_value.push_back(sol_x);
  }


  cout << "time was " << total_time << endl;
      cout << "chrono tine was " << chrono_time.count() << endl;
      cout << "ZFS size was: " << optimum << endl;
  error = GRBfreemodel(zero_forcing);
  if(error){
    cout << "ERROR at: free model" << endl;
    return 1;
  }

  GRBfreeenv(env);

  return 0;
}
