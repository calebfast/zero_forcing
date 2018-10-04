/**********************************
This program is a C++ implementation of the Wavefront Algorithm for finding
minimum zero forcing sets of a graph.

Author: Caleb Fast
Date: 2017

Compile with: g++ -std=c++11 -o wavefront wavefront.cpp
Usage: ./wavefront <path to .edg file>
**********************************/




#include<iostream>
#include<sstream>
#include<fstream>
#include<set>
#include<unordered_map>
#include<string>
#include<stdlib.h>
#include<cstdlib>
#include<ctime>
#include<vector>
#include<chrono>

#define FLOAT_TOL 0.00001

using namespace std;

class node{

public:
  int get_degree(){
    return adj_list.size();
  };
  void add_to_adj_list(int i){
    adj_list.push_back(i);
  };
  int get_adj_list_member(int i){
    return adj_list[i];
  };

private:
  std::vector<int> adj_list;

};


class node_subset{

public:
	std::vector<bool> closure_contains;
  int weight;
  std::set<int> closure;
  std::set<int> subset;

  void set_weight(){
    weight = subset.size();
  };
  int get_weight(){
    return weight;
  };
  bool is_in_closure(int i){
    return closure_contains[i];
  };

  void add_to_subset(int i){
    subset.insert(i);
  };

  void add_to_closure(int i){
    closure.insert(i);
  };

  int num_unfilled_neighbors(node& v, int& unfilled_neighbor){

    const int v_degree = v.get_degree();
    int result = 0;
    for(int i=0; i<v_degree; ++i){
      if(!closure_contains[v.get_adj_list_member(i)]){
	result++;
	unfilled_neighbor = v.get_adj_list_member(i);
      }
    }
    return result;
  };

  int get_closure_size(){
    return closure.size();
  }

  std::set<int>::iterator get_closure_begin(){
    return closure.begin();
  };
  std::set<int>::iterator get_closure_end(){
    return closure.end();
  };
  std::set<int>::iterator get_subset_begin(){
    return subset.begin();
  };
  std::set<int>::iterator get_subset_end(){
    return subset.end();
  };

  void update_closure(std::vector<node>& nodes){

    const int num_nodes = nodes.size();
    closure_contains = std::vector<bool>(num_nodes, false);
    for(std::set<int>::iterator it = subset.begin(); it != subset.end(); ++it){
      closure.insert(*it);
      closure_contains[*it] = true;
    }

    bool end_flag = false;
    while(!end_flag){
      end_flag = true;
      std::set<int> temp_set(closure);
      std::vector<bool> temp_closure_contains(closure_contains);
      for(std::set<int>::iterator it = closure.begin(); it != closure.end(); ++it){
	//Count the number of forced nodes that are adjacent to the node in the closure
	int unforced_sum = 0;
	int white_neighbor;
	const int node_degree = nodes[*it].get_degree();
	for(int i=0; i<node_degree; ++i){
	  const int neighbor = nodes[*it].get_adj_list_member(i);
	  if(!closure_contains[neighbor]){
	    unforced_sum++;
	    white_neighbor = neighbor;
	  }
	  if(unforced_sum==2){
	    break;
	  }
	}
	if(unforced_sum == 1){

	  temp_set.insert(white_neighbor);

	  temp_closure_contains[white_neighbor] = true;
	  end_flag = false;
	}
      }
      closure.swap(temp_set);
      closure_contains.swap(temp_closure_contains);
    }

  };

  size_t create_hash(){
    std::hash< std::vector<bool> > hasher;
    return hasher(closure_contains);
  };

private:


};


int main(int argc, char* argv[]){


    /******Input Check***************************************************************/
    if(argc != 2){
        cerr << "Error: usage is find_forcing_set <path to .edg file>" << endl;
        return 1;
    }




    /******File Reader*****************************************************************/
    fstream file;
    string filename = argv[1];

    file.open(argv[1]);

    int num_edges;
    int num_nodes;
    std::vector<int> result_subset;

    if(file.is_open()){

    file >> num_nodes;
    file >> num_edges;

    }
    else{
      cerr << "ERROR: could not open file: " << filename << endl;
    }
    node blank_node;
    std::vector<node> nodes(num_nodes, blank_node);

    for(int k=0; k<num_edges;++k){
      int i,j;

      file >> i;
      file >> j;
      nodes[i].add_to_adj_list(j);
      nodes[j].add_to_adj_list(i);
    }
    file.close();




    /******Wavefront (dynamic programming) forcing***************************************/
        //Data Structures
    std::unordered_map<size_t, node_subset > closures;

    int max_degree = -1;
    for(int i=0; i<num_nodes; ++i){
      if(nodes[i].get_degree() > max_degree){
        max_degree = nodes[i].get_degree();
      }
    }
     size_t best_set =0;

    node_subset a_subset;

    a_subset.set_weight();
    a_subset.update_closure(nodes);

    closures.insert(std::make_pair(a_subset.create_hash(),a_subset));

    int stop_criteria = 0;
    std::vector<node_subset>::iterator best_subset_iterator;
    int best_size = num_nodes+5;

    //Start timing
    auto start_time = std::chrono::high_resolution_clock::now();
    //clock_t start_time = clock();

    for(int R= 1; R<num_nodes; ++R){
      cout << "number of closures is " << closures.size() << endl;
      std::unordered_map<size_t, node_subset> temp_map(closures);
      for(std::unordered_map<size_t, node_subset>::iterator it=closures.begin(); it!=closures.end(); ++it){
        for(int j=0; j<num_nodes; ++j){
          int unfilled_neighbor = -1;
          const int k = (it->second).num_unfilled_neighbors(nodes[j], unfilled_neighbor);
	  /*cout << "node is "<< j << endl;
	  cout << "number of unfilled neighbors is " << k << endl;
	  cout << "unfilled neighbor is " << unfilled_neighbor << endl;
	  */

	    //Case for j in closure*************************************************
	  //cout << (it->second).is_in_closure(j) << " " << (it->second).get_weight() << endl;
	    if((it->second).is_in_closure(j) && k-1 <= R - (it->second).get_weight()){
	      /*cout << "old subset is: " << endl;
	      for(std::set<int>::iterator jt = (it->second).get_subset_begin(); jt != (it->second).get_subset_end(); ++jt){
		cout << *jt << ",";
		}*/
	      node_subset new_subset;
	      for(std::set<int>::iterator jt = (it->second).get_subset_begin(); jt != (it->second).get_subset_end(); ++jt){
          new_subset.add_to_subset(*jt);
	      }
	      for(int l=0; l<nodes[j].get_degree(); ++l){
          const int neighbor = nodes[j].get_adj_list_member(l);
          if(neighbor != unfilled_neighbor && !(it->second).is_in_closure(neighbor)){
            new_subset.add_to_subset(neighbor);
          }
	      }
	      /*cout << endl;
	      cout << "new subset is: " << endl;
	      for(std::set<int>::iterator jt = new_subset.get_subset_begin(); jt != new_subset.get_subset_end(); ++jt){
		cout << *jt << ",";
	      }
	      cout << endl;*/
	      new_subset.set_weight();
	      new_subset.update_closure(nodes);

	      //Check if new_subset is already in closures
	      const size_t new_hash = new_subset.create_hash();
	      if(closures.count(new_hash) < 1){
          //closures.insert(std::make_pair(new_hash, new_subset));
          temp_map.insert(std::make_pair(new_hash, new_subset));
	      }
	      /*else{
		if(closures[new_hash].closure_contains != new_subset.closure_contains){
		  cout << "ERROR: Hash_table Collision!!" << endl;
		}*/
	      }
	    //Case for j not in closure*************************************************
	    if(!(it->second).is_in_closure(j) && k <= R - (it->second).get_weight()){
	      node_subset new_subset;
	      for(std::set<int>::iterator jt = (it->second).get_subset_begin(); jt != (it->second).get_subset_end(); ++jt){
          new_subset.add_to_subset(*jt);
	      }
	      new_subset.add_to_subset(j);
	      for(int l=0; l<nodes[j].get_degree(); ++l){
          const int neighbor = nodes[j].get_adj_list_member(l);
          if(neighbor != unfilled_neighbor && !(it->second).is_in_closure(neighbor)){
            new_subset.add_to_subset(neighbor);
          }
	      }
	      new_subset.set_weight();
	      new_subset.update_closure(nodes);

	      //Check if new_subset is already in closures
	      const size_t new_hash = new_subset.create_hash();

	      if(closures.count(new_hash) < 1){
          //closures.insert(std::make_pair(new_hash, new_subset));
          temp_map.insert(std::make_pair(new_hash, new_subset));
	      }
	      /*else{
		if(closures[new_hash].closure_contains != new_subset.closure_contains){
		  cout << "ERROR: Hash_table Collision!!" << endl;
		}*/
	    }
        }
        if(it->second.get_weight() < R - max_degree){
          temp_map.erase(it->first);
        }
      }
      closures.swap(temp_map);

      //Check if one of closures is entire graph
      for(std::unordered_map<size_t, node_subset>::iterator it = closures.begin(); it != closures.end(); ++it){

        if((it->second).get_closure_size() == num_nodes && (it->second).get_weight() < best_size){
          stop_criteria = 1;
          best_size = (it->second).get_weight();
          best_set = it->first;
        }
      }

      if(stop_criteria == 1){
	break;
      }
    }

  auto end_time = std::chrono::high_resolution_clock::now();
    //clock_t end_time = clock();
    std::chrono::duration<double> total_time = end_time - start_time;
    //double total_time = (end_time - start_time)/static_cast<double>(CLOCKS_PER_SEC);

    cout << "Clock time is: " << total_time.count() << endl;
    cout << "The best forcing set has size: " << closures[best_set].get_weight() << endl;
  for(std::set<int>::iterator it = closures[best_set].get_subset_begin(); it != closures[best_set].get_subset_end(); ++it){
    cout << *it << endl;

  }

  ofstream fout;
  string basename = filename;
	basename = basename.substr(basename.find_last_of('/')+1);

	fout.open("wavefront_results_chrono.txt", ios::app);
	fout.precision(10);
	fout << endl << endl;
	fout << "Results for: " << basename << endl;
  fout << "ZFS size was: " << closures[best_set].get_weight() << endl;
  fout << "Time was: " << total_time.count() << endl;
  fout << "Number of closures was " << closures.size() << endl;
	fout << endl;
	fout.close();



    return 0;
}
