
/*********************************************************************************
This program uses a branch and bound style method based on generating connected
subsets to solve the connected zero-forcing problem

Author: Caleb Fast

Compilation: g++ -std=c++11 -o b_and_b_connected_forcing b_and_b_connected_forcing.cpp
Usage: ./b_and_b_connected_forcing <path to .edg file>

**********************************************************************************/


#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<set>
#include<cfloat>
#include<string>
#include<cstring>
#include<ctime>
#include<cmath>
#include<stdlib.h>
#include<algorithm>
#include<chrono>

//Determines how long to run algorithm before giving up and stopping
#define STOPPING_TIME 7200

using namespace std;

namespace forcing{
class node{

public:
  int child1, child2, child3, daddy, num_children, mark = 0;  //Stuff for trees
  std::set<int> subtree;

  int get_degree(){
    return adj_list.size();
  };
  void add_to_adj_list(int i){
    adj_list.push_back(i);
  };
  int get_adj_list_member(int i){
    return adj_list[i];
  };

  void swap_adj_lists(std::vector<int> new_adj_list){
    adj_list.swap(new_adj_list);
  }


private:
  std::vector<int> adj_list;

};

class node_subset{

public:

  node_subset(){};
  ~node_subset(){};
  node_subset(const node_subset& B){
    weight = B.weight;
    closure_contains = B.closure_contains;
    //subset_contains = B.subset_contains;
    closure = B.closure;
    subset = B.subset;
  };
  node_subset operator=(const node_subset& B){
    node_subset A;
    A.weight = B.weight;
    A.closure_contains = B.closure_contains;
    //A.subset_contains = B.subset_contains;
    A.closure = B.closure;
    A.subset = B.subset;
    return A;
  };

  void clear_subset(){
    weight = 0;
    closure_contains.clear();
    closure.clear();
    subset.clear();
  };

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
	unfilled_neighbor = v.get_adj_list_member(i);
	result++;
      }
    }
    return result;
  };

  int get_closure_size(){
    return closure.size();
  };

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
    //subset_contains = std::vector<bool>(num_nodes, false);
    for(std::set<int>::iterator it = subset.begin(); it != subset.end(); ++it){
      closure.insert(*it);
      closure_contains[*it] = true;
      //subset_contains[*it] = true;
    }

    bool end_flag = false;
    while(!end_flag){
      end_flag = true;
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
	  closure.insert(white_neighbor);
	  closure_contains[white_neighbor]=true;
	  end_flag = false;
	}
      }
    }

  };



  void update_closure_get_ends(std::vector<node>& nodes, std::set<int>& end_nodes){

    const int num_nodes = nodes.size();

    closure_contains = std::vector<bool>(num_nodes, false);
    //subset_contains = std::vector<bool>(num_nodes, false);
    for(std::set<int>::iterator it = subset.begin(); it != subset.end(); ++it){
      closure.insert(*it);
      end_nodes.insert(*it);
      closure_contains[*it] = true;
      //subset_contains[*it] = true;
    }

    bool end_flag = false;
    while(!end_flag){
      end_flag = true;
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
          closure.insert(white_neighbor);
          end_nodes.insert(white_neighbor);
          end_nodes.erase(*it);
          closure_contains[white_neighbor]=true;
          end_flag = false;
        }
      }
    }

    //cout << "subset size is " << subset.size() << endl;
    //cout << "end_nodes size is " << end_nodes.size() << endl;

  };



  size_t create_hash(){
    hash< vector<bool> > hasher;
    return hasher(closure_contains);
  };


private:
  int weight;
  std::vector<bool> closure_contains;
  //std::vector<bool> subset_contains;
  std::set<int> closure;
  std::set<int> subset;

};


class edge_list{

  public:
    int get_edge_number(int end1, int end2){
      for(int i=0; i<edges.size(); ++i){
        if(edges[i].first == end1 && edges[i].second == end2){
          return i;
        }

      }
      return -1;
    };

    int get_edge_number_v2(int end1, int end2){
      for(int i=0; i<edges.size(); ++i){
        if(edges[i].first == end1 && edges[i].second == end2){
          return i;
        }
        if(edges[i].first == end2 && edges[i].second == end1){
          return i;
        }
      }
      return -1;
    };


    //Returns 1 if edge is end1--end2.  Returns -1 if edge is end2--end1
    int get_edge_orientation(int end1, int end2){
      for(int i=0; i<edges.size(); ++i){
        if(edges[i].first == end1 && edges[i].second == end2){
          return 1;
        }
        if(edges[i].first == end2 && edges[i].second == end1){
          return -1;
        }
      }
    }

    int get_size(){
      return edges.size();
    };
    int get_end1(int i){
      return edges[i].first;
    };
    int get_end2(int i){
      return edges[i].second;
    };

    void add_edge(int end1, int end2){
      std::pair<int, int> p (end1, end2);
      //p.first = end1;
      //p.second = end2;
      edges.push_back(p);
    };

    void swap_edge_list(std::vector<std::pair<int,int> >& new_edges){
      edges.swap(new_edges);
    }

  private:
    std::vector<std::pair<int,int> > edges;

};

class graph{
public:
  std::vector<node> nodes;
  edge_list edges;

  graph operator=(graph& rhs_graph){
    graph new_graph;
    new_graph.nodes = rhs_graph.nodes;
    const int edge_list_size = rhs_graph.edges.get_size();
    for(int i=0; i<edge_list_size; ++i){
      new_graph.edges.add_edge(rhs_graph.edges.get_end1(i), rhs_graph.edges.get_end2(i));
    }
    return new_graph;
  };

    int generate_connected_subgraph(std::set<int>& vertices_to_be_considered, std::set<int>& subset_so_far, std::set<int>& neighbors, int& size_limit,
                                  std::set<int>& zero_forcing_set, clock_t start_time, double time_limit){
    double runtime = (clock() - start_time)/static_cast<double>(CLOCKS_PER_SEC);
    if(runtime > time_limit){
      return -1;
    }
    //cout << "recursing" << endl;
    int candidate;
    bool yield_set = true;
    //look for a new candidate to add to the subgraph
    if(subset_so_far.size() == 0){
      if(vertices_to_be_considered.size() >= 1){
        candidate = *vertices_to_be_considered.begin();
        //cout << "candidate is " << candidate << endl;
        yield_set = false;
      }
    }
    else{
      for(std::set<int>::iterator it = neighbors.begin(); it != neighbors.end(); ++it){
        //cout << "Looking for new candidate" << endl;
        if(vertices_to_be_considered.count(*it) == 1){
          //cout << "found new candidate" << endl;
          candidate = *it;
          yield_set = false;
          break;
        }
      }
    }
    if(yield_set){
      //cout << "yielding set" << endl;
      //Check subset to see if it is zero-forcing
      node_subset new_subset;
      for(std::set<int>::iterator it = subset_so_far.begin(); it != subset_so_far.end(); ++it){
        new_subset.add_to_subset(*it);
      }
      new_subset.set_weight();
      new_subset.update_closure(nodes);
      if(new_subset.get_closure_size() == nodes.size()){
        //Subset is zero_forcing
        const int weight = new_subset.get_weight();
        if(weight < size_limit){
          //cout << "found new incumbent with weight " << weight << endl;
          size_limit = weight;
          zero_forcing_set = subset_so_far;
        }
      }
    }
    else{
      //Recursion
      std::set<int> temp1_vertices_to_be_considered(vertices_to_be_considered);
      temp1_vertices_to_be_considered.erase(candidate);
      std::set<int> temp1_subset_so_far(subset_so_far);
      std::set<int> temp1_neighbors(neighbors);
      generate_connected_subgraph(temp1_vertices_to_be_considered, temp1_subset_so_far, temp1_neighbors, size_limit, zero_forcing_set, start_time, time_limit);
      //check to make sure we wont go over the size limit
      if(subset_so_far.size() < size_limit-1){
        std::set<int> temp2_vertices_to_be_considered(vertices_to_be_considered);
        temp2_vertices_to_be_considered.erase(candidate);
        std::set<int> temp2_subset_so_far(subset_so_far);
        temp2_subset_so_far.insert(candidate);
        std::set<int> temp2_neighbors(neighbors);
        const int candidate_degree = nodes[candidate].get_degree();
        for(int i=0; i<candidate_degree; ++i){
          temp2_neighbors.insert(nodes[candidate].get_adj_list_member(i));
        }
        //Recursion
        generate_connected_subgraph(temp2_vertices_to_be_considered, temp2_subset_so_far, temp2_neighbors, size_limit, zero_forcing_set, start_time, time_limit);
      }

    }
    return size_limit;
  };

  int b_and_b_connected_zero_forcing(std::set<int>& zero_forcing_set, const int time_limit){
     std::set<int> vertices_to_be_considered;
     std::set<int> subset_so_far;
     std::set<int> neighbors;
     for(int i=0; i<nodes.size(); ++i){
      vertices_to_be_considered.insert(i);
     }
     int size_limit = nodes.size() + 1;
     clock_t start_time = clock();
     //cout << "starting recursion" << endl;
     int return_value = generate_connected_subgraph(vertices_to_be_considered, subset_so_far, neighbors, size_limit, zero_forcing_set, start_time, time_limit);
    cout << "return value is " << return_value << endl;
     if(return_value == -1){
      return -1;
     }

     return size_limit;

  };


};

}

int main(int argc, char* argv[]){


    if(argc != 2){
        cerr << "Error: usage is find_forcing_set <path to .edg file>" << endl;
        return 1;
    }

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


    if(file.is_open()){

    file >> num_nodes;
    file >> num_edges;
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
      nodes[i].add_to_adj_list(j);
      nodes[j].add_to_adj_list(i);

      our_graph.edges.add_edge(i,j);
      our_graph.edges.add_edge(j,i);
    }
    file.close();

    our_graph.nodes = nodes;
    num_edges = our_graph.edges.get_size();

  std::set<int> zero_forcing_set;
  cout << "starting." << endl;
  clock_t time_start = clock();
   auto chrono_time_start = std::chrono::high_resolution_clock::now();
  int return_value = our_graph.b_and_b_connected_zero_forcing(zero_forcing_set, STOPPING_TIME);
     auto chrono_time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> chrono_time = chrono_time_end - chrono_time_start;
  clock_t time_end = clock();
  cout << "finished." << endl;
  double total_time = (time_end- time_start)/static_cast<double>( CLOCKS_PER_SEC );

  cout << "Forcing set size is " << return_value << endl;
  cout << "Set size is " << zero_forcing_set.size() << endl;
  cout << "Running time is " << total_time << endl;
  cout << "The zero forcing set is ";
  for(std::set<int>::iterator it = zero_forcing_set.begin(); it != zero_forcing_set.end(); ++it){
    cout << *it << ", ";
  }
  cout << endl;
//  ofstream fout;
//  string basename = filename;
//	basename = basename.substr(basename.find_last_of('/')+1);
//
//	fout.open("brute_connected_results.txt", ios::app);
//	fout.precision(10);
//	fout << endl << endl;
//	fout << "Results for: " << basename << endl;
//	fout << "ZFS size was: " << zero_forcing_set.size() << endl;
//	fout << "Best set found had size: " << return_value << endl;
//  fout << "Self-calculated time was: " << total_time << endl;
//  fout << "chrono time was: " << chrono_time.count() << endl;
//	fout << endl;
//	fout.close();



  return 0;
}
