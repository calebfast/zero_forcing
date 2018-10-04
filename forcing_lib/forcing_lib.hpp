#ifndef FORCING_LIB_H
#define FORCING_LIB_H

#include<vector>
#include<set>
#include<map>
#include<stack>
#include<iostream>
#include<unordered_map>
#include<stdlib.h>
#include<queue>
#include <ctime>
#include <chrono>
//#include"../nauty26r5/nausparse.h"


/**INCLUDE PATH TO gurobi_c++.h**/
//#include "<path to gurobi_c++.h>"
//Examples
//For my ASUS Laptop
//#include "/opt/gurobi651/linux64/include/gurobi_c++.h"
//For Harvey
#include "/opt1/opt.LOCAL//gurobi752/linux64/include/gurobi_c++.h"

#define FLOAT_TOL 0.00001
#define MIN_CLOSURE_SIZE 1
#define CUTOFF 0.9
#define MIN_CUT_VIOLATION 0.25
#define MAX_FORT_SIZE 15
#define MAX_FORT_NUMBER 100
#define SPARSIFY_PENALTY 0.0001


using namespace std;


/**Classes**********************************************/
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

  void remove_from_subset(int i){
    subset.erase(i);
    closure.clear();
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


//void identify_nodes(int node1, int node2, const int num_nodes, std::vector<double>& adj_matrix, vector<double>& new_matrix, std::map< int, std::set<int> >& cut_map);
//void legal_order(const int num_nodes, std::vector<double>& adj_matrix, std::vector<int>& ordering);
//identifies the two nodes node1 and node2
void identify_nodes(int node1, int node2, const int num_nodes, std::vector<double>& adj_matrix, vector<double>& new_matrix, std::map< int, std::set<int> >& cut_map){

  if(node2 < node1){
    std::swap(node1, node2);
  }

  //update the map
  std::copy(cut_map[node2].begin(), cut_map[node2].end(), std::inserter(cut_map[node1], cut_map[node1].end()));
  //cut_map.erase(cut_map.find(node2));

  for(int i=node2; i<num_nodes-1; ++i){
    cut_map[i] = cut_map[i+1];
  }

  cut_map.erase(num_nodes-1);
  //Add the rows and columns
  /*for(int i=0; i < num_nodes; ++i){
    if(i != node1 && i != node2){
      adj_matrix[i*num_nodes + node1] += adj_matrix[i*num_nodes+node2];
      adj_matrix[node1*num_nodes + i] += adj_matrix[node2*num_nodes+i];
    }
    }*/

  //Reset the node1xnode1 entry
  //adj_matrix[node1*num_nodes + node1] = 0;

  //construct the new smaller adjacency vector
  for(int i=0; i < num_nodes; ++i){
    for(int j=0; j < num_nodes; ++j){
      if(j != node2 && i != node2){
	if(i != node1 && j != node1){
	  new_matrix.push_back(adj_matrix[i*num_nodes+j]);
	}
	if(j == node1 && i != node1){
	  //cout << "pushing back: i is: " << i << " j is " << j << endl;
	  new_matrix.push_back(adj_matrix[i*num_nodes+j] + adj_matrix[i*num_nodes +node2]);
	}
	if(i == node1 && j != node1){
	  //cout << "pushing back: i is: " << i << " j is " << j << endl;
	  new_matrix.push_back(adj_matrix[i*num_nodes+j] + adj_matrix[node2*num_nodes+j]);
	}
	if(i == node1 && j == node1){
	  //cout << "pushing back: i is: " << i << " j is " << j << endl;
	  new_matrix.push_back(0);
	}
      }
    }
  }

  new_matrix[node1*(num_nodes-1)+node1] = 0;

}

//Finds a legal ordering of the graph represented in adj_matrix and returns it in ordering
void legal_order(const int num_nodes, std::vector<double>& adj_matrix, std::vector<int>& ordering){

  int next_node;
  std::vector<int> current_set;
  ordering.push_back(0);
  current_set.push_back(0);
  for(int k=0; k < num_nodes-1; k++){
    float max_cut_size = -1;
    const int ordering_size = ordering.size();
    for(int i = 0; i < num_nodes; ++i){
      if(find(current_set.begin(), current_set.end(), i) == current_set.end()){
	float current_cut_size = 0;
	for(int j=0; j<ordering_size; ++j){
	  current_cut_size += adj_matrix[ordering[j]*num_nodes+i];
	}
	if(current_cut_size > max_cut_size){
	  max_cut_size = current_cut_size;
	  next_node = i;
	  //cout << "next_node is: " << next_node << endl;
	}
      }
    }
    ordering.push_back(next_node);
    current_set.push_back(next_node);
  }
  /*
  cout << "Ordering is: ";
  for(std::vector<int>::iterator it=ordering.begin(); it!=ordering.end(); ++it){
    cout << *it;
  }
  cout << endl;
  */
}

//void node_identification_min_cut(const int num_nodes, std::vector<double>& adj_matrix, std::set<int>& min_cut, double& min_cut_value, std::set< std::set<int> >& cuts_under_one);
//Finds all cuts in the graph that are less that two using the node identification min cut algorithm.
void node_identification_min_cut(const int num_nodes, std::vector<double>& adj_matrix, std::set<int>& min_cut, double& min_cut_value, std::set< std::set<int> >& cuts_under_one){

  min_cut_value = FLT_MAX;
  int current_num_nodes = num_nodes;
  std::vector<double> graph_vector (adj_matrix);
  std::vector<double> temp_vector;
  std::map< int, std::set<int> > cut_map;
  for(int i=0; i<num_nodes; ++i){
    std::set<int> dummy_set;
    dummy_set.insert(i);
    cut_map.insert( std::pair<int, std::set<int> >(i, dummy_set) );
  }
  while(current_num_nodes > 1){

    //Find a legal ordering
    std::vector<int> ordering;
    legal_order(current_num_nodes, graph_vector, ordering);
    float current_cut_value = 0;
    //cout << "current n is " << ordering[current_num_nodes-1] << endl;
    //cout << "current n-1 is " << ordering[current_num_nodes-2] << endl;
    //print_matrix_vector(current_num_nodes, graph_vector);
    //Find the value of the current cut
    const int row = current_num_nodes*ordering[current_num_nodes-1];
    for(int i=0; i < current_num_nodes; ++i){
      current_cut_value += graph_vector[row+i];
    }

    if(current_cut_value < min_cut_value){
      min_cut_value = current_cut_value;
      min_cut = cut_map[ordering[current_num_nodes-1]];
    }

    //Store the cuts that have weight < 2
    if(current_cut_value < 1.0){
      //std::set<int> dummy_set (
      /*cout << "inserting dummy_set = " << endl;
      for(std::set<int>::iterator it = dummy_set.begin(); it != dummy_set.end(); ++it){
	cout << *it;
      }
      cout << endl;
      */
      cuts_under_one.insert(cut_map[ordering[current_num_nodes-1]]);
    }

    identify_nodes(ordering[current_num_nodes-1], ordering[current_num_nodes-2], current_num_nodes, graph_vector, temp_vector, cut_map);
    graph_vector = temp_vector;
    temp_vector.clear();
    current_num_nodes--;
  }

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


  int make_fort_minimal(std::set<int>& fort){
    std::map<int,int> node_degrees_in_fort;
    for(int i=0; i<fort.size(); ++i){
      const int i_degree = nodes[i].get_degree();
      for(int j=0; j<i_degree; ++j){
        const int neighbor = nodes[i].get_adj_list_member(j);
        std::pair<std::map<int,int>::iterator, bool> return_pair;
        return_pair = node_degrees_in_fort.insert(std::make_pair(neighbor, 1));
        if(!return_pair.second){
          node_degrees_in_fort[neighbor] += 1;
        }
      }
    }

    for(int i=0; i<fort.size(); ++i){
      bool can_be_deleted = true;
      const int i_degree = nodes[i].get_degree();
      for(int j=0; j<i_degree; ++j){
        const int neighbor = nodes[i].get_adj_list_member(j);
        if(node_degrees_in_fort[neighbor] < 3 && fort.count(neighbor) != 1){
          can_be_deleted = false;
        }
      }
      if(can_be_deleted){
        for(int j=0; j<i_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          node_degrees_in_fort[neighbor] -= 1;
        }
        fort.erase(i);
      }
    }

    return fort.size();
  };



  int get_single_component(std::set<int>& component, std::set<int>& candidates){


    if(candidates.size() == 0){
      return 0;
    }

    std::queue<int> bfs_queue;
    bfs_queue.push(*candidates.begin());
    candidates.erase(candidates.begin());

    while(!bfs_queue.empty()){
      const int current_node = bfs_queue.front();
      component.insert(current_node);
      const int node_degree = nodes[current_node].get_degree();
      for(int i=0; i<node_degree; ++i){
        const int neighbor = nodes[current_node].get_adj_list_member(i);
        if(candidates.count(neighbor) == 1){
          bfs_queue.push(neighbor);
          candidates.erase(neighbor);
        }
      }
      bfs_queue.pop();
    }

    return 1;

  };

  int get_connected_components(std::set<int>& separator, std::vector<std::set<int> >& components){

    const int num_nodes = nodes.size();
    std::set<int> candidates;
    for(int i=0; i<num_nodes; ++i){
      if(separator.count(i) < 1){
        candidates.insert(i);
      }
    }

    while(candidates.size() > 0){
      std::set<int> new_component;
      get_single_component(new_component, candidates);
      components.push_back(new_component);
    }

    return components.size();

  };

  int get_minimal_separator(std::set<int>& current_separator, std::vector<std::set<int> >& components, const int a, const int b){
      int num_components = get_connected_components(current_separator, components);

      //Find out which components contain a and b
      int component_with_a;
      int component_with_b;
      bool found_a = false;
      bool found_b = false;
      for(int i=0; i<num_components; ++i){
        if(!found_a && components[i].count(a) == 1){
          component_with_a = i;
        }
        if(!found_b && components[i].count(b) == 1){
          component_with_b = i;
        }
        if(found_a && found_b){
          break;
        }
      }


      //Check which nodes really need to be in separator
      std::vector<int> removed_nodes;
      for(std::set<int>::iterator it = current_separator.begin(); it != current_separator.end(); ++it){
        num_components = components.size();
        std::vector< std::set<int> > new_components;
        std::set<int> merged_component;
        bool merge_components = true;
        const int current_node = *it;
        const int node_degree = nodes[current_node].get_degree();

        //Check if node is adjacent to both a and b components
        bool adjacent_to_a = false;
        bool adjacent_to_b = false;
        for(int k=0; k<node_degree; ++k){
          if(components[component_with_a].count(nodes[current_node].get_adj_list_member(k))==1){
           adjacent_to_a = true;
          }
          if(components[component_with_b].count(nodes[current_node].get_adj_list_member(k))==1){
           adjacent_to_b = true;
          }
        }
        if(adjacent_to_a && adjacent_to_b){
          merge_components = false;
        }
        bool a_in_merged = false;
        bool b_in_merged = false;
        if(merge_components){
          int counter = 0;
          int loop_counter = 0;
          merged_component.insert(current_node);
          std::vector<std::set<int> > adjacent_components;
          //Check which components are adjacent to current node
          for(int j=0; j<num_components; ++j){

            bool adjacent = false;
            for(int k=0; k<node_degree; ++k){
              //check if neighbor is in component
              if(components[j].count(nodes[current_node].get_adj_list_member(k)) ==1){
                merged_component.insert(components[j].begin(), components[j].end());
                adjacent = true;
                if(j == component_with_a){
                  a_in_merged = true;
                }
                if(j == component_with_b){
                  b_in_merged = true;
                }
                break;
              }
            }
            if(!adjacent){
              new_components.push_back(components[j]);
              if(j == component_with_a){
                component_with_a = counter;
              }
              if(j == component_with_b){
                component_with_b = counter;
              }
              counter++;
            }
          }

        //Merge the components
          if(a_in_merged){
            component_with_a = counter;
          }
          if(b_in_merged){
            component_with_b = counter;
          }
          removed_nodes.push_back(current_node);
          new_components.push_back(merged_component);
          components.swap(new_components);
        }
      }

      for(std::vector<int>::iterator it = removed_nodes.begin(); it != removed_nodes.end(); ++it){
        current_separator.erase(*it);
      }

      return current_separator.size();
  };

  float get_weighted_cut_to_subset(const int node, std::vector<int>& subset, std::vector<float>& edge_weights){
    const int node_degree = nodes[node].get_degree();
    float result =0;
    for(int i=0; i<node_degree; ++i){
      if(std::find(subset.begin(), subset.end(), nodes[node].get_adj_list_member(i)) != subset.end()){
        result+= edge_weights[edges.get_edge_number_v2(node, nodes[node].get_adj_list_member(i))];
      }
    }
    return result;
  };



  //For find_fort_greedy.  Returns -1 if all neighbors are satisfied
  std::pair<int, int> get_unsatisfied_neighbor(node_subset& forced_nodes, std::set<int>& C, std::set<int>& F, std::set<int>& F_done, const int current_node){

    std::pair<int, int> result(-1,-1);
    const int num_neighbors = nodes[current_node].get_degree();
    for(int i=0; i<num_neighbors; ++i){
      const int neighbor = nodes[current_node].get_adj_list_member(i);
      if(F.count(neighbor)==1 || F_done.count(neighbor)==1){
        continue;
      }
      const int num_2nd_neighbors = nodes[neighbor].get_degree();
      int num_in_fort = 0;
      for(int j=0; j<num_2nd_neighbors; ++j){
        const int temp = nodes[neighbor].get_adj_list_member(j);
        if(F.count(temp)==1 || F_done.count(temp)==1){
          num_in_fort++;
        }
        else if(C.count(temp) == 1 && !forced_nodes.is_in_closure(temp)){
          result.second = temp;
        }
      }
      if(num_in_fort < 2){
        result.first = neighbor;
      }
    }

    return result;
  };

  //Finds a fort using a greedy algorithm of given subset if subset is not a forcing set
  int find_fort_greedy( std::set<std::set<int> >& forts, std::vector<double>& weights){

    const int num_nodes = nodes.size();
    node_subset forced_nodes;
    for(int i=0; i<num_nodes; ++i){
      if(weights[i] >= 1 - FLOAT_TOL){
        forced_nodes.add_to_subset(i);
      }
    }

    forced_nodes.update_closure(nodes);

    std::set<int> outside_nodes;
    if(forced_nodes.get_closure_size() < num_nodes){
      for(int i=0; i<num_nodes; ++i){
        if(!forced_nodes.is_in_closure(i)){
          outside_nodes.insert(i);
        }
      }
    }

    for(std::set<int>::iterator it = outside_nodes.begin(); it != outside_nodes.end(); ++it){
      std::set<int> C(outside_nodes);
      std::set<int> F;
      std::set<int> F_done;
      F.insert(*it);
      C.erase(*it);
      while(!F.empty()){
        //Choose and node in F
        const int current_node = *(F.begin());
        /*cout << "F is" << endl;
        for(std::set<int>::iterator kt = F.begin(); kt != F.end(); ++kt){
          cout << *kt << ",";
        }
        cout << endl;
        cout << "C is " << endl;
        for(std::set<int>::iterator kt = C.begin(); kt != C.end(); ++kt){
          cout << *kt << ",";
        }
        cout << endl;
        cout << "current_node is " << current_node << endl;
        */
        //Test whether node has a neighbor not in FUF_done that only has one neighbor in FUF_done
        std::pair<int,int> unsatisfied;
        unsatisfied = get_unsatisfied_neighbor(forced_nodes, C, F, F_done, current_node);
        //cout << "unsatisfied is " << unsatisfied.first << " " << unsatisfied.second << endl;

        if(unsatisfied.first == -1){
          F.erase(F.begin());
          F_done.insert(current_node);
        }
        else{
          if(unsatisfied.second != -1){
            F.insert(unsatisfied.second);
            C.erase(unsatisfied.second);
          }
          else{
            F.insert(unsatisfied.first);
            C.erase(unsatisfied.first);
          }

        }
      }

      forts.insert(F_done);
    }

    return 0;
  };


  //For find_fort_greedy_style_2
  int get_node_degree_in_fort(std::set<int>& fort, const int node){
    const int degree = nodes[node].get_degree();
    int num_neighbors_in_fort = 0;
    for(int i=0; i<degree; ++i){
      if(fort.count(nodes[node].get_adj_list_member(i)) == 1){
        num_neighbors_in_fort++;
      }
    }
    return num_neighbors_in_fort;
  };

  void find_fort_greedy_style_2(std::set<std::set<int> >& forts, std::vector<double>& weights){

    const int num_nodes = nodes.size();
    node_subset forced_nodes;
    for(int i=0; i<num_nodes; ++i){
      if(weights[i] >= 1 - FLOAT_TOL){
        forced_nodes.add_to_subset(i);
      }
    }

    forced_nodes.update_closure(nodes);

    std::set<int> outside_nodes;
    if(forced_nodes.get_closure_size() < num_nodes){
      for(int i=0; i<num_nodes; ++i){
        if(!forced_nodes.is_in_closure(i)){
          outside_nodes.insert(i);
        }
      }
    }


    std::stack<int> candidates;
    for(std::set<int>::iterator it = outside_nodes.begin(); it != outside_nodes.end(); ++it){
      candidates.push(*it);
    }

    while(!candidates.empty()){
        const int current_node = candidates.top();
        const int current_node_degree = nodes[current_node].get_degree();
        bool can_delete = true;
        for(int i=0; i<current_node_degree; ++i){
          const int neighbor = nodes[current_node].get_adj_list_member(i);
          if(outside_nodes.count(neighbor) !=1){
            int num_neighbors_in_fort = get_node_degree_in_fort(outside_nodes, neighbor);
            if(num_neighbors_in_fort == 2){
              can_delete = false;
              candidates.pop();
              break;
            }
          }
        }
        int num_neighbors_in_fort = get_node_degree_in_fort(outside_nodes, current_node);
        //cout << "num_neighbors_in_fort is " << num_neighbors_in_fort << endl;
        if(can_delete && num_neighbors_in_fort == 1){
          candidates.pop();
        }
        else if(can_delete){
          outside_nodes.erase(current_node);
          candidates.pop();
        }
    }

    if(outside_nodes.size() > 0){
      forts.insert(outside_nodes);
    }
    /*cout << "fort is " << endl;
    for(std::set<int>::iterator it = outside_nodes.begin(); it != outside_nodes.end(); ++it){
      cout << *it << ",";
    }
    cout << endl;*/
    return;

  };


  void find_fort_max_complement(std::set<std::set<int> >& forts, std::set<int>& forbidden_nodes){

    const int num_nodes = nodes.size();
    node_subset forced_nodes;
    //cout << "Forced nodes are: ";
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it != forbidden_nodes.end(); ++it){
      forced_nodes.add_to_subset(*it);
      //cout << *it << ", ";
    }
    //cout << endl;

    forced_nodes.update_closure(nodes);
    bool is_forcing = true;
    std::set<int> outside_nodes;
    //cout << "num_nodes is " << num_nodes << endl;
    //cout << "closure size is " << forced_nodes.get_closure_size() << endl;
    if(forced_nodes.get_closure_size() < num_nodes){
      is_forcing = false;
      for(int i=0; i<num_nodes; ++i){
        if(!forced_nodes.is_in_closure(i)){
          outside_nodes.insert(i);
        }
      }
    }

    
    std::stack<int> candidates;
    for(std::set<int>::iterator it = outside_nodes.begin(); it != outside_nodes.end(); ++it){
      candidates.push(*it);
    }
    
    while(!candidates.empty()){
        const int current_node = candidates.top();
        forced_nodes.add_to_subset(current_node);
        candidates.pop();
        forced_nodes.update_closure(nodes);
        if(forced_nodes.get_closure_size() == num_nodes){
          forced_nodes.remove_from_subset(current_node);
          forced_nodes.update_closure(nodes);
	}
        else{
          while(!candidates.empty() && forced_nodes.is_in_closure(candidates.top())){
	      candidates.pop();
	  }
	}
    }
	outside_nodes.clear();
      for(int i=0; i<num_nodes; ++i){
        if(!forced_nodes.is_in_closure(i)){
          outside_nodes.insert(i);
        }
      }
    if(!is_forcing){
      forts.insert(outside_nodes);
    }
    /*cout << "is_forciong = " << is_forcing << endl;
    // cout << "fort is " << endl;
    for(std::set<int>::iterator it = outside_nodes.begin(); it != outside_nodes.end(); ++it){
      cout << *it << ",";
    }
    cout << endl;
    for(std::set<std::set<int> >::iterator st = forts.begin(); st != forts.end(); ++st){
      for(std::set<int>::iterator it = st->begin(); it != st->end(); ++it){
        cout << *it << ",";
      }
    }
    cout << endl;*/
    return;

  };



  //Finds a fort using closure of given subset if subset is not a forcing set
  int find_fort_closure(std::set<std::set<int> >& forts, std::vector<double>& weights){

    const int num_nodes = nodes.size();
    node_subset forced_nodes;
    for(int i=0; i<num_nodes; ++i){
      if(weights[i] >= 1 - FLOAT_TOL){
        forced_nodes.add_to_subset(i);
      }
    }

    forced_nodes.update_closure(nodes);
    std::set<int> new_fort;
    if(forced_nodes.get_closure_size() < num_nodes){
      for(int i=0; i<num_nodes; ++i){
        if(!forced_nodes.is_in_closure(i)){
          new_fort.insert(i);
        }
      }
    }

    if(new_fort.size() > 0){
      make_fort_minimal(new_fort);
      forts.insert(new_fort);
    }

    return 0;
  };


  //Finds intersection of fort constraints that gives better constraint
  int find_fort_intersection_cuts(GRBenv* env, const int num_nodes, const int desired_number_forts_in_cut, std::vector<std::set<int> >& forts, std::vector<int>& forbidden_fort_indices, std::vector<int>& fort_indices, std::vector<double>& weights){
    //cout << "finding cut intersections " << endl;
    GRBmodel *model;
    int error;
    const int num_variables = forts.size();
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));
    for(int i=0; i<num_variables; ++i){
      types[i] = GRB_BINARY;
      objective_coeffs[i] = 0.0;
      for(int j=0; j<num_nodes; ++j){
        if(forts[i].count(j) ==1){
          objective_coeffs[i] += weights[j]/((double) desired_number_forts_in_cut - 1);
        }
      }
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    /*for(int i=0; i<forbidden_fort_indices.size(); ++i){
      upper_bounds[forbidden_fort_indices[i]] = 0.0;
    }*/
    //cout << "creating model" << endl;
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

    //Make sure number of forts is the same as desired
    {
    std::vector<int> cind;
    std::vector<double> cval;
    for(int i=0; i<num_variables; ++i){
      cind.push_back(i);
      cval.push_back(1.0);
    }
    //cout << "adding first constraint" << endl;
    error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, (double) desired_number_forts_in_cut, NULL);
    if(error){
      cout << "ERROR in intersection cuts at: add constraint 1: " << error << endl;
      return 1;
    }
    }

    //Add constraints that force intesection of chosen forts to be empty
    //cout << "done with first constraint" << endl;
    for(int i=0; i<num_nodes; ++i){
      std::vector<int> cind;
      std::vector<double> cval;
      for(int j=0; j<num_variables; ++j){
        //Check if node is in fort
        if(forts[j].count(i) == 1){
          cind.push_back(j);
          cval.push_back(1.0);
        }
      }
      //cout << "adding second constraint " << endl;
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, (double) desired_number_forts_in_cut-1.0, NULL);
      if(error){
        cout << "ERROR in intersection cuts at: add constraint 1: " << error << endl;
        return 1;
      }
    }

    //optimize
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: second model update" << endl;
      return 1;
    }

    //cout << "optimizing model " << endl;
    GRBoptimize(model);
    if(error){
      cout << "ERROR at: optimize" << endl;
      return 1;
    }
    //cout << "finished optimizing model " << endl;

    int optimization_status = 0;
    double optimum = -1;
    error = GRBgetintattr(model, "Status",&optimization_status);
    error = GRBgetdblattr(model, "ObjVal", &optimum);
    if(optimization_status == 3 || optimum > 2.0 - MIN_CUT_VIOLATION){
      GRBfreemodel(model);
      return 2;
    }
    //Get solution
    if(optimization_status == 2){
    for(int i=0; i<num_variables; ++i){
      double x_val;
      error = GRBgetdblattrelement(model, "X", i, &x_val);
      if(x_val >= 1- FLOAT_TOL){
        fort_indices.push_back(i);
      }
    }
      GRBfreemodel(model);
      return 0;
    }

  };

  //Finds fort constraints that are likely to give facets of the polytope
  int find_all_forts_LP(GRBenv* env, std::set<std::set<int> >& forts, std::vector<double>& weights, double max_weight, const int cutoff){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = weights[i] + SPARSIFY_PENALTY;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    //Ensure that forbidden nodes are not chosen

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);

      //Require at least 4 nodes left out
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_LESS_EQUAL, num_nodes - 4.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);

          error = GRBupdatemodel(model);
          if(error){
            cout << "ERROR at: third model update" << endl;
            return 1;
          }
        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);

    bool continue_flag = true;
    int counter = 0;
    int loop_counter = 0;
    while(continue_flag && loop_counter < cutoff){
      cout << "loop iteration " << ++counter << endl;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      if(status != GRB_OPTIMAL){
        continue_flag = false;
      }
      if(status == GRB_OPTIMAL){
      double optimum;
      error = GRBgetdblattr(model, "ObjVal", &optimum);
      for(int i=0; i<num_nodes; ++i){
        double x_val;

        error = GRBgetdblattrelement(model, "X", i, &x_val);
        //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
        if(x_val >= 1- FLOAT_TOL){
          new_fort.insert(i);
          //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
          //cout << "optimum is " << optimum << endl;
          //cind.push_back(i);
          //cval.push_back(1.0);
        }
      }

      /*if(new_fort.size() > MAX_FORT_SIZE){
        continue_flag = false;
      }
      cind.push_back(*(new_fort.begin()));
      cval.push_back(1.0);
      */
      if(new_fort.size() >= 1){
        forts.insert(new_fort);
      }


      std::vector<int> cind;
      std::vector<double> cval;

      for(std::set<int>::iterator it = new_fort.begin(); it != new_fort.end(); ++it){
        cind.push_back(*it);
        cval.push_back(1.0);
          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, (double) new_fort.size()-2 , NULL);
        if(error){
          cout << "ERROR adding constraint " << error << endl;
          return 1;
        }
      }

        loop_counter++;
      }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };

  //Finds violated fort constraints
  int find_forts_with_node_LP(GRBenv* env, std::set<std::set<int> >& forts, const int fixed_node, std::set<int>& forbidden_nodes, std::vector<double>& weights, double max_weight){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = weights[i] + SPARSIFY_PENALTY;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    lower_bounds[fixed_node] = 1.0;

    //Ensure that forbidden nodes are not chosen
    if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);

          error = GRBupdatemodel(model);
          if(error){
            cout << "ERROR at: third model update" << endl;
            return 1;
          }
        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);

    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      continue_flag = false;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);
          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            forbidden_nodes.insert(i);
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        for(std::set<int>::iterator it = new_fort.begin(); it != new_fort.end(); ++it){
          if(*it != fixed_node){
            error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
            break;
          }
        }

        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */
          if(optimum < max_weight && new_fort.size() >= 1){
            forts.insert(new_fort);
          }
          else{
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };


  //Finds violated fort constraints
  int find_minimal_border_fort_with_nodes_LP(GRBenv* env, std::set<std::set<int> >& forts, std::set<int>& required_nodes, std::set<int>& forbidden_nodes, std::vector<double>& weights, double max_weight){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = 2*num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = weights[i];
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    for(int i=num_nodes; i<2*num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = 1.0*SPARSIFY_PENALTY;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    for(std::set<int>::iterator it = required_nodes.begin(); it !=required_nodes.end(); ++it){
      lower_bounds[*it] = 1.0;
    }
    //Ensure that forbidden nodes are not chosen
    /*if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }*/

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*constraint for border variables*/
    {
      for(int i=0; i< num_nodes; ++i){
        std::vector<int> cind;
        std::vector<double> cval;
        const int node_degree = nodes[i].get_degree();
        cind.push_back(i);
        cind.push_back(num_nodes + i);
        cval.push_back((double) node_degree);
        cval.push_back(-1.0*node_degree);
        for(int j=0; j<node_degree; ++j){
          cind.push_back(nodes[i].get_adj_list_member(j));
          cval.push_back(-1.0);
        }


      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 0.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }


      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
      }
    }//*/

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);


        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      continue_flag = false;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      //cout << "status is " << status << endl;
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        //cout << "optimum is : " << optimum << endl;
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);
          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            forbidden_nodes.insert(i);
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */
              if(optimum > 1.0){
        double weight_sum = 0;
        for(int i=0; i<num_nodes; ++i){
          weight_sum += weights[i];
        }
        //cout << "sum of weights is " << weight_sum << endl;
      }

          if(optimum < max_weight && new_fort.size() >= 1){
            forts.insert(new_fort);
          }
          else{
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };


  //Finds violated fort constraints
  int find_multiple_minimal_border_fort_LP(GRBenv* env, std::set<std::set<int> >& forts, std::set<int>& forbidden_nodes, std::vector<double>& weights, double max_weight){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = 2*num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = weights[i];
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    for(int i=num_nodes; i<2*num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = 1.0*SPARSIFY_PENALTY;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    //Ensure that forbidden nodes are not chosen
    /*if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }*/

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*constraint for border variables*/
    {
      for(int i=0; i< num_nodes; ++i){
        std::vector<int> cind;
        std::vector<double> cval;
        const int node_degree = nodes[i].get_degree();
        cind.push_back(i);
        cind.push_back(num_nodes + i);
        cval.push_back((double) node_degree);
        cval.push_back(-1.0*node_degree);
        for(int j=0; j<node_degree; ++j){
          cind.push_back(nodes[i].get_adj_list_member(j));
          cval.push_back(-1.0);
        }


      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 0.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }


      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
      }
    }//*/

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);


        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      counter++;
      if(counter > 20){
        continue_flag = false;
      }
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      //cout << "status is " << status << endl;
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        //cout << "optimum is : " << optimum << endl;
        std::vector<double> cval;
        std::vector<int> cind;
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);
          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            cind.push_back(i);
            cval.push_back(1.0);
            //forbidden_nodes.insert(i);
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        error = GRBaddconstr(model,cind.size(),&cind[0],&cval[0],GRB_LESS_EQUAL, (double)cind.size()-1,NULL);
        //error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */
              if(optimum > 1.0){
        double weight_sum = 0;
        for(int i=0; i<num_nodes; ++i){
          weight_sum += weights[i];
        }
        //cout << "sum of weights is " << weight_sum << endl;
      }

          if(optimum < max_weight && new_fort.size() >= 1){
            forts.insert(new_fort);
          }
          else{
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };


  //Finds violated fort constraints
  int find_minimal_border_fort_LP(GRBenv* env, std::set<std::set<int> >& forts, std::set<int>& forbidden_nodes, std::vector<double>& weights, double max_weight, std::vector<int>& num_forts_variable_in){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = 2*num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = 0;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    for(int i=num_nodes; i<2*num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = 1.0;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    //Ensure that forbidden nodes are not chosen
    if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      //cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*constraint for border variables*/
    {
      for(int i=0; i< num_nodes; ++i){
        std::vector<int> cind;
        std::vector<double> cval;
        const int node_degree = nodes[i].get_degree();
        cind.push_back(i);
        cind.push_back(num_nodes + i);
        cval.push_back((double) node_degree);
        cval.push_back(-1.0*node_degree);
        for(int j=0; j<node_degree; ++j){
          cind.push_back(nodes[i].get_adj_list_member(j));
          cval.push_back(-1.0);
        }


      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 0.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }


      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
      }
    }//*/

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);


        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBsetintparam(GRBgetenv(model), "Threads", 1);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      continue_flag = false;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      //cout << "status is " << status << endl;
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        //cout << "optimum is : " << optimum << endl;
        //cout << "node values are: " << endl;
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);
          //cout << x_val << endl;
          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            //num_forts_variable_in[i]++;
            forbidden_nodes.insert(i);
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */
              if(optimum > 1.0){
        double weight_sum = 0;
        for(int i=0; i<num_nodes; ++i){
          weight_sum += weights[i];
        }
        //cout << "sum of weights is " << weight_sum << endl;
      }




          if(optimum < max_weight && new_fort.size() >= 1){
            //cout << "inserting new fort" << endl;
            //std::pair<std::set<std::set<int> >::iterator, bool> temp;
            forts.insert(new_fort);
            /*if(!temp.second){
              cout << "insertion failed" << endl;
              cout << "new_fort is ";
              for(std::set<int>::iterator it = new_fort.begin(); it != new_fort.end(); ++it){
                cout << *it << ",";
              }
              cout << endl;
            }*/
          }
          else{
            //cout << "did not insert new fort" << endl;
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };


  //Finds violated fort constraints
  int find_minimal_sparse_border_fort_LP(GRBenv* env, std::set<std::set<int> >& forts, std::set<int>& forbidden_nodes, std::vector<double>& weights, double max_weight, std::vector<int>& num_forts_variable_in){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = 2*num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = 0.0;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    for(int i=num_nodes; i<2*num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = 1.0*SPARSIFY_PENALTY;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    //Ensure that forbidden nodes are not chosen
    if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      //cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*constraint for border variables*/
    {
      for(int i=0; i< num_nodes; ++i){
        std::vector<int> cind;
        std::vector<double> cval;
        const int node_degree = nodes[i].get_degree();
        cind.push_back(i);
        cind.push_back(num_nodes + i);
        cval.push_back((double) node_degree);
        cval.push_back(-1.0*node_degree);
        for(int j=0; j<node_degree; ++j){
          cind.push_back(nodes[i].get_adj_list_member(j));
          cval.push_back(-1.0);
        }


      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 0.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }


      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
      }
    }//*/

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);


        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      continue_flag = false;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      //cout << "status is " << status << endl;
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        //cout << "optimum is : " << optimum << endl;
        //cout << "node values are: " << endl;
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);
          //cout << x_val << endl;
          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            //num_forts_variable_in[i]++;
            forbidden_nodes.insert(i);
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */
              if(optimum > 1.0){
        double weight_sum = 0;
        for(int i=0; i<num_nodes; ++i){
          weight_sum += weights[i];
        }
        //cout << "sum of weights is " << weight_sum << endl;
      }




          if(optimum < max_weight && new_fort.size() >= 1){
            //cout << "inserting new fort" << endl;
            //std::pair<std::set<std::set<int> >::iterator, bool> temp;
            forts.insert(new_fort);
            /*if(!temp.second){
              cout << "insertion failed" << endl;
              cout << "new_fort is ";
              for(std::set<int>::iterator it = new_fort.begin(); it != new_fort.end(); ++it){
                cout << *it << ",";
              }
              cout << endl;
            }*/
          }
          else{
            //cout << "did not insert new fort" << endl;
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };





 //Determines whether a fort constraint gives a facet
  int determine_facet_LP(GRBenv* env, const std::set<int>& fort){
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int fort_size = fort.size();
    const int num_variables = (1+fort_size)*num_nodes + fort_size;
    const int y_start = (1+fort_size)*num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));


    for(int i=0; i<num_nodes; ++i){
      if(fort.count(i) != 1){
        types[i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[i] = 0;
        upper_bounds[i] = 1.0;
        lower_bounds[i] = 0.0;
        for(int j=0; j<fort_size; ++j){
          const int index = num_nodes+j*num_nodes+i;
          types[index] = GRB_BINARY;
          objective_coeffs[index] = 0;
          upper_bounds[index] = 1;
          lower_bounds[index] = 0;
        }
      }
      else{
        types[i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[i] = 0;
        upper_bounds[i] = 0.0;
        lower_bounds[i] = 0.0;
        for(int j=0; j<fort_size; ++j){
          const int index = num_nodes+j*num_nodes+i;
          types[index] = GRB_BINARY;
          objective_coeffs[index] = 0;
          upper_bounds[index] = 1;
          lower_bounds[index] = 0;
        }

      }
    }
    for(int i=0; i<fort_size; ++i){
      types[y_start+i] = GRB_BINARY;
      objective_coeffs[y_start+i] = 1;
      upper_bounds[y_start+i] = 1.0;
      lower_bounds[y_start+i] = 0;
    }

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*First constraint requiring some outside vertex to be chosen*/
    {
      std::vector<int> cind;
      std::vector<double> cval;
      for(int i=0; i<num_nodes; ++i){
        cind.push_back(i);
        cval.push_back(1.0);
      }
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    /*Second set of constraints requiring each fort to contain a chosen v*/
    {
      for(int i=0; i<fort_size; ++i){
        std::vector<int> cind;
        std::vector<double> cval;
        for(int j=0; j<num_nodes; ++j){
          if(fort.count(j) !=1){
            cind.push_back((i+1)*num_nodes+j);
            cval.push_back(1.0);
          }
        }
        cind.push_back(y_start+i);
        cval.push_back(-1.0);

        error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 0, NULL);
        if(error){
          cout << "ERROR at: add constraint 2: " << error << endl;
          return 1;
        }
      }
    }

    /*Third set of constraints requiring each fort to chosen v's to be the same v*/
    {
      for(int i=0; i<fort_size; ++i){
        for(int j=0; j<num_nodes; ++j){
          std::vector<int> cind;
          std::vector<double> cval;
          cind.push_back((i+1)*num_nodes+j);
          cval.push_back(1.0);
          cind.push_back(j);
          cval.push_back(-1.0);
          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
        }
      }
    }


    //Fourth set of constraints requiring forts to have empty intersection
    {
      for(int i=0; i<num_nodes; ++i){
        if(fort.count(i) == 1){
          std::vector<int> cind;
          std::vector<double> cval;
          for(int j=0; j<fort_size; ++j){
            cind.push_back(num_nodes + j*num_nodes + i);
            cval.push_back(1.0);
          }

          for(int j=0; j<fort_size; ++j){
            cind.push_back(y_start+j);
            cval.push_back(-1.0);
          }


          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }

        }
      }
    }


    //Fifth set of constraints requiring each fort union v to be a fort
    {
      for(int k=0; k<fort_size; ++k){
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          std::vector<double> cval;
          std::vector<int> cind;
          //first index is v variable for original vertex
          cind.push_back(i);
          cval.push_back(-1.0);

          cind.push_back(num_nodes+k*num_nodes + i);
          cval.push_back(-1.0);
          //neighbor can either be v or in rest of fort
          cind.push_back(neighbor);
          cval.push_back(1.0);
          cind.push_back(num_nodes+k*num_nodes + neighbor);
          cval.push_back(1.0);

          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind.push_back(next_neighbor);
              cval.push_back(1.0);
              cind.push_back(num_nodes+k*num_nodes  + next_neighbor);
              cval.push_back(1.0);
            }
          }

          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
        }
	//cout << endl;
      }
      }
    }

    //Sixth set of constraints requiring positive variables only in chosen forts
    {
      for(int i=0; i<fort_size; ++i){
        for(int j=0; j<num_nodes; ++j){
          std::vector<int> cind;
          std::vector<double> cval;
          cind.push_back((i+1)*num_nodes+j);
          cval.push_back(1.0);
          cind.push_back(y_start+i);
          cval.push_back(-1.0);
          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
        }
      }

    }

     error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
     error = GRBsetintparam(GRBgetenv(model), "Threads", 1);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    GRBoptimize(model);
    if(error){
      cout << "ERROR at: optimize" << endl;
      return 1;
    }
    int status;
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status);
    int v;
    if(status == GRB_OPTIMAL){

      for(int i=0; i<num_nodes; ++i){
        double x_val;
        error = GRBgetdblattrelement(model, "X", i, &x_val);
        if(x_val > 1-FLOAT_TOL){
          v = i;
        }
      }

      /*for(int j=0; j<fort_size; ++j){
        double y_j;
        error = GRBgetdblattrelement(model, "X", y_start+j, &y_j);
        if(y_j > 1-FLOAT_TOL){
          std::set<int> new_set;
          for(int i=0; i<num_nodes; ++i){
            double x_val;
            error = GRBgetdblattrelement(model, "X", i, &x_val);
            if(x_val > 1-FLOAT_TOL){
              new_set.insert(i);
            }
          }
          new_forts.insert(new_set);

        }
      }*/
      error = GRBfreemodel(model);
      return v;
    }
    else if(status == GRB_INFEASIBLE){
      error = GRBfreemodel(model);
      return -1;
    }
    else{
      error = GRBfreemodel(model);
      return -2;
    }



  };

 //Determines whether a fort constraint is likely to be a facet
  int determine_facet_AB_LP(GRBenv* env, const std::set<int>& fort, std::set<int>& A, std::set<int>& B){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int fort_size = fort.size();
    const int num_variables = 3*num_nodes;
    const int A_start = num_nodes;
    const int B_start = A_start + num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      if(fort.count(i) != 1){
        types[i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[i] = SPARSIFY_PENALTY;
        upper_bounds[i] = 1.0;
        lower_bounds[i] = 0.0;
        types[A_start + i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[A_start + i] = SPARSIFY_PENALTY;
        upper_bounds[A_start + i] = 0.0;
        lower_bounds[A_start + i] = 0.0;
        types[B_start + i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[B_start + i] = SPARSIFY_PENALTY;
        upper_bounds[B_start + i] = 0.0;
        lower_bounds[B_start + i] = 0.0;
      }
      else{
        types[i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[i] = SPARSIFY_PENALTY;
        upper_bounds[i] = 0.0;
        lower_bounds[i] = 0.0;
        types[A_start + i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[A_start + i] = SPARSIFY_PENALTY;
        upper_bounds[A_start + i] = 1.0;
        lower_bounds[A_start + i] = 0.0;
        types[B_start + i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[B_start + i] = SPARSIFY_PENALTY;
        upper_bounds[B_start + i] = 1.0;
        lower_bounds[B_start + i] = 0.0;
      }
    }

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*First constraint requiring some outside vertex to be chosen*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring A and B to be disjoint
    {

      for(int i=0; i<num_nodes; ++i){
        if(fort.count(i) == 1){
        std::vector<int> cind;
        std::vector<double> cval;
        cind.push_back(A_start + i);
        cind.push_back(B_start + i);
        cval.push_back(1.0);
        cval.push_back(1.0);

        error = GRBaddconstr(model, 2, &cind[0], &cval[0], GRB_LESS_EQUAL, 1, NULL);
        if(error){
          cout << "ERROR at: add constraint 2: " << error << endl;
          return 1;
        }
        }

      }
    }


    //Third set of constraints requiring B union v to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          std::vector<double> cval;
          std::vector<int> cind;
          //first index is v variable for original vertex
          cind.push_back(i);
          cval.push_back(-1.0);

          cind.push_back(B_start + i);
          cval.push_back(-1.0);
          cind.push_back(neighbor);
          cval.push_back(1.0);
          cind.push_back(B_start+ neighbor);
          cval.push_back(1.0);

          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind.push_back(next_neighbor);
              cval.push_back(1.0);
              cind.push_back(B_start + next_neighbor);
              cval.push_back(1.0);
            }
          }

          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
        }
	//cout << endl;
      }
    }

        //Fourth set of constraints requiring A union v to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          std::vector<double> cval;
          std::vector<int> cind;
          //first index is v variable for original vertex
          cind.push_back(i);
          cval.push_back(-1.0);

          cind.push_back(A_start + i);
          cval.push_back(-1.0);
          cind.push_back(neighbor);
          cval.push_back(1.0);
          cind.push_back(A_start+ neighbor);
          cval.push_back(1.0);

          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind.push_back(next_neighbor);
              cval.push_back(1.0);
              cind.push_back(A_start + next_neighbor);
              cval.push_back(1.0);
            }
          }

          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
        }
	//cout << endl;
      }
    }


    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBsetintparam(GRBgetenv(model), "Threads", 1);

    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    GRBoptimize(model);
    if(error){
      cout << "ERROR at: optimize" << endl;
      return 1;
    }
    int status;
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status);
    int v;
    if(status == GRB_OPTIMAL){

      for(int i=0; i<num_nodes; ++i){
        double x_val;
        error = GRBgetdblattrelement(model, "X", i, &x_val);
        if(x_val > 1-FLOAT_TOL){
          v = i;
        }
        error = GRBgetdblattrelement(model, "X", i + num_nodes, &x_val);
        if(x_val > 1-FLOAT_TOL){
          A.insert(i);
        }
        error = GRBgetdblattrelement(model, "X", i + 2*num_nodes, &x_val);
        if(x_val > 1-FLOAT_TOL){
          B.insert(i);
        }

      }
      error = GRBfreemodel(model);
      return v;
    }
    else if(status == GRB_INFEASIBLE){
      error = GRBfreemodel(model);
      return -1;
    }
    else{
      error = GRBfreemodel(model);
      return -2;
    }


  };



 //Finds violated fort constraints
  int determine_facet_AB_LP(GRBenv* env, const std::set<int>& fort){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int fort_size = fort.size();
    const int num_variables = 3*num_nodes;
    const int A_start = num_nodes;
    const int B_start = A_start + num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      if(fort.count(i) != 1){
        types[i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[i] = SPARSIFY_PENALTY;
        upper_bounds[i] = 1.0;
        lower_bounds[i] = 0.0;
        types[A_start + i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[A_start + i] = SPARSIFY_PENALTY;
        upper_bounds[A_start + i] = 0.0;
        lower_bounds[A_start + i] = 0.0;
        types[B_start + i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[B_start + i] = SPARSIFY_PENALTY;
        upper_bounds[B_start + i] = 0.0;
        lower_bounds[B_start + i] = 0.0;
      }
      else{
        types[i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[i] = SPARSIFY_PENALTY;
        upper_bounds[i] = 0.0;
        lower_bounds[i] = 0.0;
        types[A_start + i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[A_start + i] = SPARSIFY_PENALTY;
        upper_bounds[A_start + i] = 1.0;
        lower_bounds[A_start + i] = 0.0;
        types[B_start + i] = GRB_BINARY;
        //types[i] = GRB_CONTINUOUS;
        objective_coeffs[B_start + i] = SPARSIFY_PENALTY;
        upper_bounds[B_start + i] = 1.0;
        lower_bounds[B_start + i] = 0.0;
      }
    }

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*First constraint requiring some outside vertex to be chosen*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring A and B to be disjoint
    {
      for(int i=0; i<num_nodes; ++i){
        if(fort.count(i) ==1){
        std::vector<int> cind;
        std::vector<double> cval;
        cind.push_back(A_start + i);
        cind.push_back(B_start + i);
        cval.push_back(1.0);
        cval.push_back(1.0);

        error = GRBaddconstr(model, 2, &cind[0], &cval[0], GRB_LESS_EQUAL, 1, NULL);
        if(error){
          cout << "ERROR at: add constraint 2: " << error << endl;
          return 1;
        }

        }
      }
    }


    //Third set of constraints requiring B union v to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          std::vector<double> cval;
          std::vector<int> cind;
          //first index is v variable for original vertex
          cind.push_back(i);
          cval.push_back(-1.0);

          cind.push_back(B_start + i);
          cval.push_back(-1.0);
          cind.push_back(neighbor);
          cval.push_back(1.0);
          cind.push_back(B_start+ neighbor);
          cval.push_back(1.0);

          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind.push_back(next_neighbor);
              cval.push_back(1.0);
              cind.push_back(B_start + next_neighbor);
              cval.push_back(1.0);
            }
          }

          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
        }
	//cout << endl;
      }
    }

        //Fourth set of constraints requiring A union v to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          std::vector<double> cval;
          std::vector<int> cind;
          //first index is v variable for original vertex
          cind.push_back(i);
          cval.push_back(-1.0);

          cind.push_back(A_start + i);
          cval.push_back(-1.0);
          cind.push_back(neighbor);
          cval.push_back(1.0);
          cind.push_back(A_start+ neighbor);
          cval.push_back(1.0);

          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind.push_back(next_neighbor);
              cval.push_back(1.0);
              cind.push_back(A_start + next_neighbor);
              cval.push_back(1.0);
            }
          }

          error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
        }
	//cout << endl;
      }
    }


    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBsetintparam(GRBgetenv(model), "Threads", 1);

    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    GRBoptimize(model);
    if(error){
      cout << "ERROR at: optimize" << endl;
      return 1;
    }
    int status;
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status);

    if(status == GRB_OPTIMAL){

      for(int i=0; i<num_nodes; ++i){
        double x_val;
        error = GRBgetdblattrelement(model, "X", i, &x_val);
        if(x_val > 1-FLOAT_TOL){
          error = GRBfreemodel(model);
          return i;
        }
      }
    }
    else if(status == GRB_INFEASIBLE){
      error = GRBfreemodel(model);
      return -1;
    }
    else{
      error = GRBfreemodel(model);
      return -2;
    }


  };


 //Finds violated fort constraints
  int find_minimal_fort_LP(GRBenv* env, std::set<std::set<int> >& forts, std::set<int>& forbidden_nodes, std::vector<double>& weights, double max_weight){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = 1.0; //weights[i] + SPARSIFY_PENALTY;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    //Ensure that forbidden nodes are not chosen
    if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      //cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);


        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBsetintparam(GRBgetenv(model), "Threads", 1);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      continue_flag = false;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      //cout << "status is " << status << endl;
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        //cout << "optimum is : " << optimum << endl;
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);
          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            forbidden_nodes.insert(i);
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */
              if(optimum > 1.0){
        double weight_sum = 0;
        for(int i=0; i<num_nodes; ++i){
          weight_sum += weights[i];
        }
        //cout << "sum of weights is " << weight_sum << endl;
      }

          if(optimum < max_weight && new_fort.size() >= 1){
            forts.insert(new_fort);
          }
          else{
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };

 //Finds violated fort constraints
  int find_max_intersection_fort_LP(GRBenv* env, std::set<std::set<int> >& forts, std::vector<int>& num_forts_variable_in, const int forts_size, std::vector<double>& weights, double max_weight){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      //cout << ((double)num_forts_variable_in[i])/(100.0*forts_size) << endl;
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = weights[i] + SPARSIFY_PENALTY - ((double)num_forts_variable_in[i])/(100.0*forts_size);
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    //Ensure that forbidden nodes are not chosen
    /*if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }*/

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << error << endl;
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

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);


        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      continue_flag = false;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      //cout << "status is " << status << endl;
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        //cout << "optimum is : " << optimum << endl;
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);
          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            //forbidden_nodes.insert(i);
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */
              if(optimum > 1.0){
        double weight_sum = 0;
        for(int i=0; i<num_nodes; ++i){
          weight_sum += weights[i];
        }
        //cout << "sum of weights is " << weight_sum << endl;
      }

          if(optimum < max_weight && new_fort.size() >= 1){
            forts.insert(new_fort);
          }
          else{
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };


  //Finds violated fort constraints
  int find_maximal_fort_LP(GRBenv* env, std::set<std::set<int> >& forts, std::set<int>& forbidden_nodes, std::vector<double>& weights, double max_weight){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = weights[i] - SPARSIFY_PENALTY;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    //Ensure that forbidden nodes are not chosen
    /*if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }*/

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    /*constraint requiring fort to leave out at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_LESS_EQUAL, (double) num_nodes - 1, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);


        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
    error = GRBupdatemodel(model);
    if(error){
      cout << "ERROR at: third model update" << endl;
      return 1;
    }
    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      continue_flag = false;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      //cout << "status is " << status << endl;
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        double penalty_sum = 0;
        //cout << "optimum is : " << optimum << endl;
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);

          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            forbidden_nodes.insert(i);
            penalty_sum += SPARSIFY_PENALTY;
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */

          if(optimum + penalty_sum < max_weight && new_fort.size() >= 1){
            forts.insert(new_fort);
          }
          else{
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };


  //Finds violated fort constraints
  int find_forts_LP(GRBenv* env, std::set<std::set<int> >& forts, std::set<int>& forbidden_nodes, std::vector<double>& weights, double max_weight){
    //cout << "STARTING TO FIND FORTS" << endl;
    const int num_nodes = nodes.size();
    GRBmodel *model;
    int error;
    const int num_variables = num_nodes;
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<num_nodes; ++i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      objective_coeffs[i] = 1.0; //weights[i] + SPARSIFY_PENALTY;
      upper_bounds[i] = 1.0;
      lower_bounds[i] = 0.0;
    }
    //Ensure that forbidden nodes are not chosen
    if(!forbidden_nodes.empty()){
    for(std::set<int>::iterator it = forbidden_nodes.begin(); it !=forbidden_nodes.end(); ++it){
      cout << "setting " << *it << " to 0 " << endl;
      upper_bounds[*it] = 0.0;
    }
    }

    error = GRBnewmodel(env, &model, "model", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
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

    /*First constraint requiring fort to have at least 1 vertex*/
    {
      double* cval = (double*) calloc(num_nodes, sizeof(double));
      int* cind = (int*) calloc(num_nodes, sizeof(int));
      for(int i=0; i<num_nodes; ++i){
        cind[i] = i;
        cval[i] = 1.0;
      }
      error = GRBaddconstr(model, num_nodes, cind, cval, GRB_GREATER_EQUAL, 1.0, NULL);
      if(error){
	cout << "ERROR at: add constraint 1: " << error << endl;
	return 1;
      }
      free(cind);
      free(cval);

      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: second model update" << endl;
	return 1;
      }
    }//*/

    //Second set of constraints requiring fort to be a fort
    {
      for(int i=0 ; i<num_nodes; ++i){
        const int node_degree = nodes[i].get_degree();
        //cout << "node is " << i << " with neighbors ";
        for(int j=0; j<node_degree; ++j){
          const int neighbor = nodes[i].get_adj_list_member(j);
          //cout << neighbor << " " << endl;
          const int neighbor_degree = nodes[neighbor].get_degree();
          //cout << "neighbor_degree is " << neighbor_degree << endl;
          double* cval = (double*) calloc(neighbor_degree+1, sizeof(double));
          int* cind = (int*) calloc(neighbor_degree+1, sizeof(int));
          //first index is original vertex
          cind[0] = i;
          cval[0] = -1.0;
          //second index is neighbor
          cind[1] = neighbor;
          cval[1] = 1.0;
          int counter = 2;
          //rest of vertices are neighbors neighbors
          for(int k=0; k<neighbor_degree; ++k){
            const int next_neighbor = nodes[neighbor].get_adj_list_member(k);
            if(next_neighbor != i){
              //cout << "second neighbor is " << next_neighbor << endl;
              cind[counter] = next_neighbor;
              cval[counter] = 1.0;
              counter++;
            }
          }

          error = GRBaddconstr(model, neighbor_degree+1, cind, cval, GRB_GREATER_EQUAL, 0, NULL);
          if(error){
            cout << "ERROR at: add constraint 2: " << error << endl;
            return 1;
          }
          free(cind);
          free(cval);

          error = GRBupdatemodel(model);
          if(error){
            cout << "ERROR at: third model update" << endl;
            return 1;
          }
        }
	//cout << endl;
      }
    }
    error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);

    bool continue_flag = true;
    int counter = 0;
    while(continue_flag){
      continue_flag = false;
      GRBoptimize(model);
      if(error){
	cout << "ERROR at: optimize" << endl;
	return 1;
      }
      //counter++;


      //cout << "objective value is " << optimum << endl;
      /*if(optimum > max_weight - FLOAT_TOL || counter >= MAX_FORT_NUMBER){
        continue_flag = false;
      }
      */

      std::set<int> new_fort;
      //std::vector<int> cind;
      //std::vector<double> cval;
      int status;
      error = GRBgetintattr(model, "Status", &status);
      if(status == GRB_OPTIMAL){
        double optimum;
        error = GRBgetdblattr(model, "ObjVal", &optimum);
        for(int i=0; i<num_nodes; ++i){
          double x_val;

          error = GRBgetdblattrelement(model, "X", i, &x_val);
          //optimum -= x_val*((double)num_forts_node_in[i]/1000.0);
          if(x_val >= 1- FLOAT_TOL){
            new_fort.insert(i);
            forbidden_nodes.insert(i);
            //cout << "num_forts_node_in is " << num_forts_node_in[i] << endl;
            //cout << "optimum is " << optimum << endl;
            //cind.push_back(i);
            //cval.push_back(1.0);
          }
        }
        error = GRBsetdblattrelement(model, "UB", *new_fort.begin(), 0.0);
        /*if(new_fort.size() > MAX_FORT_SIZE){
          continue_flag = false;
        }
        cind.push_back(*(new_fort.begin()));
        cval.push_back(1.0);
        */
          if(optimum < max_weight && new_fort.size() >= 1){
            forts.insert(new_fort);
          }
          else{
            continue_flag = false;
          }
        }
        else{
          continue_flag = false;
        }
      }
      //cout << "constrained node is: " << cind[0] << endl;
      /*eliminate that fort from consideration
      error = GRBaddconstr(model, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 0, NULL);
      if(error){
	cout << "ERROR at: add fort elimination constraint: " << error << endl;
	return 1;
      }
      error = GRBupdatemodel(model);
      if(error){
	cout << "ERROR at: model update" << endl;
	return 1;
      }
      //}*/
    error = GRBfreemodel(model);
    return 0;
  };


  //subgraph_cut is the incoming cut
  void get_subgraph(const std::set<int>& nodes_in_subgraph, forcing::graph& subgraph, edge_list& subgraph_cut){

    std::unordered_map<int, int> node_map;
    const int num_nodes_subgraph = nodes_in_subgraph.size();
    const int num_nodes_supergraph = nodes.size();
    const int num_edges_supergraph = edges.get_size();
    int counter = 0;
    const std::set<int>::iterator nodes_in_subgraph_end = nodes_in_subgraph.end();
    for(std::set<int>::iterator it = nodes_in_subgraph.begin(); it != nodes_in_subgraph_end; ++it){
      node current_node;
      const int current_node_former_degree = nodes[*it].get_degree();
      for(int i=0; i<current_node_former_degree; ++i){
        const int neighbor = nodes[*it].get_adj_list_member(i);
        if(nodes_in_subgraph.count(neighbor) > 0){
          current_node.add_to_adj_list(neighbor);
          subgraph.edges.add_edge(*it, neighbor);
        }
        else{
          subgraph_cut.add_edge(neighbor, *it);
        }
      }
      subgraph.nodes.push_back(current_node);
      node_map.emplace(*it, counter);
      counter++;
    }

    //Subgraph is created, now we adapt the nodes and edges to the new node numbers in the subgraph.
    /*cout << "subgraph nodes are: " << endl;
    for(int i=0; i<num_nodes_subgraph; ++i){
      const int node_degree = subgraph.nodes[i].get_degree();
      cout << "Node " << i << " : ";
      for(int j=0; j<node_degree; ++j){
	cout << subgraph.nodes[i].get_adj_list_member(j) << ", ";
      }
      cout << endl;
      }*/

    for(int i=0; i<num_nodes_subgraph; ++i){
      std::vector<int> new_adj_list;
      const int node_degree = subgraph.nodes[i].get_degree();
      for(int j=0; j<node_degree; ++j){
        new_adj_list.push_back(node_map[subgraph.nodes[i].get_adj_list_member(j)]);
      }
      subgraph.nodes[i].swap_adj_lists(new_adj_list);
    }

    /*cout << "Adjusted subgraph nodes are: " << endl;
    for(int i=0; i<num_nodes_subgraph; ++i){
      const int node_degree = subgraph.nodes[i].get_degree();
      cout << "Node " << i << " : ";
      for(int j=0; j<node_degree; ++j){
	cout << subgraph.nodes[i].get_adj_list_member(j) << ", ";
      }
      cout << endl;
      }*/

    const int num_edges_subgraph = subgraph.edges.get_size();
    std::vector<std::pair<int, int> > new_edges;
    for(int i =0; i<num_edges_subgraph; ++i){
      new_edges.push_back(std::make_pair(node_map[subgraph.edges.get_end1(i)],node_map[subgraph.edges.get_end2(i)]));
    }
    subgraph.edges.swap_edge_list(new_edges);
  };

  //subgraph_cut is the incoming cut
  void get_subgraph(const std::set<int>& nodes_in_subgraph, forcing::graph& subgraph, edge_list& subgraph_cut, std::unordered_map<int, int>& reverse_node_map){

    std::unordered_map<int, int> node_map;
    const int num_nodes_subgraph = nodes_in_subgraph.size();
    const int num_nodes_supergraph = nodes.size();
    const int num_edges_supergraph = edges.get_size();
    int counter = 0;
    const std::set<int>::iterator nodes_in_subgraph_end = nodes_in_subgraph.end();
    for(std::set<int>::iterator it = nodes_in_subgraph.begin(); it != nodes_in_subgraph_end; ++it){
      node current_node;
      const int current_node_former_degree = nodes[*it].get_degree();
      for(int i=0; i<current_node_former_degree; ++i){
        const int neighbor = nodes[*it].get_adj_list_member(i);
        if(nodes_in_subgraph.count(neighbor) > 0){
          current_node.add_to_adj_list(neighbor);
          subgraph.edges.add_edge(*it, neighbor);
        }
        else{
          subgraph_cut.add_edge(neighbor, *it);
        }
      }
      subgraph.nodes.push_back(current_node);
      node_map.emplace(*it, counter);
      reverse_node_map.emplace(counter, *it);
      counter++;
    }

    //Subgraph is created, now we adapt the nodes and edges to the new node numbers in the subgraph.
    /*cout << "subgraph nodes are: " << endl;
    for(int i=0; i<num_nodes_subgraph; ++i){
      const int node_degree = subgraph.nodes[i].get_degree();
      cout << "Node " << i << " : ";
      for(int j=0; j<node_degree; ++j){
	cout << subgraph.nodes[i].get_adj_list_member(j) << ", ";
      }
      cout << endl;
      }*/

    for(int i=0; i<num_nodes_subgraph; ++i){
      std::vector<int> new_adj_list;
      const int node_degree = subgraph.nodes[i].get_degree();
      for(int j=0; j<node_degree; ++j){
        new_adj_list.push_back(node_map[subgraph.nodes[i].get_adj_list_member(j)]);
      }
      subgraph.nodes[i].swap_adj_lists(new_adj_list);
    }

    /*cout << "Adjusted subgraph nodes are: " << endl;
    for(int i=0; i<num_nodes_subgraph; ++i){
      const int node_degree = subgraph.nodes[i].get_degree();
      cout << "Node " << i << " : ";
      for(int j=0; j<node_degree; ++j){
	cout << subgraph.nodes[i].get_adj_list_member(j) << ", ";
      }
      cout << endl;
      }*/

    const int num_edges_subgraph = subgraph.edges.get_size();
    std::vector<std::pair<int, int> > new_edges;
    for(int i =0; i<num_edges_subgraph; ++i){
      new_edges.push_back(std::make_pair(node_map[subgraph.edges.get_end1(i)],node_map[subgraph.edges.get_end2(i)]));
    }
    subgraph.edges.swap_edge_list(new_edges);
  };



};


};


struct callback_data{
  bool not_added;
  int max_timesteps;
  int edge_variable_start;
  int vertex_time_start;
  double time_in_connectivity;
  double time_in_LP;
  double previous_time;
  int min_degree;
  int a_start;
  int ae_start;
  int y_start;
  int l_vertex_start;
  int end_variable_start;
  int num_variables;
  int num_forts;
  clock_t start_time;
  std::chrono::time_point<std::chrono::high_resolution_clock> chrono_start_time;
  int num_nodes;
  int num_edges;
  int edge_forward_start;
  int edge_backward_start;
  int edge_indicator_start;
  int counter;
  int previous_size;
  std::vector<std::vector<int> > forts_node_in;
  std::vector<std::set<int> > forts;
  std::vector<std::set<int> > fort_borders;
  forcing::edge_list edges;
  std::vector<int> orbits;
  std::set<int> nodes_in_cover;
  std::vector<forcing::node> nodes;
  std::set<int> forbidden_nodes;
  std::vector<int> forbidden_fort_indices;
  std::vector<int> num_forts_variable_in;
  std::vector<forcing::node_subset> node_subsets;
  forcing::graph our_graph;
  GRBenv* env;
};


#endif
