
/*********************************************************************************
This program uses a brute force method to solve the connected zero-forcing problem

Author: Caleb Fast & Boris Brimkov

Compilation: g++ -std=c++11 -o brute_force_connected_forcing brute_force_connected_forcing.cpp
Usage: ./brute_force_connected_forcing <path to .edg file>

**********************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>
#include <climits>
#include <cfloat>
#include <math.h>
#include <time.h>
#include <fstream>
#include<chrono>

#define STOPPING_TIME 7200
#define FLOAT_TOL 0.0001
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



}


bool check_connectivity(std::set<int>& chosen_nodes, forcing::graph& our_graph){

    forcing::graph subgraph;
    forcing::edge_list subgraph_cut;
    our_graph.get_subgraph(chosen_nodes, subgraph, subgraph_cut);

    const int num_nodes = subgraph.nodes.size();
    const int num_edges = subgraph.edges.get_size();
    std::vector<double> adj_matrix(num_nodes*num_nodes, 0);

    for(int i=0; i<num_edges; ++i){
      const int end1 = subgraph.edges.get_end1(i);
      const int end2 = subgraph.edges.get_end2(i);
      adj_matrix[end1*num_nodes+end2] = 1;
    }

    std::set<int> min_cut;
    std::set< std::set<int> > cuts_under_one;
    double min_cut_value;
    forcing::node_identification_min_cut(num_nodes,adj_matrix, min_cut, min_cut_value, cuts_under_one);

    if(min_cut_value < 1 - FLOAT_TOL){
      return false;
    }
    else{
      return true;
    }

};


struct neighbor
{
  int neighbor_number;
  struct neighbor* link;
};

void add(int x, struct neighbor** vertex )
{
  struct neighbor* p= (struct neighbor*) malloc(sizeof(struct neighbor));
  (*p).neighbor_number=x;
  (*p).link=*vertex;
  *vertex=p;
}


int jthNeighbor(struct neighbor* vertex, int j)
{
  int count=0;
  struct neighbor* current=vertex;
  while(count<j)
  {
    count=count+1;
    current=(*current).link;
  }
  return (*current).neighbor_number;
}

int count(struct neighbor* vertex)
{
  int count=0;
  struct neighbor* current=vertex;
  while(current != NULL)
  {
    count=count+1;
    current=(*current).link;
  }
  return count;
}


void print(struct neighbor* vertex)
{
  struct neighbor* current=vertex;
  printf("\nThe vertex has neighbors: ");
  while(current != NULL)
  {
    printf("%d ", (*current).neighbor_number);
    current=(*current).link;
  }
  printf("\n\n");
}

int main(int argc, char * argv[]){

  if(argc!=2) printf("Error: input file");
   string basename = argv[1];
	basename = basename.substr(basename.find_last_of('/')+1);
  ofstream fout;
  FILE *edgelist_file = fopen(argv[1], "r");
  int n,M;
  fscanf(edgelist_file, "%i %i",&n, &M);

  int i;
  int num_nodes = n;
  forcing::node blank_node;
  std::vector<forcing::node> nodes(num_nodes, blank_node);

  struct neighbor ** v= (struct neighbor**) malloc(n*sizeof(struct neighbor *));

  for (i=0;i<n;i++){
    v[i]=NULL;
  }
  forcing::graph our_graph;
  int end1,end2;
  for (i=0;i<M;i++){
    fscanf(edgelist_file,"%i %i",&end1, &end2);
    add(end1, &v[end2]);
    add(end2, &v[end1]);
    nodes[end1].add_to_adj_list(end2);
    nodes[end2].add_to_adj_list(end1);
    our_graph.edges.add_edge(end1,end2);
    our_graph.edges.add_edge(end2,end1);
  }
  our_graph.nodes = nodes;
  int* inS= (int*) malloc(n*sizeof(int));
/////////////////////////////////////////////////////////////////////////////////////////////////////


  auto chrono_time_start = std::chrono::high_resolution_clock::now();
  clock_t t;
  t = clock();
  clock_t start_time = clock();


  int k;

  for(k=1;k<=n;k++){

    int* x= (int*) malloc(k*sizeof(int));
    int i,j,l;

    //initialize first combination as 1,...,k
    std::set<int> chosen_nodes;
    for(i=0;i<k;i++){
      x[i]=i+1;
     chosen_nodes.insert(i);
    }

    bool connected = check_connectivity(chosen_nodes, our_graph);
    int u,u_remember;
    if(connected){
    //first combination

    int numNofVnotinS;

    int something_got_forced=1;

    for (i=0;i<n;i++)
      inS[i]=0;

    for (i=0;i<k;i++)
      inS[x[i]-1]=1;

    while(something_got_forced==1){
    something_got_forced=0;
    for (i=0;i<n;i++){
        if (inS[i]==1){
            numNofVnotinS=0;
            for (j=0;j<count(v[i]);j++){
                u=jthNeighbor(v[i],j);
                if (inS[u]==0){
                    numNofVnotinS=numNofVnotinS+1;
                    u_remember=u;
                }
                if (numNofVnotinS==2)
                    break;
            }
            if (numNofVnotinS==1){
                 // printf("%d got forced\n",u_remember);
                inS[u_remember]=1;
                something_got_forced=1;
            }
	    }
    }
}

    int sum=0;
    for(i=0;i<n;i++)
      sum=sum+inS[i];
    if(sum==n){
  for(i=0;i<k;i++)
    printf("%d ",x[i]-1);
  printf("\n");

    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
        auto chrono_time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> chrono_time = chrono_time_end - chrono_time_start;
    printf("this took %f seconds to execute \n", time_taken);
    free(x);

    fout.open("brute_connected_results.txt", ios::app);
	fout.precision(10);
	fout << endl << endl;
	fout << "Results for: " << basename << endl;
	fout << "ZFS size was: " << k << endl;
  fout << "Self-calculated time was: " << time_taken << endl;
  fout << "Chrono time was: " << chrono_time.count() << endl;
	fout << endl;
	fout.close();

  return 0;
}
}

    while(1){
      if((clock() - start_time)/static_cast<double>(CLOCKS_PER_SEC) > STOPPING_TIME){
        t = clock() - t;
        double time_taken = ((double)t)/CLOCKS_PER_SEC;
            auto chrono_time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> chrono_time = chrono_time_end - chrono_time_start;
        printf("this took %f seconds to execute \n", time_taken);
        free(x);
        fout.open("brute_connected_results.txt", ios::app);
	fout.precision(10);
	fout << endl << endl;
	fout << "Results for: " << basename << endl;
	fout << "ZFS size was: NOT SOLVED > " << k-1 << endl;
  fout << "Self-calculated time was: " << time_taken << endl;
  fout << "Chrono time was: " << chrono_time.count() << endl;
	fout << endl;
	fout.close();
        return 0;
      }
    //x[k-1-p] ... x[k-1] are maxed out; find left-most entry p that's maxed; if p=k, stop
      int p=0;
      int flag=0;
      while(x[k-1-p]==n-p){
        p++;
        if(p>=k){
          flag=1;
          break;
        }
      }
      if(flag==1)
        break;

    //increment one left of left-most unmaxed value and set all subsequent ones in increasing order
      int c,temp;
      temp=x[k-p-1];
      for (c=0;c<=p;c++){
        x[k-p-1+c]=temp+1+c;
      }

      chosen_nodes.clear();
      for(int i=0; i<k; ++i){
        chosen_nodes.insert(x[i] - 1);
      }
      connected = check_connectivity(chosen_nodes, our_graph);
      if(connected){
    //next combinations
      int something_got_forced=1;

      for (i=0;i<n;i++)
        inS[i]=0;

      for (i=0;i<k;i++)
        inS[x[i]-1]=1;

      while(something_got_forced==1){
    something_got_forced=0;
    for (i=0;i<n;i++){
        if (inS[i]==1){
            int numNofVnotinS=0;
            for (j=0;j<count(v[i]);j++){
                u=jthNeighbor(v[i],j);
                if (inS[u]==0){
                    numNofVnotinS=numNofVnotinS+1;
                    u_remember=u;
                }
                if (numNofVnotinS==2)
                    break;
            }
            if (numNofVnotinS==1){
                 // printf("%d got forced\n",u_remember);
                inS[u_remember]=1;
                something_got_forced=1;
            }
	    }
    }
}

      int sum=0;
      for(i=0;i<n;i++)
        sum=sum+inS[i];
      if(sum==n){
  for(i=0;i<k;i++)
    printf("%d ",x[i]-1);
  printf("\n");
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("this took %f seconds to execute \n", time_taken);
     auto chrono_time_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> chrono_time = chrono_time_end - chrono_time_start;
    free(x);

//      fout.open("brute_connected_results.txt", ios::app);
//	fout.precision(10);
//	fout << endl << endl;
//	fout << "Results for: " << basename << endl;
//	fout << "ZFS size was: " << k << endl;
//  fout << "Self-calculated time was: " << time_taken << endl;
//  fout << "Chrono time was: " << chrono_time.count() << endl;
//	fout << endl;
//	fout.close();
  return 0;
}
}
    }

    free(x);
  }
}
