#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstddef>


#include <boost/graph/adjacency_list.hpp>
#include <writer_pts.hpp>
#include <reader_femg.hpp>
#include <graph_builder.hpp>
#include <graph_access.hpp>

using namespace femg;
using namespace BGLgeom;

int main( int argc, char* argv[] ) {
  //Check if number of args is correct
  if (argc != 4)
    std::cerr << "Error: invalid input arguments." << std::endl;
  //Discretization parameter N
  const int N = std::stoi(argv[1]);
  const unsigned int dim = 3;
  const unsigned int num_BC = 1;
  //Type-aliases
  using Graph = boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, Vertex_base_property<dim,num_BC>, Edge_base_property<linear_geometry<dim>, dim>>;
  using Graph_reader = reader_femg<dim, num_BC,Vertex_base_property<dim,num_BC>, Edge_base_property<linear_geometry<dim>, dim>, no_topological_data>;
  //Graph and reader declarations
  Graph G;
  Graph_reader reader(argv[2]);
  const bool check = true; //check for node uniqueness
  unsigned int count = 0; //edge index
  //Reading file .txt
  std::cout << "Acquiring graph information..." << std::endl;
  while (!reader.is_eof()) {
    Vertex_base_property<dim,num_BC> source;
    Vertex_base_property<dim,num_BC> target;
    //Edge_base_property<linear_geometry<3>, 3> edge;
    Edge_desc<Graph> e;
    Vertex_desc<Graph> s,d;
    reader.get_data();
    source = reader.get_source_data(); //acquiring data...
    target = reader.get_target_data();
    //edge = reader.get_edge_data();
    s = new_vertex(source, G, check);
    d = new_vertex(target, G, check); //modifying the graph
    e = new_linear_edge(s, d, G);
    G[e].mesh.uniform_mesh(N, G[e].geometry); //meshing the edge
    G[e].index = count; //edge index
    count++;
  }
  //Exporting file .pts
  std::cout << "Exporting file .pts..." << std::endl;
  std::string out_filename = argv[3];
  writer_pts<Graph, dim, num_BC> W(out_filename);
  W.export_pts(G);

  std::ifstream ist{out_filename};

return 0;
};
