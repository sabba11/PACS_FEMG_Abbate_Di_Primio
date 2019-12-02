/*!
	\file ptsbuilder_h.cpp
	\author Stefano Abbate
	\author Andrea Di Primio
	\brief File used to build the meshed graph *.pts with a step of length h.
*/

//Standard libraries
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstddef>

//BGLgeom libraries
#include <boost/graph/adjacency_list.hpp>
#include <writer_pts.hpp>
#include <reader_femg.hpp>
#include <graph_builder.hpp>
#include <graph_access.hpp>

using namespace getfem;
using namespace BGLgeom;

//! Main function.
/*!
  This function takes three input parameters in the command line:
  - The length of the step h in which the user wants to divide the edges.
    We will divide the length L in N step, where N equals:
    \f[
       N = \lfloor \frac{L}{h} \rfloor
    \f]
  - The *.txt containing a matrix with the edges informations.
    Respectively starting, ending points coordinates and then their BCs.
  - The address and name at which the user wants to save the file.
    Use *.pts format.
*/

int main( int argc, char* argv[] ) {
  //Check if number of args is correct
  if (argc != 4)
    std::cerr << "Error: invalid input arguments." << std::endl;
  //Discretization parameter h
  const double h = std::stod(argv[1]);
  const unsigned int num_BC = 1;
  //Checking .txt point dimension
  std::fstream ost_check{argv[2]};
  double coord_check;
  unsigned counter = 0;
  while (ost_check>>coord_check)
    counter++;
  // 2 dim
  if (counter == 4){
    const unsigned int dim = 2;
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
      Edge_desc<Graph> e;
      Vertex_desc<Graph> s,d;
      reader.get_data();
      source = reader.get_source_data(); //acquiring data...
      target = reader.get_target_data();
      s = new_vertex(source, G, check);
      d = new_vertex(target, G, check); //modifying the graph
      e = new_linear_edge(s, d, G);
      // Calculation for correct mesh step for each edge
      double l = G[e].geometry.length();
      if (l<h)
         std::cerr << "Error: mesh step longer than the edge." << std::endl;
      int N = l/h;
      if ((l-N*h)>0.5)
         N++;
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
  //3 dim
  }else if (counter == 6){
    const unsigned int dim3 = 3;
    //Type-aliases
    using Graph3 = boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, Vertex_base_property<dim3,num_BC>, Edge_base_property<linear_geometry<dim3>, dim3>>;
    using Graph_reader3 = reader_femg<dim3, num_BC,Vertex_base_property<dim3,num_BC>, Edge_base_property<linear_geometry<dim3>, dim3>, no_topological_data>;
    //Graph and reader declarations
    Graph3 G;
    Graph_reader3 reader(argv[2]);
    const bool check = true; //check for node uniqueness
    unsigned int count = 0; //edge index
    //Reading file .txt
    std::cout << "Acquiring graph information..." << std::endl;
    while (!reader.is_eof()) {
      Vertex_base_property<dim3,num_BC> source;
      Vertex_base_property<dim3,num_BC> target;
      Edge_desc<Graph3> e;
      Vertex_desc<Graph3> s,d;
      reader.get_data();
      source = reader.get_source_data(); //acquiring data...
      target = reader.get_target_data();
      s = new_vertex(source, G, check);
      d = new_vertex(target, G, check); //modifying the graph
      e = new_linear_edge(s, d, G);
      // Calculation for correct mesh step for each edge
      double l = G[e].geometry.length();
      if (l<h)
         std::cerr << "Error: mesh step longer than the edge." << std::endl;
      int N = l/h;
      if ((l-N*h)>0.5)
         N++;
      G[e].mesh.uniform_mesh(N, G[e].geometry); //meshing the edge
      G[e].index = count; //edge index
      count++;
    }
    //Exporting file .pts
    std::cout << "Exporting file .pts..." << std::endl;
    std::string out_filename = argv[3];
    writer_pts<Graph3, dim3, num_BC> W(out_filename);
    W.export_pts(G);
    std::ifstream ist{out_filename};
  }else
    std::cerr << "Error: invalid input arguments." << std::endl;
return 0;
};
