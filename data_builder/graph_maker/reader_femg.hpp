/*!
	\file reader_femg.hpp
	\author Stefano Abbate
	\author Andrea Di Primio
	\brief File to read*.txt with edges informations into a Graph object.
*/

#ifndef HH_READER_FEMG_HH
#define HH_READER_FEMG_HH

//Standard libraries
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <cctype>
#include <vector>

//BGLgeom libraries
#include <reader_ASCII.hpp>
#include <boundary_conditions.hpp>
#include <point.hpp>
#include <linear_geometry.hpp>
#include <base_properties.hpp>

using namespace BGLgeom;
using BGLgeom::operator>>;

namespace getfem {
  //! Reader FEMG.
  /*!
    This template class is the children of a BGLgeom class: reader_ASCII.
    It is the specialization of the code that reads from a *.txt files the
    informations of the edges written as source and target coordinates and
    then their boundary conditions.
  */

  template <unsigned int dim, unsigned int num_BC, typename Vertex_prop, typename Edge_prop, typename Topological_data = no_topological_data>
  class reader_femg : public reader_ASCII<Vertex_prop, Edge_prop, Topological_data> {
  public:
    //! Default Constructor.
    reader_femg() : reader_ASCII<Vertex_prop, Edge_prop, Topological_data> () {}
    //! Constructor. Takes as an input the name of the file as a string.
    reader_femg(const std::string & _filename)  : reader_ASCII<Vertex_prop, Edge_prop, Topological_data> (_filename) {}
    //!Destructor.
    virtual ~reader_femg() {};
    /*
		+------------------------------------------------------+
		| 1. Methods						              	   |
		+------------------------------------------------------+
		*/
    //! Function to return all the informations of the edge.
    virtual void get_data() override;
    //! Function to return the informations regarding the geometry of the edge.
    virtual Edge_prop get_edge_data() override;
    //! Function to return the informations regarding the source: BC and coordinates.
    virtual Vertex_prop get_source_data() override;
    //! Function to return the informations regarding the target: BC and coordinates.
    virtual Vertex_prop get_target_data() override;
    //! Dummy function to return the informations about the topology. Useless in this reader.
    virtual Topological_data get_topological_data() override {
      Topological_data t;
      return t;
    }

  private:
    /*
    +------------------------------------------------------+
    | 2. Stored informations						              	   |
    +------------------------------------------------------+
    */
    //! Vector containg boundary conditions informations.
    std::vector<boundary_condition> BCs = std::vector<boundary_condition>(2*num_BC);
    //! Vector containing the edge i.e. the source and the target.
    std::vector<point<dim>> extrema = std::vector<point<dim>>(2);
  };

  template <unsigned int dim, unsigned int num_BC, typename Vertex_prop, typename Edge_prop, typename Topological_data>
  void reader_femg<dim, num_BC, Vertex_prop, Edge_prop, Topological_data>::get_data() {
    this->in_file >> extrema[0] >> extrema[1];
    for (unsigned j = 0; j < 2*num_BC; j++) {
      this->in_file >> BCs[j];
    } //reads from txt
    return;
  };

  //Gets edge data. This method is never called.
  template <unsigned int dim, unsigned int num_BC, typename Vertex_prop, typename Edge_prop, typename Topological_data>
  Edge_prop reader_femg<dim, num_BC, Vertex_prop, Edge_prop, Topological_data>::get_edge_data() {
    BGLgeom::linear_geometry<dim> lin_geom(extrema[0], extrema[1]);
    Edge_prop E(lin_geom);
    return E;
  };

  //Gets source data.
  template <unsigned int dim, unsigned int num_BC, typename Vertex_prop, typename Edge_prop, typename Topological_data>
  Vertex_prop reader_femg<dim, num_BC, Vertex_prop, Edge_prop, Topological_data>::get_source_data() {
    std::array<BGLgeom::boundary_condition, num_BC> BC;
    for (unsigned j = 0; j < num_BC; j++) {
      BC[j] = BCs[j];
    } //initialize bc
    Vertex_prop src(extrema[0], BC); //creates vertex
    return src;
  };

  //Gets target data (same code as source)
  template <unsigned int dim, unsigned int num_BC, typename Vertex_prop, typename Edge_prop, typename Topological_data>
  Vertex_prop reader_femg<dim, num_BC, Vertex_prop, Edge_prop, Topological_data>::get_target_data() {
    std::array<BGLgeom::boundary_condition, num_BC> BC;
    for (unsigned j = 0; j < num_BC; j++) {
      BC[j] = BCs[num_BC + j];
    } //initialize bc
    Vertex_prop tgt(extrema[1], BC);
    return tgt;
  };
}

#endif
