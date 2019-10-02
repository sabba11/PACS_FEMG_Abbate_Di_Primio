#ifndef HH_READER_FEMG_HH
#define HH_READER_FEMG_HH

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <cctype>
#include <vector>

#include <reader_ASCII.hpp>
#include <boundary_conditions.hpp>
#include <point.hpp>
#include <linear_geometry.hpp>
#include <base_properties.hpp>

using namespace BGLgeom;
using BGLgeom::operator>>;

namespace femg {
  template <unsigned int dim, unsigned int num_BC, typename Vertex_prop, typename Edge_prop, typename Topological_data = no_topological_data>
  class reader_femg : public reader_ASCII<Vertex_prop, Edge_prop, Topological_data> {
  public:
    //Constructors
    reader_femg() : reader_ASCII<Vertex_prop, Edge_prop, Topological_data> () {}
    reader_femg(const std::string & _filename)  : reader_ASCII<Vertex_prop, Edge_prop, Topological_data> (_filename) {}
    //Destructor
    virtual ~reader_femg() {};
    //Methods
    virtual void get_data() override;
    virtual Edge_prop get_edge_data() override;
    virtual Vertex_prop get_source_data() override;
    virtual Vertex_prop get_target_data() override;
    virtual Topological_data get_topological_data() override {
      Topological_data t;
      return t;
    }
  private:
    std::vector<boundary_condition> BCs = std::vector<boundary_condition>(2*num_BC);
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
