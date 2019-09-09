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
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_fem.h>
#include <getfem/getfem_assembling.h>
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_superlu.h>
#include <getfem/getfem_interpolated_fem.h>
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_iter_solvers.h>
#include <node.hpp>
#include <mesh1d_prova_new_order.hpp>
#include <gmm/gmm_dense_qr.h>
using namespace femg;
using namespace BGLgeom;

int main( int argc, char* argv[] ) {
  //Check if number of args is correct
  if (argc != 2)
    std::cerr << "Error: invalid input arguments." << std::endl;

  //Reading file .pts
  std::string out_filename = argv[1];
  std::ifstream ist{out_filename};

  if (!ist)
  	std::cerr << "Something wrong in accessing input file." << std::endl;
  getfem::mesh mh1d;
  std::vector<getfem::node> BCList;
  std::vector<size_t> Nn;
  std::string s = "";
  std::cout << "[GetFEM++] building mesh..." << std::endl;
  getfem::import_pts_file<std::vector<size_t>>(ist, mh1d, BCList, Nn, s);
  bgeot::size_type n_points = mh1d.nb_points();

  //Test for mesh exporting
  //std::string out_prova = "prova.txt";
  //std::cout << "Writing on file " << out_prova << std::endl;
  //mh1d.write_to_file(out_prova);

  //Definition of FEM Method
  getfem::pfem pf = getfem::fem_descriptor("FEM_PK(1,1)");
  getfem::mesh_fem mf(mh1d);
  mf.set_finite_element(pf);

  //getfem::mesh_fem mf_data(mh1d);
  //mf_data.set_finite_element(pf);

  //Definition of Integration Method
  getfem::pintegration_method ppi = getfem::int_method_descriptor("IM_GAUSS1D(7)");
  getfem::mesh_im mim(mh1d);
  mim.set_integration_method(ppi);

  //Computation of mass matrix
  std::cout << "[GetFEM++] computing Mass Matrix..." << std::endl;
  // ACHTUNG! stiamo passando matrice da leggere all'assembly meglio passare
  // wsvector e poi usare copy in rsvector
  gmm::row_matrix<gmm::rsvector<bgeot::scalar_type>> MM(n_points,n_points);
  getfem::asm_mass_matrix(MM, mim, mf);
  gmm::csc_matrix<double> M2(n_points,n_points);
  gmm::clean(MM, 1E-12);
  gmm::copy(MM, M2);
  gmm::MatrixMarket_save("mass_matrix.txt", M2);

  //Computation of stiffness matrix
  gmm::row_matrix<gmm::rsvector<bgeot::scalar_type>> SM(n_points,n_points);

  std::cout << "[GetFEM++] computing Stiffness Matrix..." << std::endl;
  std::vector<double> A(n_points, 1);
  getfem::asm_stiffness_matrix_for_laplacian(SM, mim, mf, mf, A);
  gmm::csc_matrix<double> M3(n_points,n_points);
  gmm::clean(SM, 1E-12);
  gmm::copy(SM, M3);
  gmm::MatrixMarket_save("stiff_matrix.txt", M3);

  std::cout << "[GetFEM++] printing BC..." << std::endl;
  for (int i=0; i<BCList.size(); ++i){
     std::cout << BCList[i] << std::endl;
  }

// PROVA SEGA
// gmm::col_matrix< gmm::wsvector<double> > M_prova(5, 20);
// M_prova(3, 4) = 5.0;
// M_prova(3, 2) = 2.2;
// M_prova(2, 2) = 2.2;
// M_prova(4, 3) = 2.2;
// M_prova(4, 4) = 2.2;
// std::cout << gmm::sub_matrix(M_prova, gmm::sub_interval(2, 3), gmm::sub_interval(2, 3))
//           << std::endl;
// gmm::col_matrix< gmm::wsvector<double> > M_sub(3, 3);
// gmm::copy(gmm::sub_matrix(M_prova, gmm::sub_interval(2, 3), gmm::sub_interval(2, 3)), M_sub);
//
// std::cout << M_prova <<std::endl;
// std::cout << M_sub <<std::endl;


// Calcolo gli autovalori del problema del paper in maniera burina,
// ovvero inverto brutalmente la matrice del RHS e la moltiplico a sinistra in entrambi
// i membri poi calcolo gli autovalori della nuova matrice

// PROBLEMA: è ovvio che ora è un sistema 21 per 21 quindi trovo 21 autovalori cosa devo fare per
// ricondurmi ai cinque che trova corrispondenti ai vertici reali?
//Cosa hanno fatto loro?
std::vector<double> qr_eigenvalues(n_points, 1);

gmm::dense_matrix<double> Inverse_Mass(n_points, n_points);
gmm::copy(MM,Inverse_Mass);

std::cout << gmm::condition_number(Inverse_Mass)<<std::endl;
gmm::lu_inverse(Inverse_Mass);
gmm::dense_matrix<double> HM(n_points, n_points);
gmm::add(SM,MM,HM);
gmm::dense_matrix<double> eig_M(n_points, n_points);
gmm::mult(Inverse_Mass,HM,eig_M);
double tol = 1E-16;
std::cout << "ci siamo?..." << std::endl;

gmm::implicit_qr_algorithm(eig_M, qr_eigenvalues, tol);
//gmm::extract_eig(eig_M, qr_eigenvalues,tol);
std::sort(qr_eigenvalues.begin(),qr_eigenvalues.end());

std::cout << "[GetFEM++] printing the eigenvalues of the problem in the paper calculated with QR method..." << std::endl;

for(int i=0; i<boost::num_vertices(G); ++i){
   std::cout << qr_eigenvalues[i] << std::endl;
}

//Qua ho voluto calcolare gli autovalori della matrice di Stifness per controllarli
// su matlab e vengono identici (non sono ordinati e non li ho guardati proprio
// uno per uno, ma sono davvero tutti uguali fino a buone cifre dopo la virgola)

std::vector<double> qr_lpl_eig(n_points,1);
gmm::implicit_qr_algorithm(SM, qr_lpl_eig,tol);

std::sort(qr_lpl_eig.begin(),qr_lpl_eig.end());

std::cout << "[GetFEM++] printing the eigenvalues of the stifness matrix (Laplacian) calculated with QR method..." << std::endl;

for(int i=0; i<boost::num_vertices(G); ++i){
   std::cout << qr_lpl_eig[i] << std::endl;
}

return 0;
};
