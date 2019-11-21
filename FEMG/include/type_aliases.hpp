// Header with custom definitions. To be completed.

#ifndef FEMG_TYPE_ALIASES_HPP
#define FEMG_TYPE_ALIASES_HPP

#include <complex>
#include <functional>
#include <gmm/gmm.h>
#include <getfem/getfem_import.h>
#include <getfem/getfem_assembling.h>

namespace getfem {

	//bgeot namespace
	using bgeot::scalar_type;
	using bgeot::size_type;

	typedef std::complex<scalar_type> complex_scalar_type;
	typedef std::vector<size_type> vector_size_type;
	typedef std::vector<scalar_type> vector_type;
	typedef std::vector<complex_scalar_type> complex_vector_type;
	typedef std::function<scalar_type(base_node)> function_type;
	typedef std::vector<std::function<scalar_type(base_node)>> vector_function_type;
	typedef std::vector<std::string> vector_string_type;

	//gmm namespace
	typedef gmm::rsvector<scalar_type> sparse_vector_type;
	typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
	typedef gmm::dense_matrix<scalar_type> dense_matrix_type;

	//getfem namespace
	//could be useful to avoid namespace getfem and putting namespace femg
}

#endif
