// Header with custom definitions. To be completed.

#ifndef FEMG_TYPE_ALIASES_HPP
#define FEMG_TYPE_ALIASES_HPP

namespace getfem {

	//bgeot namespace
	using bgeot::scalar_type;
	using bgeot::size_type;

	typedef std::vector<size_type> vector_size_type;
	typedef std::vector<scalar_type> vector_type;

	//gmm namespace
	typedef gmm::rsvector<scalar_type> sparse_vector_type;
	typedef gmm::row_matrix<scalar_type> sparse_matrix_type;
	typedef gmm::dense_matrix<scalar_type> dense_matrix_type;

	//getfem namespace
	//could be useful to avoid namespace getfem and putting namespace femg
}

#endif
