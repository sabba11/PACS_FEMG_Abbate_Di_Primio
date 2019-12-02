/*!
	\file eigen_problem.cpp
	\author Stefano Abbate
	\author Andrea Di Primio
	\brief Header containing the elliptic_problem derived class.
*/

// GetFEM++ libraries
#include <getfem/getfem_assembling.h>
#include <gmm/gmm_condition_number.h>
#include <gmm/gmm_except.h>
#include <gmm/gmm_iter_solvers.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm.h>

// Standard libraries
#include <algorithm>
#include <map>
#include <complex>
#include <iostream>
#include <fstream>
#include <time.h>
#include <boost/filesystem.hpp>
#include <vector>
#include <math.h>

// Project headers
#include "eigen_problem.hpp"

extern "C" {void dggev_( char* jobvl, char* jobvr, long unsigned int* na, double* a, long unsigned int* nb, double* b,
                         long unsigned int* lda, double* wr, double* wi, double* wd, double* vl, long unsigned int* ldvl,
                         double* vr, long unsigned int* ldvr, double* work, long int* lwork, long int* info );}
// extern void dggev_( char* jobvl, char* jobvr, long unsigned int* na, double* a, long unsigned int* nb, double* b,
//                 long unsigned int* lda, double* wr, double* wi, double* wd, double* vl, long unsigned int* ldvl,
//                 double* vr, long unsigned int* ldvr, double* work, long int* lwork, long int* info );

namespace getfem {

	void
	eigen_problem::import_data()
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "[eigen_problem] Importing problem descriptors..." << std::endl;
	    #endif

	    descr.import_all(INPUT);

	    #ifdef FEMG_VERBOSE_
	    std::cout << descr;
	    #endif
	}

	void
	eigen_problem::build_mesh()
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "[eigen_problem] Importing mesh from data file..." << std::endl;
	    #endif

	    std::ifstream ifs(descr.MESH_FILEG);
	    GMM_ASSERT1(ifs.good(),"Unable to read from file " << descr.MESH_FILEG);
      	//ACHTUNG!: for the moment we don't give mesh step as an input parameter
			// if we reimplement it back you should give it
		std::ifstream rad;
		if (descr.IMPORT_RADIUS) {
			#ifdef FEMG_VERBOSE_
			std::cout << "[eigen_problem] Importing radii from data file... " << std::endl;
			#endif

			rad.open(descr.RFILE);
			GMM_ASSERT1(rad.good(), "Unable to read from file " << descr.RFILE);
		}
	    import_pts_file(ifs, rad, descr.IMPORT_RADIUS);

  		n_totalvert = meshg.nb_points();

	    ifs.close();
		if (descr.IMPORT_RADIUS)
			rad.close();
		return;
	}

	void
	eigen_problem::set_im_and_fem()
	{
		#ifdef FEMG_VERBOSE_
	    std::cout << "[eigen_problem] Setting GetFEM++ integration method..."<< std::endl;
	    #endif

	    pintegration_method pim_g = int_method_descriptor(descr.IM_TYPEG);
	    mimg.set_integration_method(meshg.convex_index(), pim_g);

	    #ifdef FEMG_VERBOSE_
	    std::cout << "[eigen_problem] Setting GetFEM++ Finite Element method..." << std::endl;
	    #endif

	    bgeot::pgeometric_trans pgt_g = bgeot::geometric_trans_descriptor(descr.MESH_TYPEG);

	    pfem pf_Ug = fem_descriptor(descr.FEM_TYPEG);
	    pfem pf_coeffg = fem_descriptor(descr.FEM_TYPEG_DATA);

	    mf_Ug.set_finite_element(meshg.convex_index(), pf_Ug);
	    mf_coeffg.set_finite_element(meshg.convex_index(), pf_coeffg);
		return;
	}

	void
	eigen_problem::set_default_coefficients(void)
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Setting coefficients at default values..." << std::endl;
		#endif

		left_weights = vector_type(n_totalvert, 1.0);
		right_weights = vector_type(n_totalvert, 1.0);
		potential = vector_type(n_totalvert, 0.0);
		return;
	}

  	//Builds coefficient vectors
	void
	eigen_problem::set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec, const unsigned & n_mean_points)
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Computing coefficients on mesh nodes (via means)..." << std::endl;
		#endif

		GMM_ASSERT1((f_vec.size() == s_vec.size()) && (f_vec.size() <= 3) && (!f_vec.empty()), "Input data has invalid size.");
		for (unsigned i = 0; i < f_vec.size(); i++) {
			vector_size_type weight_den(n_totalvert, 0);
	  		vector_type weight(n_totalvert, 0.0);
			std::set<unsigned> branch_idx;
			unsigned k = 0; //aux counter
			unsigned thresh = n_totalvert - n_origvert; //distinguish between idxs of real vertices and nodes
			for (size_type b=0; b<n_branches; ++b) {
				for (getfem::mr_visitor mrv(meshg.region(b)); !mrv.finished(); ++mrv) {
					unsigned idx = mrv.cv(); //get convex index
					vector_size_type pts = meshg.ind_points_of_convex(idx); //get point indexes of convex idx
					for (size_type m = 0; m < pts.size(); m++)
						branch_idx.insert(pts[m]); //collect in set to eliminate duplicates
				}
				for (auto it = branch_idx.begin(); it != branch_idx.end(); it++) {
					weight_den[*it]++; //update the counter of ith point
					weight[*it] += compute_circular_mean(n_mean_points, radii[k], meshg.points()[*it], tg_vectors[k], f_vec[i]); //update value
					k++; //update the aux counter
				}
				branch_idx.clear(); //clear the set
			}
			// Check last region with real vertices
			for (getfem::mr_visitor mrv(meshg.region(n_branches)); !mrv.finished(); ++mrv) {
				unsigned idx = mrv.cv(); //get convex index
				vector_size_type pts = meshg.ind_points_of_convex(idx); //get point indexes of convex idx
				for (size_type m = 0; m < pts.size(); m++) {
					if (pts[m] >= thresh) {
						weight_den[pts[m]]++; //update counter
						weight[pts[m]] += compute_circular_mean(n_mean_points, radii[k], meshg.points()[pts[m]], tg_vectors[k], f_vec[i]); //update value
						k++; //update counter
					}// if real point
				}
			}
			for (size_type n = 0; n < n_totalvert; n++)
				weight[n] = weight[n]/weight_den[n]; //all points will be counted at least once so weight_den cannot be 0
			GMM_ASSERT1(s_vec[i] == "left" || s_vec[i] == "right" || s_vec[i] == "potential", "Invalid string input data.");
			if (s_vec[i] == "left")
				weight.swap(left_weights);
			else if (s_vec[i] == "right")
				weight.swap(right_weights);
			else if (s_vec[i] == "potential")
				weight.swap(potential);
		}
		return;
	}

	void
	eigen_problem::set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec)
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Computing coefficients on mesh nodes (via direct evaluation)..." << std::endl;
		#endif

		GMM_ASSERT1((f_vec.size() == s_vec.size()) && (f_vec.size() <= 3) && (!f_vec.empty()), "Input data has invalid size.");
		for (unsigned i = 0; i < f_vec.size(); i++) {
			vector_type weight(n_totalvert);
			for (unsigned j = 0; j < meshg.points().size(); j++)
				weight[j] = f_vec[i](meshg.points()[j]);
		 	GMM_ASSERT1(s_vec[i] == "left" || s_vec[i] == "right" || s_vec[i] == "potential", "Invalid string input data.");
		 	if (s_vec[i] == "left")
				weight.swap(left_weights);
		  	else if (s_vec[i] == "right")
				weight.swap(right_weights);
		  	else if (s_vec[i] == "potential")
				weight.swap(potential);
		}
		return;
	}

	scalar_type
	eigen_problem::compute_circular_mean
		(
		const unsigned & n_mean_points,
		const scalar_type & radius,
		const base_node & point,
		const vector_type & tg_vector,
		const function_type & f
		)
	{
		std::vector<base_node> mean_points(n_mean_points);
		for (size_type i = 0; i < n_mean_points; i++){
	  		mean_points[i].resize(3);
	  		mean_points[i][0] = radius*cos(2*M_PI*i/n_mean_points);
	  		mean_points[i][1] = 0.0;
	  		mean_points[i][2] = radius*sin(2*M_PI*i/n_mean_points);
		}
		dense_matrix_type Rotx_1(3,3);
		dense_matrix_type Rotx_2(3,3);
		dense_matrix_type Rotz_1(3,3);
		dense_matrix_type Rotz_2(3,3);

		Rotx_1(0,0) = 1;
		Rotx_1(0,1) = 0;
		Rotx_1(0,2) = 0;
		Rotx_1(1,0) = 0;
		Rotx_1(1,1) = sqrt(1-tg_vector[2]*tg_vector[2]);
		Rotx_1(1,2) = -tg_vector[2];
		Rotx_1(2,0) = 0;
		Rotx_1(2,1) = tg_vector[2];
		Rotx_1(2,2) = sqrt(1-tg_vector[2]*tg_vector[2]);

		Rotx_2(0,0) = 1;
		Rotx_2(0,1) = 0;
		Rotx_2(0,2) = 0;
		Rotx_2(1,0) = 0;
		Rotx_2(1,1) = -sqrt(1-tg_vector[2]*tg_vector[2]);
		Rotx_2(1,2) = -tg_vector[2];
		Rotx_2(2,0) = 0;
		Rotx_2(2,1) = tg_vector[2];
		Rotx_2(2,2) = -sqrt(1-tg_vector[2]*tg_vector[2]);

		Rotz_1(0,0) = sqrt(1-tg_vector[0]*tg_vector[0]);
		Rotz_1(0,1) = -tg_vector[0];
		Rotz_1(0,2) = 0;
		Rotz_1(1,0) = tg_vector[0];
		Rotz_1(1,1) = sqrt(1-tg_vector[0]*tg_vector[0]);
		Rotz_1(1,2) = 0;
		Rotz_1(2,0) = 0;
		Rotz_1(2,1) = 0;
		Rotz_1(2,2) = 1;

		Rotz_2(0,0) = -sqrt(1-tg_vector[0]*tg_vector[0]);
		Rotz_2(0,1) = -tg_vector[0];
		Rotz_2(0,2) = 0;
		Rotz_2(1,0) = tg_vector[0];
		Rotz_2(1,1) = -sqrt(1-tg_vector[0]*tg_vector[0]);
		Rotz_2(1,2) = 0;
		Rotz_2(2,0) = 0;
		Rotz_2(2,1) = 0;
		Rotz_2(2,2) = 1;
		vector_type v,w;
		v.resize(3); w.resize(3);
		for (size_type i = 0; i < n_mean_points; i++){
			v[0] = mean_points[i][0]; v[1] = mean_points[i][1]; v[2] = mean_points[i][2];
		  	if (tg_vector[1]>0) {
		    	gmm::mult(Rotx_1,v,w);
		    	gmm::mult(Rotz_1,w,v);
		  	}
		  	else {
		    gmm::mult(Rotx_2,v,w);
		    gmm::mult(Rotz_2,w,v);
		  	}
		  	mean_points[i][0] = v[0] + point[0]; mean_points[i][1] =  v[1] + point[1]; mean_points[i][2] = v[2] + point[2];
		}
		scalar_type mean = 0.0;
		for (size_type i = 0; i < n_mean_points; i++){
		  	mean += f(mean_points[i]);
		}
		mean = mean/n_mean_points;
		return mean;
	}

	void
	eigen_problem::assembly(void)
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Assembling matrices..." << std::endl;
		#endif

		assembly_matL();
		assembly_matM();
		assembly_matA();

		return;
	}

	void
	eigen_problem::assembly_matL()
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Assembling stiffness matrix L..." << std::endl;
		#endif

		L.resize(n_totalvert, n_totalvert);
		getfem::asm_stiffness_matrix_for_laplacian(L, mimg, mf_Ug, mf_coeffg, left_weights);
		gmm::clean(L, 1E-10);

		return;
	}

	void
	eigen_problem::assembly_matM()
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Assembling mass matrix M..." << std::endl;
		#endif

		M.resize(n_totalvert, n_totalvert);
		getfem::asm_mass_matrix_param(M, mimg, mf_Ug, mf_coeffg, right_weights);
		gmm::clean(M, 1E-10);

		return;
	}

	void
	eigen_problem::assembly_matA()
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Assembling operator matrix A..." << std::endl;
		#endif

		V.resize(n_totalvert, n_totalvert);
		getfem::asm_mass_matrix_param(V, mimg, mf_Ug, mf_coeffg, potential);
		A.resize(n_totalvert, n_totalvert);
		gmm::add(L, V, A);
		//vector_type R_dir(n_totalvert, 0);
		//vector_type zeros(n_totalvert, 0);
		//for (unsigned i = 0; i < BCg.size(); i++)
		//	if (BCg[i].label == "DIR")
		//		R_dir[BCg[i].idx] = BCg[i].value;
		//getfem::assembling_Dirichlet_condition(A, zeros	, mf_Ug, n_branches + 2, R_dir);
		return;
	}

	bool
	eigen_problem::solve(void)
	{
		if (descr.COMP_METHOD == "QR") {
			#ifdef FEMG_VERBOSE_
			std::cout << "[eigen_problem] Solving the problem via QR algorithm..." << std::endl;
			#endif
			try {
				clock_t time_inverse, time_eig, time_tot, time_mult;
				time_tot = clock();
				dense_matrix_type eigvects(n_totalvert, n_totalvert);
				// Initial guess for all eigenvalues is 1.
				complex_vector_type eigvals(n_totalvert, 1);
				dense_matrix_type inverse_mass(n_totalvert, n_totalvert);
				dense_matrix_type aux_H(n_totalvert, n_totalvert);
				gmm::copy(M, inverse_mass);
				gmm::copy(A, aux_H);
				scalar_type cond_number = gmm::condition_number(M);
				std::cout << "[eigen_problem] The mass matrix M has condition number " << cond_number << std::endl;
				time_inverse = clock();
				gmm::lu_inverse(inverse_mass);
				time_inverse = clock() - time_inverse;
				log_data.push_back(std::make_pair("Time for matrix inversion (seconds): ", static_cast<float>(time_inverse)/CLOCKS_PER_SEC));
        		time_mult = clock();
				gmm::mult(inverse_mass, aux_H, A);
				time_mult = clock() - time_mult;
				log_data.push_back(std::make_pair("Time for matrix multiplication (seconds): ", static_cast<float>(time_mult)/CLOCKS_PER_SEC));
				std::cout << static_cast<float>(time_inverse)/CLOCKS_PER_SEC << std::endl;
				scalar_type tol = descr.TOL;
				std::cout << "[eigen_problem] Starting QR routine..." << std::endl;
				dense_matrix_type dense_A(n_totalvert, n_totalvert);
				gmm::copy(A, dense_A);
				time_eig = clock();
				gmm::geev_interface_right(dense_A, eigvals, eigvects);
				time_eig = clock() - time_eig;
				log_data.push_back(std::make_pair("Time to QR convergence (seconds): ", static_cast<float>(time_eig)/CLOCKS_PER_SEC));
				#ifdef FEMG_VERBOSE_
				for (unsigned i = 0; i < log_data.size(); i++)
					std::cout << "[eigen_problem] " << log_data[i].first << log_data[i].second << std::endl;
				#endif
				for (unsigned i = 0; i < n_totalvert; i++) {
					vector_type aux_v;
					for (unsigned j = 0; j < n_totalvert; j++)
						aux_v.push_back(eigvects(j, i));
					if (abs(eigvals[i].imag()) < 1e-10) {
						auto eigpair = std::make_pair(eigvals[i].real(), aux_v);
						eigpairs.insert(eigpair);
					}
					else
						std::cout << "[eigen_problem] Warning: complex eigenvalue ignored. Value: " << eigvals[i].real() << " + " << eigvals[i].imag() << "i" << std::endl;
				}
				time_tot = clock() - time_tot;
				log_data.push_back(std::make_pair("Time to compute eigencouples (seconds): ", static_cast<float>(time_tot)/CLOCKS_PER_SEC));
			}
			GMM_STANDARD_CATCH_ERROR;
		}
		else if (descr.COMP_METHOD == "QZ") {
			#ifdef FEMG_VERBOSE_
			std::cout << "[eigen_problem] Solving the problem via QZ algorithm..." << std::endl;
			#endif
			try {
				clock_t time_eig,time_tot;
				time_tot = clock();
				dense_matrix_type eigvects(n_totalvert, n_totalvert);
				// Initial guess for all eigenvalues is 1.
				complex_vector_type eigvals(n_totalvert, 1);
				scalar_type tol = descr.TOL;
				std::cout << "[eigen_problem] Starting QZ routine..." << std::endl;
        		long info(0), lwork(-1);
				double work1;
        		gmm::dense_matrix<double> A_blas(n_totalvert, n_totalvert);
				gmm::copy(A, A_blas);
				gmm::dense_matrix<double> M_blas(n_totalvert, n_totalvert);
				gmm::copy(M, M_blas);
        		std::vector<double> eigvr(n_totalvert), eigvi(n_totalvert), eigvden(n_totalvert);
				time_eig = clock();
				char boolleft = 'V';
				char boolright = 'N';
        		dggev_(&boolleft, &boolright, &n_totalvert, &A_blas(0,0), &n_totalvert, &M_blas(0,0),
				      &n_totalvert, &eigvr[0], &eigvi[0], &eigvden[0],
							&eigvects(0,0), &n_totalvert, &eigvects(0,0), &n_totalvert, &work1, &lwork, &info);
        		lwork = long(gmm::real(work1));
        		std::vector<double> work(lwork);
				dggev_(&boolleft, &boolright, &n_totalvert, &A_blas(0,0), &n_totalvert, &M_blas(0,0),
				      &n_totalvert, &eigvr[0], &eigvi[0], &eigvden[0],
							&eigvects(0,0), &n_totalvert, &eigvects(0,0), &n_totalvert, &work[0], &lwork, &info);
				GMM_ASSERT1(!info, "QZ algorithm failed");

				for (unsigned i=0; i<n_totalvert; i++){
          			eigvi[i] = eigvi[i]/eigvden[i];
          			eigvr[i] = eigvr[i]/eigvden[i];
				}
        		gmm::copy(eigvr, gmm::real_part(const_cast<complex_vector_type &>(eigvals)));
        		gmm::copy(eigvi, gmm::imag_part(const_cast<complex_vector_type &>(eigvals)));
				time_eig = clock() - time_eig;
				log_data.push_back(std::make_pair("Time to QZ convergence (seconds): ", static_cast<float>(time_eig)/CLOCKS_PER_SEC));

				#ifdef FEMG_VERBOSE_
				for (unsigned i = 0; i < log_data.size(); i++)
					std::cout << "[eigen_problem] " << log_data[i].first << log_data[i].second << std::endl;
				#endif
				for (unsigned i = 0; i < n_totalvert; i++) {
					vector_type aux_v;
					for (unsigned j = 0; j < n_totalvert; j++)
						aux_v.push_back(eigvects(j, i));
					if (abs(eigvals[i].imag()) < 1e-10) {
						auto eigpair = std::make_pair(eigvals[i].real(), aux_v);
						eigpairs.insert(eigpair);
					}
					else
						std::cout << "[eigen_problem] Warning: complex eigenvalue ignored. Value: " << eigvals[i].real() << " + " << eigvals[i].imag() << "i" << std::endl;
				}
				time_tot = clock() - time_tot;
				log_data.push_back(std::make_pair("Time to compute eigencouples (seconds): ", static_cast<float>(time_tot)/CLOCKS_PER_SEC));
			}
			GMM_STANDARD_CATCH_ERROR;
		}
		else {
			std::cout << "[eigen_problem] Invalid computational method descriptor." << std::endl;
			return false;
		}
		return true;
	}

	void
	eigen_problem::sol_export(void)
	{
	// Temporary code.
		std::cout << "[eigen_problem] Printing eigenvalues..." << std::endl;
		for(auto it = eigpairs.begin(); it != eigpairs.end(); it++)
		   std::cout << it->first << std::endl;
		// Number of vertices in the graph without the mesh
		std::cout << "Number of original vertices: " << n_origvert << std::endl;
		// Number of vertices in the extended graph
		std::cout << "Number of total vertices: " << n_totalvert << std::endl;

		// Exporting routine to MATLAB interface.
		// Executing this code requires GetFEM++-MATLAB interface and the boost::filesystem library.
		// sudo nano ~/.bashrc
		// ultima riga export LD_LIBRARY_PATH=path/to/boost/lib:$LD_LIBRARY_PATH
		std::string output = descr.OUTPUT;
		unsigned n_digits = std::to_string(n_totalvert).size();
		unsigned i = 1;
		std::ostringstream dir_name_builder;
		dir_name_builder << output << "export/" << std::to_string(n_totalvert) + " point-mesh/" << descr.COMP_METHOD;
		boost::filesystem::path dir(dir_name_builder.str().c_str());
		boost::system::error_code ec;
		bool status = boost::filesystem::create_directories(dir, ec);
		if( status )
			std::cout << "[eigen_problem] Directory " << dir_name_builder.str() << " created." << std::endl;
		else {
			if ( ec.message() == "Success" )
				std::cout << "[eigen_problem] Directory " << dir_name_builder.str() << " already exists." << std::endl;
			else {
				std::cerr << "[eigen_problem] Error while creating export directory: " << ec.message() << std::endl;
				return;
			}
		}
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Exporting eigenvalues and eigenvector files..." << std::endl;
		#endif
		std::ostringstream eigval_name_builder, log_name_builder;
		std::string s(n_digits, '0');
		eigval_name_builder << dir_name_builder.str() << "/" << s << "-eigenvalues.U";
		log_name_builder << dir_name_builder.str() << "/" << "export.log";
		std::fstream log(log_name_builder.str(), std::ios::out);
		if (!log) {
			std::cerr << "[eigen_problem] Unable to open log file " << log_name_builder.str() << std::endl;
		}
		std::fstream eig(eigval_name_builder.str(), std::ios::out);
		if (!eig) {
			log << "Exporting failed.\n";
			log << "Unable to open file " << eigval_name_builder.str() << " to export eigenvalues.\n";
			return;
		}
		for (auto it = eigpairs.begin(); it != eigpairs.end(); it++) {
			unsigned n_zeros = n_digits - std::to_string(i).size();
			std::string zeros(n_zeros, '0');
			scalar_type eigen = it->first;
			std::ostringstream file_name_builder;
	 		file_name_builder << dir_name_builder.str() << "/" << zeros << i << "-eigenvector-" << eigen << ".U";
	 		std::fstream f(file_name_builder.str(), std::ios::out);
			if (!f) {
				log << "Exporting failed.\n";
				log << "Unable to open file " << file_name_builder.str() << " to export eigenvector " << i << "\n";
				return;
			}
	 		for (unsigned j = 0; j < gmm::vect_size(it->second); ++j)
	   			f << it->second[j] << "\n";
			eig << eigen << "\n";
			i++;
			f.close();
		}
		eig.close();
		// Exports log file if files are created successfully.
		log << "Exporting successful.\n";
		for (auto it = log_data.begin(); it != log_data.end(); it++) {
			log << it->first << it->second << "\n";
		}
		// Exports mesh and mesh_fem objects to text files.
		meshg.write_to_file(dir_name_builder.str() + "/mesh.mh" );
	 	// when the 2nd arg is true, the mesh is saved with the |mf|
		mf_Ug.write_to_file(dir_name_builder.str() + "/solution.mf", true);
		return;
	}
}// end of namespace
