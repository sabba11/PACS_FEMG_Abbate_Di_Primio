// Source file containing definitions for class eigen_problem.

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
#include "elliptic_problem.hpp"
#include "mesh_1d.hpp"


namespace getfem {

	void
	elliptic_problem::init(int argc, char *argv[])
	{
	     //1. Read the .param filename from standard input
	     INPUT.read_command_line(argc, argv);

	     //2. Import data (algorithm specifications, boundary conditions,...)
	     import_data();

	     //3. Build mesh for the graph
	     build_mesh();

	     //4. Set finite elements and integration methods
	     set_im_and_fem();

	     //5. Build problem parameters
	     build_param();

	  	 //6. Set default values for coefficients
	  	 weights = vector_type(n_totalvert, 1.0);
	  	 potential = vector_type(n_totalvert, 0.0);

	     //7. Build the lists of the data of the vertices
	     // build_vertices_lists();
	}

	void
	elliptic_problem::import_data()
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "Importing descriptors for the problem..."<< std::endl;
	    #endif

	    descr.import_all(INPUT);

	    #ifdef FEMG_VERBOSE_
	    std::cout << descr;
	    #endif
	}

	void
	elliptic_problem::build_mesh()
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "Importing the mesh for the graph..."<< std::endl;
	    #endif

	    std::ifstream ifs(descr.MESH_FILEG);
	    GMM_ASSERT1(ifs.good(),"Unable to read from file " << descr.MESH_FILEG);
      //ACHTUNG!: for the moment we don't give mesh step as an input parameter
			// if we reimplement it back you should give it
	    import_pts_file(ifs, meshg, BCg, n_origvert, dim_prob, tg_vectors, n_vertices, descr.MESH_TYPEG);

	    n_branches = n_vertices.size();
  		n_totalvert = meshg.nb_points();

	    ifs.close();
	}

	void
	elliptic_problem::set_im_and_fem()
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "Setting Integration Methods for the discrete problem..."<< std::endl;
	    #endif

	    pintegration_method pim_g = int_method_descriptor(descr.IM_TYPEG);
	    mimg.set_integration_method(meshg.convex_index(), pim_g);

	    #ifdef FEMG_VERBOSE_
	    std::cout << "Setting Finite Element Methods for the discrete problem..."<< std::endl;
	    #endif

	    bgeot::pgeometric_trans pgt_g = bgeot::geometric_trans_descriptor(descr.MESH_TYPEG);

	    pfem pf_Ug = fem_descriptor(descr.FEM_TYPEG);
	    pfem pf_coeffg = fem_descriptor(descr.FEM_TYPEG_DATA);

	    mf_Ug.set_finite_element(meshg.convex_index(), pf_Ug);
	    mf_coeffg.set_finite_element(meshg.convex_index(), pf_coeffg);
	}

	void
	elliptic_problem::build_param(void)
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "Building parameters for the problem..." << std::endl;
	    #endif
      	if (descr.IMPORT_RADIUS) {
        	cout << "Importing radius values from file " << descr.RFILE << " ..." << endl;
        	std::ifstream ist(descr.RFILE);
        	if (!ist) cerr << "Impossible to read from file " << descr.RFILE << endl;
        	import_network_radius(radii, n_branches, ist, meshg, tg_vectors);
        	ist.close();
      	}
      return;
	}

  //Builds coefficient vectors
  void
  elliptic_problem::set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec, const unsigned & n_mean_points)
  {
	GMM_ASSERT1((f_vec.size() == s_vec.size()) && (f_vec.size() <= 2) && (!f_vec.empty()), "Input data has invalid size.");
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
		GMM_ASSERT1(s_vec[i] == "weight" || s_vec[i] == "potential", "Invalid string input data.");
		std::cout << "weights vector" << std::endl;
		for (unsigned l = 0; l < weight.size(); l++) {
			std::cout << weight[l] << " " << weight_den[l] << std::endl;
		}
		if (s_vec[i] == "weight")
			weight.swap(weights);
		else if (s_vec[i] == "potential")
			weight.swap(potential);
	}
	return;
  }


  void
  elliptic_problem::set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec)
  {
    GMM_ASSERT1((f_vec.size() == s_vec.size()) && (f_vec.size() <= 2) && (!f_vec.empty()), "Input data has invalid size.");
    for (unsigned i = 0; i < f_vec.size(); i++) {
    vector_type weight(n_totalvert);
      for (unsigned j = 0; j < meshg.points().size(); j++)
        weight[j] = f_vec[i](meshg.points()[j]);
      GMM_ASSERT1(s_vec[i] == "weight" || s_vec[i] == "potential", "Invalid string input data.");
      if (s_vec[i] == "weight")
        weight.swap(weights);
      else if (s_vec[i] == "potential")
        weight.swap(potential);
    }
    return;
  }


  //Builds source vector
  void
  elliptic_problem::set_source(const function_type & f, const unsigned & n_mean_points)
  {
		vector_size_type F_den(n_totalvert, 0);
    F_source.resize(n_totalvert);
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
				F_den[*it]++; //update the counter of ith point
				F_source[*it] += compute_circular_mean(n_mean_points, radii[k], meshg.points()[*it], tg_vectors[k], f); //update value
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
					F_den[pts[m]]++; //update counter
					F_source[pts[m]] += compute_circular_mean(n_mean_points, radii[k], meshg.points()[pts[m]], tg_vectors[k], f); //update value
					k++; //update counter
				}// if real point
			}
		}
		for (size_type n = 0; n < n_totalvert; n++)
			F_source[n] = F_source[n]/F_den[n]; //all points will be counted at least once so weight_den cannot be 0
	return;
  }

  void
  elliptic_problem::set_source(const function_type & f)
  {
	  F_source.resize(n_totalvert);
   	for (unsigned i = 0; i < meshg.points().size(); i++)
   	 	F_source[i] = f(meshg.points()[i]);
    return;
  }

  scalar_type
  elliptic_problem::compute_circular_mean(
    const unsigned & n_mean_points,
    const scalar_type & radius,
    const base_node & point,
    const vector_type & tg_vector,
    const function_type & f)
  {
    std::vector<base_node> mean_points(n_mean_points);
    for(size_type i=0; i<n_mean_points; i++){
	  mean_points[i].resize(3);
      mean_points[i][0] = radius*cos(2*M_PI*i/n_mean_points);
      mean_points[i][1] = 0.0;
      mean_points[i][2] = radius*sin(2*M_PI*i/n_mean_points);
    }
    dense_matrix_type Rotx_1(3,3);
    dense_matrix_type Rotx_2(3,3);
    dense_matrix_type Rotz_1(3,3);
    dense_matrix_type Rotz_2(3,3);

	//indexes must start from 0!!
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
	for(size_type i=0; i<n_mean_points; i++){
      v[0] = mean_points[i][0]; v[1] = mean_points[i][1]; v[2] = mean_points[i][2];
      if (tg_vector[1]>0){
        gmm::mult(Rotx_1,v,w);
        gmm::mult(Rotz_1,w,v);
      }
      else{
        gmm::mult(Rotx_2,v,w);
        gmm::mult(Rotz_2,w,v);
      }
      mean_points[i][0] = v[0] + point[0]; mean_points[i][1] =  v[1] + point[1]; mean_points[i][2] = v[2] + point[2];
    }
    scalar_type mean = 0.0;
    for(size_type i=0; i<n_mean_points; i++){
      mean += f(mean_points[i]);
    }
    mean = mean/n_mean_points;
    return mean;
  }


	void
	elliptic_problem::assembly(void)
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[elliptic_problem] Assembling matrices..." << std::endl;
		#endif

		assembly_matL();
		assembly_matA();

		return;
	}

	void
	elliptic_problem::assembly_matL()
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[elliptic_problem] Assembling stiffness matrix L..." << std::endl;
		#endif
		L.resize(n_totalvert, n_totalvert);
		getfem::asm_stiffness_matrix_for_laplacian(L, mimg, mf_Ug, mf_coeffg, weights);
		gmm::clean(L, 1E-10);
		return;
	}

	void
	elliptic_problem::assembly_matA()
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[elliptic_problem] Assembling operator matrix A..." << std::endl;
		#endif
		V.resize(n_totalvert, n_totalvert);
		getfem::asm_mass_matrix_param(V, mimg, mf_Ug, mf_coeffg, potential);
		A.resize(n_totalvert, n_totalvert);
		gmm::add(L, V, A);
  	return;
	}

	bool
	elliptic_problem::solve(void)
	{
		bool check_source = false;
		for (auto i:F_source){
			if (i!=0)
			 check_source = true;
		 }
		GMM_ASSERT1(check_source, "Invalid or empty source term, please set it again")
		if (descr.COMP_METHOD == "GMRES") {
			#ifdef FEMG_VERBOSE_
			std::cout << "[elliptic_problem] Solving the problem via GMRES algorithm..." << std::endl;
			#endif
			try {
				clock_t time_solve,time_tot;
				time_tot = clock();
				U.resize(n_totalvert,1);
				scalar_type tol = descr.TOL;
				std::cout << "[elliptic_problem] Starting GMRES routine..." << std::endl;
				gmm::iteration iter(tol);
				if (descr.BY_ITERATION)
					iter.set_maxiter(descr.ITER);
				size_t restart = 50;
				gmm::identity_matrix PR;
				time_solve = clock();

				gmm::gmres(A, U, F_source, PR, restart, iter);
				time_solve = clock() - time_solve;
				log_data.push_back(std::make_pair("Time to GMRES convergence (seconds): ", static_cast<float>(time_solve)/CLOCKS_PER_SEC));
				n_iteration = iter.get_iteration();
				converged_by_tol = iter.converged();
				time_tot = clock() - time_tot;
				log_data.push_back(std::make_pair("Time to compute solution (seconds): ", static_cast<float>(time_tot)/CLOCKS_PER_SEC));
				log_data.push_back(std::make_pair("Number of iteration: ", static_cast<float>(n_iteration)));
			  log_data.push_back(std::make_pair("Converged by tolerance condition: ", static_cast<float>(converged_by_tol)));
			}
			GMM_STANDARD_CATCH_ERROR;
		}
		else if (descr.COMP_METHOD == "QMR") {
			#ifdef FEMG_VERBOSE_
			std::cout << "[elliptic_problem] Solving the problem via QMR algorithm..." << std::endl;
			#endif
			try {
				clock_t time_solve,time_tot;
				time_tot = clock();
				U.resize(n_totalvert,1);
				scalar_type tol = descr.TOL;
				std::cout << "[elliptic_problem] Starting QMR routine..." << std::endl;
				gmm::iteration iter(tol);
				if (descr.BY_ITERATION)
					iter.set_maxiter(descr.ITER);
				gmm::identity_matrix PR;
				time_solve = clock();
				gmm::qmr(A, U, F_source, PR, iter);
				time_solve = clock() - time_solve;
				log_data.push_back(std::make_pair("Time to QMR convergence (seconds): ", static_cast<float>(time_solve)/CLOCKS_PER_SEC));
				n_iteration = iter.get_iteration();
				converged_by_tol = iter.converged();
				time_tot = clock() - time_tot;
				log_data.push_back(std::make_pair("Time to compute solution (seconds): ", static_cast<float>(time_tot)/CLOCKS_PER_SEC));
				log_data.push_back(std::make_pair("Number of iteration: ", static_cast<float>(n_iteration)));
				log_data.push_back(std::make_pair("Converged by tolerance condition: ", static_cast<float>(converged_by_tol)));
			}
			GMM_STANDARD_CATCH_ERROR;
		}
		else if (descr.COMP_METHOD == "CG") {
			#ifdef FEMG_VERBOSE_
			std::cout << "[elliptic_problem] Solving the problem via CG type algorithm..." << std::endl;
			#endif
			try {
				clock_t time_solve,time_tot;
				time_tot = clock();
				U.resize(n_totalvert,1);
				scalar_type tol = descr.TOL;
				std::cout << "[elliptic_problem] Starting CG type routine..." << std::endl;
				gmm::iteration iter(tol);
				if (descr.BY_ITERATION)
					iter.set_maxiter(descr.ITER);
				bool check_sym = false;
				dense_matrix_type A_trans(n_totalvert,n_totalvert);
				gmm::copy(A, A_trans);
        gmm::transposed(A_trans);
				dense_matrix_type A_sym(n_totalvert,n_totalvert);
				gmm::add(A, gmm::scaled(A_trans, -1.0), A_sym);
				if (gmm::mat_maxnorm(A_sym)<1E-12)
					check_sym = true;
				if (check_sym){
					#ifdef FEMG_VERBOSE_
					std::cout << "[elliptic_problem] Using Standard CG..." << std::endl;
					#endif
					gmm::add(A,A_trans,A_sym);
					gmm:scaled(A_sym,0.5);
					time_solve = clock();
					gmm::least_squares_cg(A_sym, U, F_source, iter);
					time_solve = clock() - time_solve;
					log_data.push_back(std::make_pair("Time to Standard CG convergence (seconds): ", static_cast<float>(time_solve)/CLOCKS_PER_SEC));
					n_iteration = iter.get_iteration();
					converged_by_tol = iter.converged();
					time_tot = clock() - time_tot;
					log_data.push_back(std::make_pair("Time to compute solution (seconds): ", static_cast<float>(time_tot)/CLOCKS_PER_SEC));
					log_data.push_back(std::make_pair("Number of iteration: ", static_cast<float>(n_iteration)));
					log_data.push_back(std::make_pair("Converged by tolerance condition: ", static_cast<float>(converged_by_tol)));
				} else{
					#ifdef FEMG_VERBOSE_
					std::cout << "[elliptic_problem] Using BICGSTAB..." << std::endl;
					#endif
			  	gmm::identity_matrix PR;
					time_solve = clock();
          gmm::bicgstab(A, U, F_source, PR, iter);
					time_solve = clock() - time_solve;
					log_data.push_back(std::make_pair("Time to BICGSTAB convergence (seconds): ", static_cast<float>(time_solve)/CLOCKS_PER_SEC));
					n_iteration = iter.get_iteration();
					converged_by_tol = iter.converged();
					time_tot = clock() - time_tot;
					log_data.push_back(std::make_pair("Time to compute solution (seconds): ", static_cast<float>(time_tot)/CLOCKS_PER_SEC));
	        log_data.push_back(std::make_pair("Number of iteration: ", static_cast<float>(n_iteration)));
					log_data.push_back(std::make_pair("Converged by tolerance condition: ", static_cast<float>(converged_by_tol)));
				}
			}
			GMM_STANDARD_CATCH_ERROR;
		}
		else if (descr.COMP_METHOD == "LU") {
			#ifdef FEMG_VERBOSE_
			std::cout << "[elliptic_problem] Solving the problem via SuperLU algorithm..." << std::endl;
			#endif
			try {
				clock_t time_solve,time_tot;
				time_tot = clock();
				U.resize(n_totalvert,1);
				std::cout << "[elliptic_problem] Starting SuperLU routine..." << std::endl;
				time_solve = clock();
				SuperLU_solve(A, U, F_source, condest);
				condest = 1/condest;
				time_solve = clock() - time_solve;
				log_data.push_back(std::make_pair("Time for SuperLU algorithm (seconds): ", static_cast<float>(time_solve)/CLOCKS_PER_SEC));
				#ifdef FEMG_VERBOSE_
				std::cout << "[elliptic_problem]Condition number: " << condest << std::endl;
				#endif
				time_tot = clock() - time_tot;
				log_data.push_back(std::make_pair("Time to compute solution (seconds): ", static_cast<float>(time_tot)/CLOCKS_PER_SEC));
				log_data.push_back(std::make_pair("Inverse conditioning number estimate: ", condest));
			}
			GMM_STANDARD_CATCH_ERROR;
		}
		else {
			std::cout << "[elliptic_problem] Invalid computational method descriptor." << std::endl;
		}
		return true;
	}

	void
	elliptic_problem::sol_export(void)
	{
	// Temporary code.
		std::cout << "[elliptic_problem] Printing solution..." << std::endl;
		for(auto it = U.begin(); it != U.end(); it++)
		   std::cout << *it << std::endl;
		// Number of vertices in the graph without the mesh
		std::cout << "Number of original vertices: " << n_origvert << std::endl;
		// Number of vertices in the extended graph
		std::cout << "Number of total vertices: " << n_totalvert << std::endl;

		// Exporting routine to MATLAB interface.
		// Executing this code requires GetFEM++-MATLAB interface and the boost::filesystem library.
		// sudo nano ~/.bashrc
		// ultima riga export LD_LIBRARY_PATH=path/to/boost/lib:$LD_LIBRARY_PATH
		std::string output = descr.OUTPUT;
		std::ostringstream dir_name_builder;
		dir_name_builder << output << "export/" << std::to_string(n_totalvert) + " point-mesh/" << descr.COMP_METHOD;
		boost::filesystem::path dir(dir_name_builder.str().c_str());
		boost::system::error_code ec;
		bool status = boost::filesystem::create_directories(dir, ec);
		if( status )
			std::cout << "[elliptic_problem] Directory " << dir_name_builder.str() << " created." << std::endl;
		else {
			if ( ec.message() == "Success" )
				std::cout << "[elliptic_problem] Directory " << dir_name_builder.str() << " already exists." << std::endl;
			else {
				std::cerr << "[elliptic_problem] Error while creating export directory: " << ec.message() << std::endl;
				return;
			}
		}
		#ifdef FEMG_VERBOSE_
		std::cout << "[elliptic_problem] Exporting solution files..." << std::endl;
		#endif
		std::ostringstream solution_name_builder, log_name_builder;
		solution_name_builder << dir_name_builder.str() << "/" << "solution.U";
		log_name_builder << dir_name_builder.str() << "/" << "export.log";
		std::fstream log(log_name_builder.str(), std::ios::out);
		if (!log) {
			std::cerr << "[elliptic_problem] Unable to open log file " << log_name_builder.str() << std::endl;
		}
		std::fstream sol(solution_name_builder.str(), std::ios::out);
		if (!sol) {
			log << "Exporting failed.\n";
			log << "Unable to open file " << solution_name_builder.str() << " to export eigenvalues.\n";
			return;
		}
		for (unsigned j = 0; j < gmm::vect_size(U); ++j)
				sol << U[j] << "\n";
		sol.close();
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
