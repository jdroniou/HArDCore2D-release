#ifndef PARALLEL_FOR
#define PARALLEL_FOR

#include <thread>
#include <Eigen/Sparse>

namespace HArDCore2D {
  
  /*!
   *	\addtogroup Common
   * @{
   */
   
  /// Generic function to execute threaded processes
  static inline void parallel_for(unsigned nb_elements,
                                  std::function<void(size_t start, size_t end)> functor,
                                  bool use_threads = true)
  {
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

    unsigned batch_size = nb_elements / nb_threads;
    unsigned batch_remainder = nb_elements % nb_threads;

    std::vector<std::thread> my_threads(nb_threads);

    if (use_threads) {
      // Multithread execution
      for (unsigned i = 0; i < nb_threads; ++i) {
        int start = i * batch_size;
        my_threads[i] = std::thread(functor, start, start + batch_size);
      }
    } else {
      // Single thread execution (for easy debugging)
      for(unsigned i = 0; i < nb_threads; ++i) {
        int start = i * batch_size;
        functor(start, start + batch_size);
      }
    }

    // Execute the elements left
    int start = nb_threads * batch_size;
    functor(start, start + batch_remainder);

    // Wait for the other thread to finish their task
    if (use_threads) {
      std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
    }
  }
  
  
  /// Function to assemble global matrix and right-hand side from a procedure that computes local triplets and rhs contributions
  static inline std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd>
  parallel_assembly_system(
                           size_t nb_elements, //< nb of elements over which the threading will be done
                           size_t size_system, //< size of the system to assemble
                           std::function<void(size_t start, size_t end, std::list<Eigen::Triplet<double>> * triplets, Eigen::VectorXd * rhs)> batch_local_assembly, //< procedure to compute all the local contributions to the matrix (put in the triplets) and rhs (put in the vector) between start and end points (e.g. elements indices) 
                           bool use_threads = true //< determine if threaded process is used or not
                           )
  {
    // Matrix and rhs
    Eigen::SparseMatrix<double> A(size_system, size_system);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(size_system);
  
    if (use_threads) {
      // Select the number of threads
      unsigned nb_threads_hint = std::thread::hardware_concurrency();
      unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

      // Compute the batch size and the remainder
      unsigned batch_size = nb_elements / nb_threads;
      unsigned batch_remainder = nb_elements % nb_threads;
      
      // Create vectors of triplets and vectors
      std::vector<std::list<Eigen::Triplet<double> > > triplets(nb_threads + 1);
      std::vector<Eigen::VectorXd> rhs(nb_threads + 1);

      for (unsigned i = 0; i < nb_threads + 1; i++) {
        rhs[i] = Eigen::VectorXd::Zero(size_system);
      } // for i

      // Assign a task to each thread
      std::vector<std::thread> my_threads(nb_threads);
      for (unsigned i = 0; i < nb_threads; ++i) {
        int start = i * batch_size;
        my_threads[i] = std::thread(batch_local_assembly, start, start + batch_size, &triplets[i], &rhs[i]);
      }

      // Execute the elements left
      int start = nb_threads * batch_size;
      batch_local_assembly(start, start + batch_remainder, &triplets[nb_threads], &rhs[nb_threads]);

      // Wait for the other threads to finish their task
      std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));

      // Create matrix from triplets
      size_t n_triplets = 0;
      for (auto triplets_thread : triplets) {
        n_triplets += triplets_thread.size();
      }
      std::vector<Eigen::Triplet<double> > all_triplets(n_triplets);
      auto triplet_index = all_triplets.begin();
      for (auto triplets_thread : triplets) {
        triplet_index = std::copy(triplets_thread.begin(), triplets_thread.end(), triplet_index);
      }
      A.setFromTriplets(all_triplets.begin(), all_triplets.end());
      for (auto rhs_thread : rhs) {
        b += rhs_thread;
      }
    } else {
      std::list<Eigen::Triplet<double> > triplets;
      batch_local_assembly(0, nb_elements, &triplets, &b);
      A.setFromTriplets(triplets.begin(), triplets.end());
    }

    return std::make_pair(A, b);  
  }
  
  /// Function to assemble a two global matrices and one vector (such as: system matrix+vector and matrix for BC) from a procedure that computes local triplets and rhs contributions
  static inline std::tuple<Eigen::SparseMatrix<double>, Eigen::VectorXd, Eigen::SparseMatrix<double>>
  parallel_assembly_system(
                           size_t nb_elements, //< nb of elements over which the threading will be done
                           size_t size_system1, //< size of the system 1 to assemble (must be square, corresponds to first matrix and vector=rhs)
                           std::pair<size_t, size_t> size_Mat2, //< size of the second matrix to assemble (can be rectangular)
                           std::function<void(size_t start, size_t end, std::list<Eigen::Triplet<double>> * triplets1, Eigen::VectorXd * vec1, std::list<Eigen::Triplet<double>> * triplets2)> batch_local_assembly, //< procedure to compute all the local contributions to the matrices (put in the triplets) and vectors between start and end points (e.g. elements indices) 
                           bool use_threads = true //< determine if threaded process is used or not
                           )
  {
    // Matrices and vectors
    Eigen::SparseMatrix<double> A1(size_system1, size_system1);
    Eigen::VectorXd b1 = Eigen::VectorXd::Zero(size_system1);
    Eigen::SparseMatrix<double> A2(size_Mat2.first, size_Mat2.second);
  
    if (use_threads) {
      // Select the number of threads
      unsigned nb_threads_hint = std::thread::hardware_concurrency();
      unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

      // Compute the batch size and the remainder
      unsigned batch_size = nb_elements / nb_threads;
      unsigned batch_remainder = nb_elements % nb_threads;
      
      // Create vectors of triplets and vectors
      std::vector<std::list<Eigen::Triplet<double> > > triplets1(nb_threads + 1);
      std::vector<Eigen::VectorXd> vec1(nb_threads + 1);
      std::vector<std::list<Eigen::Triplet<double> > > triplets2(nb_threads + 1);

      for (unsigned i = 0; i < nb_threads + 1; i++) {
        vec1[i] = Eigen::VectorXd::Zero(size_system1);
      } // for i

      // Assign a task to each thread
      std::vector<std::thread> my_threads(nb_threads);
      for (unsigned i = 0; i < nb_threads; ++i) {
        int start = i * batch_size;
        my_threads[i] = std::thread(batch_local_assembly, start, start + batch_size, &triplets1[i], &vec1[i], &triplets2[i]);
      }

      // Execute the elements left
      int start = nb_threads * batch_size;
      batch_local_assembly(start, start + batch_remainder, &triplets1[nb_threads], &vec1[nb_threads], &triplets2[nb_threads]);

      // Wait for the other threads to finish their task
      std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));

      // Create system 1 from triplets
      size_t n_triplets1 = 0;
      for (auto triplets_thread : triplets1) {
        n_triplets1 += triplets_thread.size();
      }
      std::vector<Eigen::Triplet<double> > all_triplets1(n_triplets1);
      auto triplet1_index = all_triplets1.begin();
      for (auto triplets_thread : triplets1) {
        triplet1_index = std::copy(triplets_thread.begin(), triplets_thread.end(), triplet1_index);
      }
      A1.setFromTriplets(all_triplets1.begin(), all_triplets1.end());
      for (auto vec_thread : vec1) {
        b1 += vec_thread;
      }

      // Create system 2 from triplets
      size_t n_triplets2 = 0;
      for (auto triplets_thread : triplets2) {
        n_triplets2 += triplets_thread.size();
      }
      std::vector<Eigen::Triplet<double> > all_triplets2(n_triplets2);
      auto triplet2_index = all_triplets2.begin();
      for (auto triplets_thread : triplets2) {
        triplet2_index = std::copy(triplets_thread.begin(), triplets_thread.end(), triplet2_index);
      }
      A2.setFromTriplets(all_triplets2.begin(), all_triplets2.end());


    } else {
      std::list<Eigen::Triplet<double> > triplets1;
      std::list<Eigen::Triplet<double> > triplets2;
      batch_local_assembly(0, nb_elements, &triplets1, &b1, &triplets2);
      A1.setFromTriplets(triplets1.begin(), triplets1.end());
      A2.setFromTriplets(triplets2.begin(), triplets2.end());      
    }

    return std::make_tuple(A1, b1, A2);  
  }


  /// Function to assemble a two global matrices and vectors (such as: system and static condensation operator) from a procedure that computes local triplets and rhs contributions
  static inline std::tuple<Eigen::SparseMatrix<double>, Eigen::VectorXd, Eigen::SparseMatrix<double>, Eigen::VectorXd>
  parallel_assembly_system(
                           size_t nb_elements, //< nb of elements over which the threading will be done
                           size_t size_system1, //< size of the system 1 to assemble (must be square, corresponds to first matrix and vector=rhs)
                           std::pair<size_t, size_t> size_Mat2, //< size of the second matrix to assemble (can be rectangular)
                           size_t size_b2, //< size of second vector
                           std::function<void(size_t start, size_t end, std::list<Eigen::Triplet<double>> * triplets1, Eigen::VectorXd * vec1, std::list<Eigen::Triplet<double>> * triplets2, Eigen::VectorXd * vec2)> batch_local_assembly, //< procedure to compute all the local contributions to the matrices (put in the triplets) and vectors between start and end points (e.g. elements indices) 
                           bool use_threads = true //< determine if threaded process is used or not
                           )
  {
    // Matrices and vectors
    Eigen::SparseMatrix<double> A1(size_system1, size_system1);
    Eigen::VectorXd b1 = Eigen::VectorXd::Zero(size_system1);
    Eigen::SparseMatrix<double> A2(size_Mat2.first, size_Mat2.second);
    Eigen::VectorXd b2 = Eigen::VectorXd::Zero(size_b2);
  
    if (use_threads) {
      // Select the number of threads
      unsigned nb_threads_hint = std::thread::hardware_concurrency();
      unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

      // Compute the batch size and the remainder
      unsigned batch_size = nb_elements / nb_threads;
      unsigned batch_remainder = nb_elements % nb_threads;
      
      // Create vectors of triplets and vectors
      std::vector<std::list<Eigen::Triplet<double> > > triplets1(nb_threads + 1);
      std::vector<Eigen::VectorXd> vec1(nb_threads + 1);
      std::vector<std::list<Eigen::Triplet<double> > > triplets2(nb_threads + 1);
      std::vector<Eigen::VectorXd> vec2(nb_threads + 1);

      for (unsigned i = 0; i < nb_threads + 1; i++) {
        vec1[i] = Eigen::VectorXd::Zero(size_system1);
        vec2[i] = Eigen::VectorXd::Zero(size_b2);
      } // for i

      // Assign a task to each thread
      std::vector<std::thread> my_threads(nb_threads);
      for (unsigned i = 0; i < nb_threads; ++i) {
        int start = i * batch_size;
        my_threads[i] = std::thread(batch_local_assembly, start, start + batch_size, &triplets1[i], &vec1[i], &triplets2[i], &vec2[i]);
      }

      // Execute the elements left
      int start = nb_threads * batch_size;
      batch_local_assembly(start, start + batch_remainder, &triplets1[nb_threads], &vec1[nb_threads], &triplets2[nb_threads], &vec2[nb_threads]);

      // Wait for the other threads to finish their task
      std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));

      // Create system 1 from triplets
      size_t n_triplets1 = 0;
      for (auto triplets_thread : triplets1) {
        n_triplets1 += triplets_thread.size();
      }
      std::vector<Eigen::Triplet<double> > all_triplets1(n_triplets1);
      auto triplet1_index = all_triplets1.begin();
      for (auto triplets_thread : triplets1) {
        triplet1_index = std::copy(triplets_thread.begin(), triplets_thread.end(), triplet1_index);
      }
      A1.setFromTriplets(all_triplets1.begin(), all_triplets1.end());
      for (auto vec_thread : vec1) {
        b1 += vec_thread;
      }

      // Create system 2 from triplets
      size_t n_triplets2 = 0;
      for (auto triplets_thread : triplets2) {
        n_triplets2 += triplets_thread.size();
      }
      std::vector<Eigen::Triplet<double> > all_triplets2(n_triplets2);
      auto triplet2_index = all_triplets2.begin();
      for (auto triplets_thread : triplets2) {
        triplet2_index = std::copy(triplets_thread.begin(), triplets_thread.end(), triplet2_index);
      }
      A2.setFromTriplets(all_triplets2.begin(), all_triplets2.end());
      for (auto vec_thread : vec2) {
        b2 += vec_thread;
      }


    } else {
      std::list<Eigen::Triplet<double> > triplets1;
      std::list<Eigen::Triplet<double> > triplets2;
      batch_local_assembly(0, nb_elements, &triplets1, &b1, &triplets2, &b2);
      A1.setFromTriplets(triplets1.begin(), triplets1.end());
      A2.setFromTriplets(triplets2.begin(), triplets2.end());      
    }

    return std::make_tuple(A1, b1, A2, b2);  
  }

  //@}
  
} // end of namespace HArDCore2D
#endif
