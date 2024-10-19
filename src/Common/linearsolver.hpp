// Class to select and apply solver
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)


#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef WITH_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#include <mkl.h>
#endif

#ifdef WITH_PASTIX
#include <PastixInterface.hpp>
#endif

#include <boost/algorithm/string.hpp>

namespace HArDCore2D
{
  /*!
   *	\addtogroup Common
   * @{
   */
  
  /// Enumeration of all available solvers
  enum SolverName { EigenLU, EigenBiCGStab, PardisoLU, UMFPACK, PaStiXLU, PaStiXLLT };

  /// To parse optional Solver parameter from command line
  std::map<std::string, int> ParamName{
  #ifdef WITH_PASTIX
  {"eps_refinement",0},{"eps_ctrl",1}
  #endif
  };

  /// Create space and stores default values of optional solver parameters
  std::vector<double> SolverParam{
  #ifdef WITH_PASTIX
    -1.0, -1.0 // eps_refinement, eps_ctrl with default values
  #endif
  };

  /// Create a solver specific description in boost::desc.options
  const char *SolverParamHelper{
  #ifdef WITH_PASTIX
  "Set the solver parameters; Usage [with v1, v2 as doubles]: \n"
   "--solver-parameters param1=v1 param2=v2\n"
   "--solver-parameters param1=v1 --solver-parameters param2=v2 \n"
   "-sparm param=v \n"
   "Possible parameters are [all optionals]:\n eps_refinement \n eps_ctrl\n"
  #else
  "Set the solver parameters, not implemented with this solver, see src/Common/linearsolver.hpp for modifications\n"
  #endif
  };

  // Types of Eigen wrappers for each solver (if the module is not available, we define bools just to have a name of that type)
  template<typename MatrixType>
  using EigenLUType = Eigen::SparseLU<MatrixType>;
  template<typename MatrixType>
  using EigenBiCGStabType = Eigen::BiCGSTAB<MatrixType, Eigen::IncompleteLUT<double> >;

  template<typename MatrixType>
  #ifdef WITH_MKL
  using PardisoLUType = Eigen::PardisoLU<MatrixType>;
  #else
  using PardisoLUType = bool;
  #endif
  
  template<typename MatrixType>
  #ifdef WITH_UMFPACK
  using UMFPACKType = Eigen::UmfPackLU<MatrixType>;
  #else
  using UMFPACKType = bool;
  #endif
  
  #ifdef WITH_PASTIX
  using PastixLUType = Eigen::PastixLU<Eigen::SparseMatrix<double>>;
  using PastixLLTType = Eigen::PastixLLT<Eigen::SparseMatrix<double>, Eigen::Lower>;
  #else
  using PastixLUType = bool;
  using PastixLLTType = bool;  
  #endif
  

  /// Map to associate to each lowercase name a solver
  /** The user passes a string and this map identifies the solver (this is case-independent)*/
  std::map<std::string, SolverName> map_solver = {{"eigenlu",EigenLU},
                                            {"eigenbicgstab",EigenBiCGStab},
                                            {"pardisolu",PardisoLU},
                                            {"umfpack",UMFPACK},
                                            {"pastixlu",PaStiXLU},
                                            {"pastixllt",PaStiXLLT},};

  /// Map to associate to each solver its proper name
  std::map<SolverName, std::string> map_realname = {{EigenLU,"Eigen LU"},
                                            {EigenBiCGStab,"Eigen BiCGStab"},
                                            {PardisoLU,"Pardiso LU"},
                                            {UMFPACK,"UMFPACK"},
                                            {PaStiXLU,"PaStiX LU"},
                                            {PaStiXLLT,"PaStiX LLT"}};
                                            
  std::map<SolverName, size_t> map_id = {{EigenLU,0},
                                            {EigenBiCGStab,1},
                                            {PardisoLU,2},
                                            {UMFPACK,3},
                                            {PaStiXLU,4},
                                            {PaStiXLLT,5}};

  
  // A tuple of unique pointers to list all available pointers (or "int" as a placeholder in case a pointer is not available),
  // in the order in which they appear in SolverName. If "name" is a SolverName, a get<name> accesses the pointer corresponding to that solver.
  template<typename MatrixType>
  using ListSolvers = std::tuple<std::unique_ptr<EigenLUType<MatrixType>>,
                                 std::unique_ptr<EigenBiCGStabType<MatrixType>>,
                                 std::unique_ptr<PardisoLUType<MatrixType>>,
                                 std::unique_ptr<UMFPACKType<MatrixType>>,
                                 std::unique_ptr<PastixLUType>,
                                 std::unique_ptr<PastixLLTType>>;
  
  
  /// This structure is a wrapper to allow a given code to select various linear solvers
  template<typename MatrixType>
  class LinearSolver
  {
    public:
    /// Constructor
    /** The name of the solver is not caps-dependent
    */
    LinearSolver(const std::string & namesolver)
         : m_name(map_solver[boost::algorithm::to_lower_copy(namesolver)])
           {
            // We adjust the solver if the selected one was not available
            #ifndef WITH_MKL
            if (m_name == PardisoLU){
              m_name = EigenLU;
              std::cout << "[linearsolver] Pardiso not available, reverting to Eigen LU" << std::endl;
            };
            #endif
            #ifndef WITH_UMFPACK
            if (m_name == UMFPACK){
              m_name = EigenLU;
              std::cout << "[linearsolver] UMFPACK not available, reverting to Eigen LU" << std::endl;
            };
            #endif
            #ifndef WITH_PASTIX
            if (m_name == PaStiXLU || m_name == PaStiXLLT){
              m_name = EigenLU;
              std::cout << "[linearsolver] PaStiX not available, reverting to Eigen LU" << std::endl;               
            };
            #endif
            
            // Instantiate selected solver, put all other pointers to NULL
            switch(m_name){
              case EigenLU:
              {
                std::get<0>(m_list_solvers).reset( new EigenLUType<MatrixType> );
                break;
              }
              case EigenBiCGStab:
              {
                std::get<1>(m_list_solvers).reset( new EigenBiCGStabType<MatrixType> );                  
                break;
              }
              #ifdef WITH_MKL
              case PardisoLU:
              {
                unsigned nb_threads_hint = std::thread::hardware_concurrency();
                mkl_set_dynamic(0);
                mkl_set_num_threads(nb_threads_hint);
                std::get<2>(m_list_solvers).reset( new PardisoLUType<MatrixType> );
                break;
              }
              #endif
              #ifdef WITH_UMFPACK
              case UMFPACK:
              {
                std::get<3>(m_list_solvers).reset( new UMFPACKType<MatrixType>);
                break;
              }
              #endif
              #ifdef WITH_PASTIX
              case PaStiXLU:
              {
                std::get<4>(m_list_solvers).reset( new PastixLUType );
                std::get<4>(m_list_solvers).get()->init(SolverParam[0],SolverParam[1]);
                break;
              }
              case PaStiXLLT:
              {
                std::get<5>(m_list_solvers).reset( new PastixLLTType );
                break;
              }
              #endif
              default:
                break;
            } // end switch
         }; // end constructor
    
    /// Returns the name of the solver
    std::string name() const
    {
      return map_realname[m_name];
    }
    
    /// Returns the information message after the "factorize" step
    Eigen::ComputationInfo info_factorize() const
    {
      return m_info_factorize;
    }
    
    /// Returns the information message after the "solve" step
    Eigen::ComputationInfo info_solve() const
    {
      return m_info_solve;
    }

    /// Analyze the pattern of the matrix
    void analyzePattern(MatrixType & A  ///< Matrix of the system to solve
                       )
    {
      switch(m_name){
        case EigenLU:
          {
            std::get<0>(m_list_solvers).get()->analyzePattern(A);
            break;
          }
        case EigenBiCGStab:
          {
            std::get<1>(m_list_solvers).get()->analyzePattern(A);
            break;
          }
        #ifdef WITH_MKL
        case PardisoLU:
          {
            std::get<2>(m_list_solvers).get()->analyzePattern(A);
            break;
          }
        #endif
        
        #ifdef WITH_UMFPACK
        case UMFPACK:
          {
            std::get<3>(m_list_solvers).get()->analyzePattern(A);
            break;
          }
        #endif

        #ifdef WITH_PASTIX
        case PaStiXLU:
          {
            std::get<4>(m_list_solvers).get()->analyzePattern(A);
            break;
          }

        case PaStiXLLT:
          {
            std::get<5>(m_list_solvers).get()->analyzePattern(A);
            break;
          }
        #endif
                
        default:
            break;
      } // end switch for m_name 
      
    }
    
    /// Factorize the matrix
    void factorize(MatrixType & A  ///< Matrix of the system to solve
                   )
    {
      switch(m_name){
        case EigenLU:
          {
            std::get<0>(m_list_solvers).get()->factorize(A);
            m_info_factorize = std::get<0>(m_list_solvers).get()->info();
            break;
          }
        case EigenBiCGStab:
          {
            std::get<1>(m_list_solvers).get()->factorize(A);
            m_info_factorize = std::get<1>(m_list_solvers).get()->info();
            break;
          }
        #ifdef WITH_MKL
        case PardisoLU:
          {
            std::get<2>(m_list_solvers).get()->factorize(A);
            m_info_factorize = std::get<2>(m_list_solvers).get()->info();
            break;
          }
        #endif
        
        #ifdef WITH_UMFPACK
        case UMFPACK:
          {
            std::get<3>(m_list_solvers).get()->factorize(A);
            m_info_factorize = std::get<3>(m_list_solvers).get()->info();
            break;
          }
        #endif

        #ifdef WITH_PASTIX
        case PaStiXLU:
          {
            std::get<4>(m_list_solvers).get()->factorize(A);
            m_info_factorize = std::get<4>(m_list_solvers).get()->info();
            break;
          }

        case PaStiXLLT:
          {
            std::get<5>(m_list_solvers).get()->factorize(A);
            m_info_factorize = std::get<5>(m_list_solvers).get()->info();
            break;
          }
        #endif
                
        default:
            break;
      } // end switch for m_name 

      _check_info(m_info_factorize, "factorize");

    }

    /// Analyze and factorize the matrix
    void compute(MatrixType & A  ///< Matrix of the system to solve
                 )
    {
      analyzePattern(A);
      factorize(A);
    }

    /// Solve the system Ax=b using the selected solver (after analysis of pattern and computation of matrix)
    template<typename VectorType>
    VectorType solve(VectorType & b   ///< Right-hand side of the system to solve
                     )
    {
      VectorType x;
      
      switch(m_name){
        case EigenLU:
          {
            x = std::get<0>(m_list_solvers).get()->solve(b);
            m_info_solve = std::get<0>(m_list_solvers).get()->info();
            break;
          }
        case EigenBiCGStab:
          {
            x = std::get<1>(m_list_solvers).get()->solve(b);
            m_info_solve = std::get<1>(m_list_solvers).get()->info();
            break;
          }
        #ifdef WITH_MKL
        case PardisoLU:
          {
            x = std::get<2>(m_list_solvers).get()->solve(b);
            m_info_solve = std::get<2>(m_list_solvers).get()->info();
            break;
          }
        #endif
        
        #ifdef WITH_UMFPACK
        case UMFPACK:
          {
            x = std::get<3>(m_list_solvers).get()->solve(b);
            m_info_solve = std::get<3>(m_list_solvers).get()->info();
            break;
          }
        #endif

        #ifdef WITH_PASTIX
        case PaStiXLU:
          {
            x = std::get<4>(m_list_solvers).get()->solve(b);
            m_info_solve = std::get<4>(m_list_solvers).get()->info();
            break;
          }

        case PaStiXLLT:
          {
            x = std::get<5>(m_list_solvers).get()->solve(b);
            m_info_solve = std::get<5>(m_list_solvers).get()->info();
            break;
          }
        #endif
                
        default:
            break;
      } // end switch for m_name 

      _check_info(m_info_solve, "solve");

      return x;
    };

    /// Perform all operators to solve Ax=b
    template<typename VectorType>
    VectorType compute_and_solve(MatrixType & A,  ///< Matrix of the system to solve
                                 VectorType & b   ///< Right-hand side of the system
                                 )
    {
      analyzePattern(A);
      factorize(A);
      return solve(b);
    }

    /// Check relative infinity norm of residual: ||Ax-b||/||b||
    template <typename VectorType>
    inline double residual(MatrixType & A, VectorType & b, VectorType & x) 
      { return (A*x - b).template lpNorm<Eigen::Infinity>() / (b.template lpNorm<Eigen::Infinity>()+1e-8); };

    private:
  
    // Checks an info, and stop the execution if there is an error
    void _check_info(Eigen::ComputationInfo const & info, std::string step) const
    {
      if (info != Eigen::Success) {
        std::cerr << "[linearsolver] ERROR in step: " << step << std::endl;
        std::cerr << "[linearsolver] Info: " << info << " (see Eigen::ComputationInfo for the meaning)" << std::endl;
        exit(1);
      }
    }

    // Members
    SolverName m_name;
    ListSolvers<MatrixType> m_list_solvers;
    Eigen::ComputationInfo m_info_factorize;
    Eigen::ComputationInfo m_info_solve;
  }; // end definition class


  //@}
  
} // end namespace HArDCore2D


#endif
