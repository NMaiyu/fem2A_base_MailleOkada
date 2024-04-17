#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quad = false;
    const bool t_elmapping = false;
    const bool t_shapefunct = false;
    const bool t_elemM = false;
    const bool t_loc_to_glob = false;
    const bool t_bdr_cond = false;
    const bool t_elemV = false;
    const bool t_loc_glob_vector = true;

    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if( t_quad ) Tests::test_quadrature();
    if (t_elmapping ) Tests::test_elementMapping();
    if (t_shapefunct) Tests::test_ShapeFunction();
    if (t_elemM) Tests::test_ElementaryMatrix();
    if (t_loc_to_glob) Tests::test_LocalToGlobal();
    if (t_bdr_cond) Tests::test_BdrConditions();
    if (t_elemV) Tests::test_ElementaryVector();
    if (t_loc_glob_vector) Tests::test_LocToGlobVector();
}

void run_simu()
{

    const bool simu_pure_dirichlet = false;
    const bool simu_source_dirichlet = true ; 

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) {
        Simu::pure_dirichlet_pb("data/square.mesh", verbose);
    }
    if( simu_source_dirichlet){
        Simu::source_dirichlet_pb("data/square.mesh", verbose);
    }
}

int main( int argc, const char * argv[] )
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
