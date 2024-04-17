#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }
        
        double sinus_fct( vertex v )
        {
            return 2 * pow(2, M_PI)*sin(M_PI * v.x) * sin(M_PI * v.y);
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            
            // LOAD MESH AND SET ATTRIBUTES
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct, 0, true );
            
            // CHOOSE QUADRATURE AND SHAPE FUNCTIONS
            ShapeFunctions fonctions(2,1);
            Quadrature quadrat = Quadrature::get_quadrature(0,false);
            
            // CREATE EMPTY K AND F
            int t_max;
            t_max = M.nb_vertices() * quadrat.nb_points();
            DenseMatrix Ke ;
            SparseMatrix K(t_max);
            std::vector< double > F(t_max,0);
            
            // CREATE K
            for (int t=0 ; t<M.nb_triangles(); ++t)
            {
                ElementMapping element(M, false, t);
                // K
                assemble_elementary_matrix(element, fonctions, quadrat, unit_fct, Ke);
                local_to_global_matrix(M, t, Ke, K);
            }
            
            // CHOOSE AND APPLY DIRICHLET CONDITIONS            
            std::vector< bool > attribute_bool(1, true);
            std::vector< double > values(M.nb_vertices());
            for (int i =0 ; i<M.nb_vertices(); ++i)
            {
                values[i] = M.get_vertex(i).x + M.get_vertex(i).y; 
            }
            apply_dirichlet_boundary_conditions(M, attribute_bool, values, K, F);
            
            // SOLVE
            std::vector< double > x(t_max, 0);
            solve(K,F, x);
            
            
            // PRINT SOLUTION
            std::cout<<"Vector x \n";
            for (double i :x)
            {
                std::cout << i << ' ';
            }
            
            std::cout<<"Vector F \n";
            for (double i :F)
            {
                std::cout << i << ' ';
            }            
            std::cout <<'\n';
            
            save_solution(x, "data/square_fine.bb");
            
        }
        
        
        
        void source_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a Dirichlet problem with source" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            
            // LOAD MESH AND SET ATTRIBUTES
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct, 0, true );
            
            // CHOOSE QUADRATURE AND SHAPE FUNCTIONS
            ShapeFunctions fonctions(2,1);
            Quadrature quadrat = Quadrature::get_quadrature(0,false);
            
            // CREATE EMPTY K, F
            int t_max;
            t_max = M.nb_vertices() * quadrat.nb_points();
            DenseMatrix Ke ;
            SparseMatrix K(t_max);
            std::vector< double > Fe;
            std::vector< double > F(t_max,0);
            
            
            // CREATE MATRIX K AND VECTOR F
            for (int t=0 ; t<M.nb_triangles(); ++t)
            {
                ElementMapping element(M, false, t);
                // K
                assemble_elementary_matrix(element, fonctions, quadrat, unit_fct, Ke);
                local_to_global_matrix(M, t, Ke, K);
                
                // F
                assemble_elementary_vector(element, fonctions, quadrat, unit_fct, Fe);
                local_to_global_vector(M, false,t, Fe, F);
            }            
            
            // CHOOSE AND APPLY DIRICHLET CONDITIONS
            std::vector< bool > attribute_bool(1, true);
            std::vector< double > values(M.nb_vertices(),0);
            apply_dirichlet_boundary_conditions(M, attribute_bool, values, K, F);
            
            // SOLVE
            std::vector< double > x(M.nb_vertices(), 0);
            solve(K,F, x);
            
            // PRINT SOLUTION
            std::cout<<"Vector x \n";
            for (double i :x)
            {
                std::cout << i << ' ';
            }
            std::cout <<'\n';
            
            save_solution(x, "data/square.bb");
            
        }
        
        
        
        
        void sinus_bump_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a sinus bump problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            
            // LOAD MESH AND SET ATTRIBUTES
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct, 0, true );
            
            // CHOOSE QUADRATURE AND SHAPE FUNCTIONS
            ShapeFunctions fonctions(2,1);
            Quadrature quadrat = Quadrature::get_quadrature(0,false);
            
            // CREATE EMPTY K, F
            int t_max;
            t_max = M.nb_vertices() * quadrat.nb_points();
            DenseMatrix Ke ;
            SparseMatrix K(t_max);
            std::vector< double > Fe;
            std::vector< double > F(t_max,0);
            
            
            // CREATE MATRIX K AND VECTOR F
            for (int t=0 ; t<M.nb_triangles(); ++t)
            {
                ElementMapping element(M, false, t);
                // K
                assemble_elementary_matrix(element, fonctions, quadrat, unit_fct, Ke);
                local_to_global_matrix(M, t, Ke, K);
                
                // F
                assemble_elementary_vector(element, fonctions, quadrat, sinus_fct, Fe);
                local_to_global_vector(M, false,t, Fe, F);
            }            
            
            // CHOOSE AND APPLY DIRICHLET CONDITIONS
            std::vector< bool > attribute_bool(1, true);
            std::vector< double > values(M.nb_vertices(),0);
            apply_dirichlet_boundary_conditions(M, attribute_bool, values, K, F);
            
            // SOLVE
            std::vector< double > x(M.nb_vertices(), 0);
            solve(K,F, x);
            
            // PRINT SOLUTION
            std::cout<<"Vector x \n";
            for (double i :x)
            {
                std::cout << i << ' ';
            }
            std::cout <<'\n';
            
            save_solution(x, "data/square.bb");
            
        }

    }

}
