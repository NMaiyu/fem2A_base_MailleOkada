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
        
        double sinus2_fct( vertex v )
        {
            return 2 * pow(2, M_PI)*sin(M_PI * v.x) * sin(M_PI * v.y);
        }
        
        double sin_fct(vertex v)
        {
            return sin(M_PI * v.y);
        }
        
        double edge_top( vertex v)
        {
            return floor(v.y);
        }
        
        double edge_bottom( vertex v)
        {
            return -1*(floor(v.y-1));
        }

        double edge_right( vertex v)
        {
            return floor(v.x);
        }
        
        double edge_left( vertex v)
        {
            return (floor(v.x)-1)*-1;
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
            Quadrature quadrat = Quadrature::get_quadrature(2,false);
            
            // CREATE EMPTY K AND F
            int t_max;
            t_max = M.nb_vertices();
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
            
            // CREATE SOLUTION FILE
            std::string solution_filename;
            solution_filename.assign(mesh_filename.begin(),mesh_filename.end()-4);
            solution_filename.append("bb");
            save_solution(x, solution_filename);
            
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
            t_max = M.nb_vertices() ;
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
            
            // CREATE SOLUTION FILE
            std::string solution_filename;
            solution_filename.assign(mesh_filename.begin(),mesh_filename.end()-4);
            solution_filename.append("bb");
            save_solution(x, solution_filename);
            
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
            t_max = M.nb_vertices() ;
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
                assemble_elementary_vector(element, fonctions, quadrat, sinus2_fct, Fe);
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
            
            // CREATE SOLUTION FILE
            std::string solution_filename;
            solution_filename.assign(mesh_filename.begin(),mesh_filename.end()-4);
            solution_filename.append("bb");
            save_solution(x, solution_filename);
            
        }
        
        
        
        void neumann_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a neumann problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            
            
            // LOAD MESH AND SET ATTRIBUTES
            Mesh M;
            M.load(mesh_filename);
            M.set_attribute( unit_fct, 1, true ); // Null Neumann
            M.set_attribute( edge_right, 0, true ); // Dirichlet
            
            // CHOOSE QUADRATURE AND SHAPE FUNCTIONS
            ShapeFunctions fonctions(2,1);
            ShapeFunctions fonctions_1D(1,1);
            Quadrature quadrat = Quadrature::get_quadrature(0,false);
            Quadrature quadrat_1D = Quadrature::get_quadrature(0,true);
            
            // CREATE EMPTY K, F
            int t_max;
            t_max = M.nb_vertices() ;
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
            
            
            for (int b=0 ; b<M.nb_edges(); ++b)
            {

                ElementMapping element_1D(M,true,b);
                assemble_elementary_neumann_vector(element_1D, fonctions_1D, quadrat_1D, sin_fct,Fe);
                local_to_global_vector(M, true,b, Fe, F);
            }
                       
            
            
            // CHOOSE AND APPLY DIRICHLET CONDITIONS
            std::vector< bool > attribute_bool(2, false);
            attribute_bool[0]=true;
            std::vector< double > values(M.nb_vertices(),0);
            apply_dirichlet_boundary_conditions(M, attribute_bool, values, K, F);
            
            std::cout << "F\n";
            for (double i :F)
            {
                std::cout << i << ' ';
            }
            std::cout <<'\n';
            
            
            
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
            
            // CREATE SOLUTION FILE
            std::string solution_filename;
            solution_filename.assign(mesh_filename.begin(),mesh_filename.end()-4);
            solution_filename.append("bb");
            save_solution(x, solution_filename);
            
        }

    }

}
