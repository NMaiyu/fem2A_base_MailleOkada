#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

double unit_fct( FEM2A::vertex v )
{
    return 1.;
}

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quadrature()
        {
            // Creation de la quadrature d'ordre 2 dans un triangle
            std::cout <<"QUADRATURE TEST\n";
            Quadrature quadrat_test = Quadrature::get_quadrature(2);
            
            // Calcul de l'integrale sur le triangle
            double poids=0;
            for(int i=0; i<quadrat_test.nb_points() ;++i){
                poids += quadrat_test.weight(i);
            }
            std::cout <<"Triangle area "<<poids <<std::endl;
            return true;
        }
        
        bool test_elementMapping()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            int element_index = 4;
            ElementMapping element(mesh, false, element_index,true);
            
            vertex point;
            point.x =0.2; 
            point.y = 0.4;
            element.transform(point, true);
            
            element.jacobian(point, true);
             
            return true; 
        }
        
        bool test_ShapeFunction()
        {
            ShapeFunctions fonctions(2,1);
            
            vertex point;
            point.x =0.2; 
            point.y = 0.4;
            
            std::cout<<fonctions.nb_functions()<<std::endl;
            std::cout<<fonctions.evaluate(2,point)<<std::endl;
            fonctions.evaluate_grad(2,point, true);
            
            return true;
        }
        
        bool test_ElementaryMatrix()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            int element_index = 4;
            ElementMapping element(mesh, false, element_index);
            ShapeFunctions fonctions(2,1);
            Quadrature quadrat = Quadrature::get_quadrature(6,false);
            DenseMatrix Ke;
            

            
            assemble_elementary_matrix(element, fonctions, quadrat, unit_fct, Ke);
            return true;
        }
        
        bool test_LocalToGlobal()
        {           
            Mesh mesh;
            mesh.load("data/square.mesh");
            ShapeFunctions fonctions(2,1);
            Quadrature quadrat = Quadrature::get_quadrature(0,false);
            int t_max;
            t_max = mesh.nb_vertices() * quadrat.nb_points();
            std::cout << t_max<<std::endl;
            DenseMatrix Ke ;
            SparseMatrix K(t_max);
                        
            
            for (int t=0 ; t<mesh.nb_triangles(); ++t)
            {
                ElementMapping element(mesh, false, t);
                assemble_elementary_matrix(element, fonctions, quadrat, unit_fct, Ke);
                local_to_global_matrix(mesh, t, Ke, K);
            }
            
            K.print();
            
            
            return true;
        }
        
        bool test_BdrConditions()
        {
            Mesh M;
            M.load("data/square.mesh");
            
            M.set_attribute( unit_fct, 0, true );
            
            ShapeFunctions fonctions(2,1);
            Quadrature quadrat = Quadrature::get_quadrature(0,false);
            int t_max;
            t_max = M.nb_vertices() * quadrat.nb_points();
            DenseMatrix Ke ;
            SparseMatrix K(t_max);
            
            
            for (int t=0 ; t<M.nb_triangles(); ++t)
            {
                ElementMapping element(M, false, t);
                assemble_elementary_matrix(element, fonctions, quadrat, unit_fct, Ke);
                local_to_global_matrix(M, t, Ke, K);
            }
            
            std::vector< double > F(M.nb_vertices(), 1);
            
            std::vector< bool > attribute_bool(1, true);
            std::vector< double > values(M.nb_vertices());
            
            for (int i =0 ; i<M.nb_vertices(); ++i)
            {
                values[i] = M.get_vertex(i).x + M.get_vertex(i).y; 
            }
            
            apply_dirichlet_boundary_conditions(M, attribute_bool, values, K, F);
            
            std::vector< double > x(M.nb_vertices(), 0);
            solve(K,F, x);
            
            for (double i :x)
            {
                std::cout << i << ' ';
            }
            std::cout <<'\n';
            
            return true;
        }
        
        bool test_ElementaryVector()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            int element_index = 4;
            ElementMapping element(mesh, false, element_index);
            ShapeFunctions fonctions(2,1);
            Quadrature quadrat = Quadrature::get_quadrature(6,false);
            std::vector< double > Fe;
            assemble_elementary_vector(element, fonctions, quadrat, unit_fct, Fe);
            
            std::cout << "Fe (pour f=1)\n";
            for (double i :Fe)
            {
                std::cout << i << " ";
            }
            std::cout <<"\n";
            return true;
        }
        
        bool test_LocToGlobVector()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            ShapeFunctions fonctions(2,1);
            Quadrature quadrat = Quadrature::get_quadrature(0,false);
            int t_max;
            t_max = mesh.nb_vertices() * quadrat.nb_points();
            std::cout << t_max<<std::endl;
            std::vector< double > Fe(3,0);
            std::vector< double > F(t_max,0);
                        
            
            for (int t=0 ; t<mesh.nb_triangles(); ++t)
            {
                ElementMapping element(mesh, false, t);
                assemble_elementary_vector(element, fonctions, quadrat, unit_fct, Fe);
                local_to_global_vector(mesh, false,t, Fe, F);
            }
            
            std::cout << "F (pour f=1)\n";
            for (double i :F)
            {
                std::cout << i << " ";
            }
            std::cout <<"\n";
            return true;
        }
        
        
        bool test_Neumann()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");
            int element_index = 4;
            ElementMapping element_1D(mesh, true, element_index);
            ShapeFunctions fonctions_1D(1,1);
            Quadrature quadrat_1D = Quadrature::get_quadrature(0,true);
            std::vector< double > Fe;
            assemble_elementary_neumann_vector(element_1D, fonctions_1D, quadrat_1D, unit_fct,Fe);
            
            std::cout << "Fe (pour f=1)\n";
            for (double i :Fe)
            {
                std::cout << i << " ";
            }
            std::cout <<"\n";
            return true;        
        }
    }
}
