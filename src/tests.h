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
            ElementMapping element(mesh, false, element_index);
            std::cout << "Creation of element mapping\n";
            
            vertex point;
            point.x =0.2; 
            point.y = 0.4;
            element.transform(point);
            
            
            element.jacobian_matrix(point);
             
            return true; 
        }
    }
}
