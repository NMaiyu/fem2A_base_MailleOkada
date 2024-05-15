#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        //std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border)
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            //std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i)
        : border_( border )
    {
        if ( border ) 
        {
            // Gets the vertices of the element
            vertices_.push_back(M.get_edge_vertex(i, 0));
            vertices_.push_back(M.get_edge_vertex(i,1));
        }
        
        else 
        {
            vertices_.push_back(M.get_triangle_vertex(i, 0));
            vertices_.push_back(M.get_triangle_vertex(i,1));
            vertices_.push_back(M.get_triangle_vertex(i,2));
        }
    }
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i, bool verbose)
        : border_( border )
    {
        if(verbose){std::cout << "[ElementMapping] constructor for element " << i << " ";}
        // TODO
        if ( border ) 
        {
            //std::cout << "(edge border)\n";
            
            // Gets the vertices of the element
            vertices_.push_back(M.get_edge_vertex(i, 0));
            vertices_.push_back(M.get_edge_vertex(i,1));
            
            // Prints the vertices
            if (verbose){ 
                std::cout<<"\nBorder vertices \nx y\n";
                std::cout<<vertices_[0].x<< " "<<vertices_[0].y<< std::endl;
                std::cout<<vertices_[1].x<< " "<<vertices_[1].y<< std::endl;
                std::cout <<"\n";
                }
        }
        
        else 
        {
            //std::cout << "(triangle border)\n";
            
            // Gets the vertices of the element
            if (verbose){ 
                std::cout<<"\nTriangle vertices \n";
                vertices_.push_back(M.get_triangle_vertex(i, 0));
                vertices_.push_back(M.get_triangle_vertex(i,1));
                vertices_.push_back(M.get_triangle_vertex(i,2));
                
            
                // Prints the vertices
                std::cout<<"x y\n";
                std::cout<<vertices_[0].x<< " "<<vertices_[0].y<< std::endl;
                std::cout<<vertices_[1].x<< " "<<vertices_[1].y<< std::endl;
                std::cout<<vertices_[2].x<< " "<<vertices_[2].y<< std::endl;
                std::cout <<"\n";
                }
        }
    }

    vertex ElementMapping::transform( vertex x_r , bool verbose) const
    {
        if(verbose){
        std::cout << "\n[ElementMapping] transform reference to world space" << '\n';
        }
        // TODO
        vertex r ;
        
        if (border_)
        {
            r.x = (1-x_r.x) * vertices_[0].x + x_r.x*vertices_[1].x;
            r.y = (1-x_r.x) * vertices_[0].y + x_r.x*vertices_[1].y;
        }
        
        else
        {
            r.x = (1 - x_r.x - x_r.y)*vertices_[0].x + x_r.x*vertices_[1].x + x_r.y*vertices_[2].x;
            r.y = (1 - x_r.x - x_r.y)*vertices_[0].y + x_r.x*vertices_[1].y + x_r.y*vertices_[2].y; 
        }
        
        
        // Prints x_r and r
        if(verbose){
            std::cout << "Reference vertex "<< x_r.x << " "<<x_r.y<<std::endl;
            std::cout << "World space vertex "<< r.x << " " << r.y<<"\n"<<std::endl;
            }
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r , bool verbose) const
    {
        if(verbose){
            std::cout << "\n[ElementMapping] compute jacobian matrix" << '\n';
            std::cout << "Reference vertex (xi="<<x_r.x<<", eta="<<x_r.y<<")\n";
        }
        // TODO
        DenseMatrix J;
        
        if(border_){
            J.set_size(2,1);
            J.set(0,0,vertices_[1].x - vertices_[0].x);
            J.set(1,0, vertices_[1].y - vertices_[0].y);
        }
        else
        {
            if(verbose){std::cout<<"triangle"<<std::endl;}
            J.set_size(2,2);
            J.set(0,0,vertices_[1].x - vertices_[0].x);
            J.set(0,1,vertices_[2].x - vertices_[0].x);
            J.set(1,0,vertices_[1].y - vertices_[0].y);
            J.set(1,1,vertices_[2].y - vertices_[0].y);
        }
        
       if(verbose){
           std::cout << "Matrice jacobienne\n";
           std::cout <<J.get(0,0)<<" "<<J.get(0,1)<<" \n"<<J.get(1,0)<<" "<<J.get(1,1)<<"\n"<<std::endl;
       }
       return J ;
    }

    double ElementMapping::jacobian( vertex x_r , bool verbose) const
    {
        if(verbose){
            std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        }
        // TODO
        DenseMatrix J;
        J = jacobian_matrix(x_r); 
        double det;
        if(border_){
            det = std::sqrt(J.get(0,0)*J.get(0,0) + J.get(1,0)*J.get(1,0));
            }
        else {
            det = J.det_2x2();
            }
        
       if(verbose){std::cout << "Determinant\n" << det<<"\n"<<std::endl;}
        return det ;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        //std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
        // TODO
        assert(order==1);
        assert(dim ==1 or dim ==2);
        
    }

    int ShapeFunctions::nb_functions() const
    {
        //std::cout << "[ShapeFunctions] number of functions" << '\n';
        // TODO
        
        return dim_+1 ;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r, bool verbose) const
    {
        if (verbose){
            std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
            std::cout << "Reference vertex (xi="<<x_r.x<<", eta="<<x_r.y<<")\n";
        }
        // TODO
        double nb;
        
        // Pour un edge
        if(dim_==1){
            if(i ==0){nb = 1-x_r.x;}
            else{nb = x_r.x;}
        }
        
        
        // Pour un triangle
        else{
            switch (i){
            case 0:
                nb = 1 - x_r.x - x_r.y;
                break;
            case 1:
                nb = x_r.x;
                break;
            case 2:
                nb = x_r.y;
            }
        }
        
        return nb ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r, bool verbose) const
    {
        if(verbose){
            std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
            std::cout << "Reference vertex (xi="<<x_r.x<<", eta="<<x_r.y<<")\n";
            }
        // TODO
        vec2 g ;
        
        if(dim_==1){
            if(i ==0){
                g.x = -1;
                g.y = 0;
                }
            else{
                g.x = 1;
                g.y = 0;
                }
        }
        
        
        // Pour un triangle
        else{
            switch (i){
            case 0:
                g.x = -1;
                g.y = -1;
                break;
            case 1:
                g.x = 1;
                g.y = 0;
                break;
            case 2:
                g.x = 0;
                g.y = 1;
                break;
            }
        }
        
        // Prints the gradient
        //std::cout<<g.x<< " "<<g.y<< std::endl;
        //std::cout <<"\n";
        
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke,
        bool verbose )
    {
        if(verbose){std::cout << "\ncompute elementary matrix" << '\n';}
        
        int ijmax;
        DenseMatrix je_inv;
        vec2 nabla_shape_i;
        vec2 nabla_shape_j;
        vec2 nabla_base_i;
        vec2 nabla_base_j;
        
        
        ijmax = reference_functions.nb_functions();
        Ke.set_size(ijmax, ijmax);
        
        for (int i =0; i<ijmax; ++i)
            {
            for (int j=0 ; j<ijmax; ++j)
                {
                double keij;
                keij =0;
                for(int q = 0; q< quadrature.nb_points(); ++q)
                    {
                    je_inv = elt_mapping.jacobian_matrix(quadrature.point(q)).invert_2x2().transpose();
                    nabla_shape_i = reference_functions.evaluate_grad(i,quadrature.point(q));
                    nabla_shape_j = reference_functions.evaluate_grad(j,quadrature.point(q));
                        
                    nabla_base_i = je_inv.mult_2x2_2(nabla_shape_i);
                    nabla_base_j = je_inv.mult_2x2_2(nabla_shape_j);
                        
                    keij+= quadrature.weight(q)*coefficient(elt_mapping.transform(quadrature.point(q)))* dot(nabla_base_i, nabla_base_j) * elt_mapping.jacobian(quadrature.point(q));
                        
                    }
                Ke.set(i,j,keij);

                }
            }
            
    }
       

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K,
        bool verbose )
    {
        if(verbose){std::cout << "Ke -> K for triangle no." << t << " \n";}
        int global_i;
        int global_j;
        
        // TODO
        for (int i=0 ; i<Ke.height() ; ++i)
        {
            global_i = M.get_triangle_vertex_index(t, i);
            for (int j = 0 ; j<Ke.height() ; j++)
            {
                global_j = M.get_triangle_vertex_index(t, j);
                K.add(global_i,global_j , Ke.get(i,j));
                if(verbose)
                {
                    std::cout<<"Le point de Ke en ("<< i <<","<<j<<") correspond à ("<<global_i<<","<<global_j<<") dans K" <<std::endl;
                }

            }
        }
        
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe,
        bool verbose )
    {
        if (verbose)
        {std::cout << "\ncompute elementary vector (source term)" << '\n';}
        // TODO
        
        int imax;
        double shape_i;
        
        
        imax = reference_functions.nb_functions();
        Fe.resize(imax,0);
        if(verbose){std::cout << "create Fe" << '\n';}
        
        for (int i =0; i<imax; ++i)
            {
            Fe[i] =0;
            for(int q = 0; q< quadrature.nb_points(); ++q)
                {
                shape_i = reference_functions.evaluate(i,quadrature.point(q));
                Fe[i]+= quadrature.weight(q)*source(elt_mapping.transform(quadrature.point(q)))* shape_i * elt_mapping.jacobian(quadrature.point(q));
                }
            }
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe,
        bool verbose )
    {
        if (verbose)
        {
        std::cout << "\ncompute elementary neumann vector" << '\n';
        }
        
        int imax;
        double shape_i;
        
        imax = reference_functions_1D.nb_functions();
        Fe.resize(imax,0);
        if(verbose){std::cout << "create Fe" << '\n';}

        for (int i =0; i<imax; ++i)
            {
            Fe[i] =0;
            for(int q = 0; q< quadrature_1D.nb_points(); ++q)
                {
                    vertex point_q = quadrature_1D.point(q);
                    shape_i = reference_functions_1D.evaluate(i,point_q);
                    double w = quadrature_1D.weight(q);
                    Fe[i]+= w*neumann(elt_mapping_1D.transform(point_q))* shape_i * elt_mapping_1D.jacobian(point_q, verbose);
                }
                
            }
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F,
        bool verbose)
    {
        if(verbose){std::cout << "Fe -> F no." << i << " \n";}

        int global_j;
        F.resize(M.nb_vertices());
        
        // TODO
        for (int j=0 ; j<Fe.size() ; ++j)
        {
            if (not border) {global_j = M.get_triangle_vertex_index(i, j);}
            else {global_j = M.get_edge_vertex_index(i, j);}
            if(verbose){std::cout << "Numero global du point " << j<< " : "<< global_j<<std::endl;}
            F[global_j]+=Fe[j];
        }
        
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F, 
        bool verbose)
    {
        if(verbose){std::cout << "\napply dirichlet boundary conditions" << '\n';}
        double P = 10000.;
        std::vector< bool > passage(values.size(), false);
        
        for(int i= 0 ; i<M.nb_edges();++i)
        {
            if(attribute_is_dirichlet[M.get_edge_attribute(i)])
            {
                
                if(not passage[M.get_edge_vertex_index(i,0)])
                {
                    K.add(M.get_edge_vertex_index(i,0), M.get_edge_vertex_index(i,0),P);
                    F[M.get_edge_vertex_index(i,0)] += P*values[M.get_edge_vertex_index(i,0)] ;
                    passage[M.get_edge_vertex_index(i,0)] = true;
                }
                if(not passage[M.get_edge_vertex_index(i,1)])
                {
                    K.add(M.get_edge_vertex_index(i,1), M.get_edge_vertex_index(i,1),P);
                    F[M.get_edge_vertex_index(i,1)] += P*values[M.get_edge_vertex_index(i,1)];
                    passage[M.get_edge_vertex_index(i,1)] = true;
                }    
            }
        }
    }



    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            bool verbose)
    {

        std::cout << "solve poisson problem" << '\n';
            
        std::cout <<"\n[CREATING QUADRATURE AND SHAPE FUNCTIONS]\n";
        ShapeFunctions fonctions(2,1);
        ShapeFunctions fonctions_1D(1,1);
        Quadrature quadrat = Quadrature::get_quadrature(2,false);
        Quadrature quadrat_1D = Quadrature::get_quadrature(2,true);
        
        std::cout <<"\n[CREATE EMPTY K, F]\n";
        int t_max;
        t_max = M.nb_vertices() ;
        SparseMatrix K(t_max);
        std::vector< double > F(t_max,0);
        
        std::cout<<"\n[CREATE MATRIX K AND VECTOR F]\n";
        for (int t=0 ; t<M.nb_triangles(); ++t)
        {
            ElementMapping element(M, false, t);
            // K
            DenseMatrix Ke ;
            assemble_elementary_matrix(element, fonctions, quadrat, diffusion_coef, Ke, verbose);
            local_to_global_matrix(M, t, Ke, K, verbose);
            
            // F
            std::vector< double > Fe;
            assemble_elementary_vector(element, fonctions, quadrat, source_term, Fe, verbose);
            local_to_global_vector(M, false,t, Fe, F, verbose);
        }
        
        std::cout<<"\n[APPLY DIRICHLET CONDITIONS]\n";
        std::vector< bool > attribute_bool(2, false);
        attribute_bool[1]=true;
        std::vector< double > values(M.nb_vertices());
        for(int i=0;i<M.nb_vertices(); i++)
        {
            values[i] = dirichlet_fct(M.get_vertex(i));
        }
        apply_dirichlet_boundary_conditions(M, attribute_bool, values, K, F, verbose);
        
        
        std::cout<<"\n[APPLY NEUMANN CONDITIONS]\n";
        for (int b=0 ; b<M.nb_edges(); ++b)
        {
            if (M.get_edge_attribute(b)==2)
            {
                ElementMapping element_1D(M,true,b);
                std::vector< double > Fe;
                assemble_elementary_neumann_vector(element_1D, fonctions_1D, quadrat_1D, neumann_fct,Fe, verbose);
                local_to_global_vector(M, true,b, Fe, F, verbose);
            }
        }
        
        
        if(verbose)
        {
            std::cout << "\nF\n";
            for (double i :F)
            {
                std::cout << i << ' ';
            }
            std::cout <<'\n';
        }
        

        std::cout<<"\n[SOLVING PROBLEM]\n";
        solve(K,F, solution);
        
        // PRINT SOLUTION
        if(verbose)
        {   
            std::cout<<"\n[SOLUTION]\n";
            for (double sol:solution)
            {
                std::cout << sol << ' ';
            }
            std::cout <<'\n';
        }

    }
    
}
