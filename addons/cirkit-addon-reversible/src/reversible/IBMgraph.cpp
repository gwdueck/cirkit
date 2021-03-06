/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2017  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "IBMgraph.hpp"
#include <reversible/io/print_circuit.hpp>
#include <reversible/functions/move_qubit.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/gate.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/rotation_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <boost/lexical_cast.hpp>
#include <reversible/functions/ibm_helper.hpp>
#include <reversible/functions/clear_circuit.hpp>
#include <reversible/io/write_qc.hpp>


namespace cirkit
{
    using matrix = std::vector<std::vector<int>>;
    // Here the memory error stopped
    std::vector<TransPath> path_list;
    bool **graph_adjacency = NULL;  // graph structure
    int **trans_cost = NULL;        // the cost of each cnot
    TransPath **trans_path = NULL;  // the transformation path or a given cnot
    
    bool read_graph( const std::string& filename )
    {
        int v,w;
        std::fstream graphfile;
        graphfile.open ( filename, std::fstream::in );
        if ( !graphfile.is_open() )
        {
            return false;
        }
        graphfile >> graph_size;
        graph_adjacency = new bool*[graph_size];
        for(int i = 0; i < graph_size ; i++ )
        {
            graph_adjacency[i] = new bool[graph_size];
            for(int j = 0; j < graph_size; j++)
            {
                graph_adjacency[i][j] = false;
            }
        }
        while ( !graphfile.eof() )
        {
            graphfile >> v >> w;
            graph_adjacency[v][w] = true;
        }
        graphfile.close();
        return true;
    }

    // read the matrix and the transformations from a file
    bool read_from_file ( const std::string& filename )
    { 
        std::vector<std::string> type_name = { "cab", "cba", "tab", "tba", "cabi", "cbai", "tabi", "tbai", "nop", "flip", "cnot3"};
        std::string tmp;
        TransPath tp;
        int a,b,c,v,w;
        std::ifstream graphfile ( filename );
        if ( !graphfile.is_open() )
            return false;
        
        graphfile >> graph_size; // get the number of qubits
        graph_adjacency = new bool*[graph_size];
        allocate_data_stuctures();

        // read the graph
        for(int i = 0; i < graph_size ; i++ ){
            graph_adjacency[i] = new bool[graph_size];
            for(int j = 0; j < graph_size; j++){
                graphfile >> tmp;
                if(tmp == "X")
                    graph_adjacency[i][j] = true;
                else
                    graph_adjacency[i][j] = false;    
            }
        }
        // read the matrix with the transformations costs
        for( v = 0; v < graph_size; v++){ 
            for( w = 0; w < graph_size; w++)
                graphfile >> trans_cost[v][w];
        }
        
        v = 0; // for the variable trans_path
        w = 1; // starting with cnot(0,1)
        while ( !graphfile.eof() ){
            graphfile >> tmp;
            if(tmp == "cost"){ // it means that the transformation has already been read
                trans_path[v][w] = tp;  // update trans_path,
                tp.clear();             // and clear for the next one
                ++w;
                if(v == w) // cnot(x,x) doesn't exist
                    ++w;
                if(w == graph_size){    // finished all w possibilities
                    ++v;                // and start the next one
                    w = 0;
                }
            }
            auto it = std::find(type_name.begin(), type_name.end(), tmp);   
            if(it != type_name.end()){                                  // check if it is a movement
                move_qubit_type pos = (move_qubit_type) std::distance(type_name.begin(), it);    // get the type of movement
                if(tmp == "cnot3"){     // cnot3 needs three parameters
                    graphfile >> a >> b >> c;
                    tp.add( MoveQubit( pos, a, b, c ));
                }
                else{                   // the other movements only two
                    graphfile >> a >> b;
                    tp.add( MoveQubit( pos, a, b ));
                }
            }
        }
        graphfile.close();
        return true;
    }

    // write the graph, the matrix and the transformations to a file
    bool write_to_file ( const std::string& filename )
    {
        std::ofstream graphfile ( filename );
        if ( !graphfile.is_open() )
            return false;

        graphfile << graph_size << std::endl;

        // writing the graph
        for( int i = 0; i < graph_size; i++ ){ 
            for( int j = 0; j < graph_size; j++ ){
                graphfile<< (graph_adjacency[i][j] ? "X " : "- ");
            }
            graphfile << std::endl;
        }
        // writing the matrix
        for( int v = 0; v < graph_size; v++){   
            for( int w = 0; w < graph_size; w++)
                graphfile << trans_cost[v][w] << " ";
            graphfile << std::endl;
        }
        // writing the transformations
        for( int v = 0; v < graph_size; v++){
            for( int w = 0; w < graph_size; w++){
                if( v != w ){
                    graphfile << "cnot(" << v << "," << w << ") => ";
                    trans_path[v][w].print( graphfile );
                }
            }
        }
        graphfile.close();
        return true;
    }

   void print_graph( ){
        for( int i = 0; i < graph_size; i++ ){
            for( int j = 0; j < graph_size; j++ ){
                std::cout << (graph_adjacency[i][j] ? "X " : "- ");
            }
            std::cout << std::endl;
        }
    }

    // print the matrix and the transformations
    void print_matrix( )
    {
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
                std::cout << trans_cost[v][w] << " ";
            std::cout << std::endl;
        }
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
            {
                if( v != w )
                {
                    std::cout << "cnot(" << v << "," << w << ") => ";
                    trans_path[v][w].print( );
                }
            }
        }
    }

    void delete_graph( ){
        for( int i = 0; i < graph_size; i++ ){
            delete [] graph_adjacency[i];
            delete [] trans_path[i];
            delete [] trans_cost[i];
            
        }
        delete [] graph_adjacency;
        delete [] trans_path;
        delete [] trans_cost;
        graph_size = 0;
    }
    
    // When a path is found add it to the global vector of paths.
    void find_all_paths( int v, int w, TransPath &tp, bool *visited )
    {
        bool done = false;
        for( int i = 0; i < graph_size; i++ )
        {
            if( graph_adjacency[v][w] )
            {
                tp.add( MoveQubit( nop, v, w ));
                path_list.push_back( tp );
                tp.remove_last();
                done = true;
            }
            if( graph_adjacency[w][v] )
            {
                tp.add( MoveQubit( flip, v, w ));
                path_list.push_back( tp );
                tp.remove_last();
                done = true;
            }
            if( done ) return; // no need to go further
            if( graph_adjacency[v][i] && !visited[i] )
            {
                tp.add( MoveQubit( cab, v, i ));
                visited[i] = true;
                find_all_paths( i,  w, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[i][v] && !visited[i] )
            {
                tp.add( MoveQubit( cba, v, i ));
                visited[i] = true;
                find_all_paths( i,  w, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[w][i] && !visited[i] )
            {
                tp.add( MoveQubit( tab, w, i ));
                visited[i] = true;
                find_all_paths( v,  i, tp, visited );           
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[i][w] && !visited[i] )
            {
                tp.add( MoveQubit( tba, w, i ));
                visited[i] = true;
                find_all_paths( v,  i, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
        }
    }

    /* Precondition: the path_list contains all paths from v to w
     Postcondition: the best path and its cost will be stored in:
        - trans_cost[v][w]
        - trans_path[v][w]
     */
    
    void set_best_path(int v, int w )
    {
        TransPath tp, best_tp;
        int best_cost;
        best_tp = path_list[0];
        best_cost = best_tp.costPlus();

        for ( auto &p : path_list )
        {
            p.movCnot3();
            if(p.costPlus() < best_cost )
            {
                best_cost = p.costPlus();
                best_tp = p;
            }
        }
        trans_cost[v][w] = best_cost;
        best_tp.addInverse();
        trans_path[v][w] = best_tp;
    }

    void allocate_data_stuctures(){
        trans_cost = new int*[graph_size];
        trans_path = new TransPath*[graph_size];
        for(int i = 0; i < graph_size ; i++ )
        {
            trans_cost[i] = new int[graph_size];
            trans_path[i] = new TransPath[graph_size];
        }
    }
    
    void create_trans( bool verbose )
    {
        allocate_data_stuctures();
        path_list.clear();
        bool visited[graph_size];
        for( int i = 0; i < graph_size ; i++ )
        {
            visited[i] = false;
        }

        TransPath tp;
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
            {
                if( v == w )
                {
                    trans_cost[v][w] = 0;
                }
                else{
                    for( int i = 0; i < graph_size ; i++ ) visited[i] = false;
                    visited[v] = true;
                    visited[w] = true;
                    path_list.clear();
                    tp.clear();
                    find_all_paths( v,  w, tp, visited );
                    set_best_path( v,  w);
                }
            }
        }
        if ( verbose )
        {
            for( int v = 0; v < graph_size; v++)
            {
                for( int w = 0; w < graph_size; w++)
                {
                    trans_cost[v][w] = trans_cost[v][w];
                    // trans_cost[v][w] = trans_cost[v][w]-trans_path[v][w].opt();
                    std::cout << trans_cost[v][w] << " ";
                }
                std::cout << std::endl;
            }
            // std::cout << "== Optimization ==" << std::endl;
            // for( int v = 0; v < graph_size; v++)
            // {
            //     for( int w = 0; w < graph_size; w++)
            //     {
            //         std::cout << trans_cost[v][w]-trans_path[v][w].opt() << " ";
            //     }
            //     std::cout << std::endl;
            // }
            for( int v = 0; v < graph_size; v++)
            {
                for( int w = 0; w < graph_size; w++)
                {
                    if( v != w )
                    {
                        std::cout << "cnot(" << v << "," << w << ") => ";
                        trans_path[v][w].print();
                        // std::cout << "Can reduce: " << trans_path[v][w].opt() << std::endl;
                    }
                }
            }
        }
    }

    void print_movements( const circuit& circ_in )
    {
        unsigned target, control;
        for ( const auto& gate : circ_in )
        {
            if ( is_toffoli( gate ) && gate.controls().size() == 1 )
            {
                target = gate.targets().front();
                control = gate.controls().front().line();
                for ( auto &p : trans_path[control][target].tpath )
                {
                    switch ( p.getType() )
                    {
                        case cab : 
                            std::cout << " cab(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case cba : 
                            std::cout << " cba(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case tab : 
                            std::cout << " tab(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case tba : 
                            std::cout << " tba(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case cabi : 
                            std::cout << " cabi(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case cbai : 
                            std::cout << " cbai(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case tabi : 
                            std::cout << " tabi(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case tbai : 
                            std::cout << " tbai(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case nop : 
                            std::cout << " nop(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case flip : 
                            std::cout << " flip(" << p.getA() << "," << p.getB() << ")"; //<< std::endl;
                            break;
                        case cnot3 :
                            std::cout << " cnot3(" << p.getA() << "," << p.getB() << "," << p.getC() << ")"; //<< std::endl;
                            break;
                        case cnot3i :
                            std::cout << " cnot3i(" << p.getA() << "," << p.getB() << "," << p.getC() << ")"; //<< std::endl;
                            break;                                
                        default : std::cout << "ERROR expand_cnots" << std::endl;
                    }
                }
                std::cout << std::endl;
            }
        }
    }

    circuit transform_tof_clif( const circuit& circ )
    {
        circuit circ_out;
        copy_metadata(circ, circ_out);
        for ( const auto& gate : circ )
        {
            if ( is_toffoli( gate ) )
            {
                if( gate.controls().size() == 1 || gate.controls().empty() )
                {
                    circ_out.append_gate() = gate;
                }
                else if( gate.controls().size() == 2 )
                {
                    unsigned ca, cb, target, aux;
                    bool pa, pb;
                    if( gate.controls().front().line() < gate.controls().back().line())
                    {
                        ca = gate.controls().front().line();
                        cb = gate.controls().back().line();
                        pa = gate.controls().front().polarity();
                        pb = gate.controls().back().polarity();
                    }
                    else
                    {
                        cb = gate.controls().front().line();
                        ca = gate.controls().back().line();
                        pb = gate.controls().front().polarity();
                        pa = gate.controls().back().polarity();
                    }
                    target = gate.targets().front();
                    
                    unsigned tbc, tac, tab, t1, t2, t3;
                    bool ta1, ta2, ta3, ta4;
                    bool tb, tc1, tc2;

                    if(trans_cost[ca][cb] < trans_cost[cb][ca])
                        tab = 2*trans_cost[ca][cb];
                    else
                        tab = 2*trans_cost[cb][ca];
                    
                    if(trans_cost[cb][target] < trans_cost[target][cb])
                        tbc = 2*trans_cost[cb][target];
                    else
                        tbc = 2*trans_cost[target][cb];
                    if(trans_cost[ca][target] < trans_cost[target][ca])
                        tac = 2*trans_cost[ca][target];
                    else
                        tac = 2*trans_cost[target][ca];

                    t1 = 2*trans_cost[ca][target] + 2*trans_cost[cb][target] + tab;
                    t2 = 2*trans_cost[target][ca] + 2*trans_cost[cb][ca] + tbc;
                    t3 = 2*trans_cost[target][cb] + 2*trans_cost[ca][cb] + tac;
                    // std::cout << "t1: " << t1 << " t2: " << t2 << " t3: " << t3 << std::endl;  
                    std::vector<unsigned> controla, controlb, controlt;
                    if( t1 < t2 && t1 < t3 )
                    {
                        if(pa == true && pb == true)
                        {
                            ta2 = ta4 = tc1 = true;
                            ta1 = tb = ta3 = tc2 = false;
                        }
                        else if(pa == true && pb == false)
                        {
                            ta1 = ta4 = tc1 = tc2 = true;
                            tb = ta2 = ta3 = false;
                        }
                        else if(pa == false && pb == true)
                        {
                            tb = ta2 = tc1 = tc2 = true;
                            ta1 = ta3 = ta4 = false;
                        }
                        else if(pa == false && pb == false)
                        {
                            ta2 = ta3 = ta4 = tc2 = true;
                            ta1 = tb = tc1 = false;
                        }
                        controla.push_back(ca);
                        controlb.push_back(cb);
                        controlt.push_back(target);
                        append_hadamard( circ_out, target );
                        append_pauli( circ_out,  ca, pauli_axis::Z, 4u, ta1 );
                        append_pauli( circ_out,  cb, pauli_axis::Z, 4u, tb );
                        append_toffoli( circ_out, controla, target );
                        append_pauli( circ_out,  target, pauli_axis::Z, 4u, ta2 );
                        append_toffoli( circ_out, controlb, target );
                        append_pauli( circ_out,  target, pauli_axis::Z, 4u, ta3 );
                        append_toffoli( circ_out, controla, target );
                        append_pauli( circ_out,  target, pauli_axis::Z, 4u, ta4 );
                        append_toffoli( circ_out, controlb, target );

                        if(trans_cost[ca][cb] < trans_cost[cb][ca])
                        {
                            append_toffoli( circ_out, controla, cb );
                            append_pauli( circ_out,  cb, pauli_axis::Z, 4u, tc1 );
                            append_toffoli( circ_out, controla, cb );
                        }
                        else
                        {
                            append_toffoli( circ_out, controlb, ca );
                            append_pauli( circ_out,  ca, pauli_axis::Z, 4u, tc1 );
                            append_toffoli( circ_out, controlb, ca );
                        }
                        append_pauli( circ_out,  target, pauli_axis::Z, 4u, tc2 );
                        append_hadamard( circ_out, target );
                    }
                    else
                    {
                        if(pa == true && pb == true)
                        {
                            ta1 = ta3 = tb = tc2 = true;
                            ta2 = ta4 = tc1 = false;
                        }
                        else if(pa == false && pb == true)
                        {
                            ta2 = ta4 = tb = tc2 = true;
                            ta1 = ta3 = tc1 = false;
                        }
                        else if(pa == true && pb == false)
                        {
                            ta1 = ta4 = tc1 = tc2 = true;
                            ta2 = ta3 = tb = false;
                        }
                        else if(pa == false && pb == false)
                        {
                            ta2 = ta3 = tc1 = tc2 = true;
                            ta1 = ta4 = tb = false;
                        }
                        if(t3 < t2)
                        {
                            aux = cb;
                            cb = ca;
                            ca = aux;
                            if(pa == true && pb == false)
                            {
                                ta2 = ta4 = tb = tc2 = true;
                                ta1 = ta3 = tc1 = false;
                            }
                            else if(pa == false && pb == true)
                            {
                                ta1 = ta4 = tc1 = tc2 = true;
                                ta2 = ta3 = tb = false;
                            }
                        }
                        controla.push_back(ca);
                        controlb.push_back(cb);
                        controlt.push_back(target);

                        append_hadamard( circ_out, target );
                        append_pauli( circ_out,  ca, pauli_axis::Z, 4u, ta1 );
                        append_pauli( circ_out,  cb, pauli_axis::Z, 4u, tb );
                        append_toffoli( circ_out, controlt, ca );
                        append_pauli( circ_out,  ca, pauli_axis::Z, 4u, ta2 );
                        append_toffoli( circ_out, controlb, ca );
                        append_pauli( circ_out,  ca, pauli_axis::Z, 4u, ta3 );
                        append_toffoli( circ_out, controlt, ca );
                        append_pauli( circ_out,  ca, pauli_axis::Z, 4u, ta4 );
                        append_toffoli( circ_out, controlb, ca );

                        if(trans_cost[cb][target] < trans_cost[target][cb])
                        {
                            append_toffoli( circ_out, controlb, target );
                            append_pauli( circ_out,  target, pauli_axis::Z, 4u, tc1 );
                            append_toffoli( circ_out, controlb, target );
                        }
                        else
                        {
                            append_toffoli( circ_out, controlt, cb );
                            append_pauli( circ_out,  cb, pauli_axis::Z, 4u, tc1 );
                            append_toffoli( circ_out, controlt, cb );
                        }
                        append_pauli( circ_out,  target, pauli_axis::Z, 4u, tc2 );
                        append_hadamard( circ_out, target );
                    }
                }
                else
                {
                    assert( false );
                }
            }
            else
                circ_out.append_gate() = gate;
        }
        return circ_out;
    }

    circuit transform_v_clif( const circuit& circ )
    {
        circuit circ_out;
        copy_metadata(circ, circ_out);
        std::vector<unsigned int> controls;
        unsigned target, control;
        for ( const auto& gate : circ )
        {
            if ( is_v( gate ) )
            {
                const auto& tag = boost::any_cast<v_tag>( gate.type() );
                target = gate.targets().front();
                control = gate.controls().front().line();
                

                if (trans_cost[control][target] < trans_cost[target][control])
                {
                    controls.clear();
                    controls.push_back( control );   
                    append_hadamard( circ_out, target );
                    append_toffoli( circ_out, controls, target );
                    append_pauli( circ_out, target, pauli_axis::Z, 4u, !tag.adjoint );
                    append_toffoli( circ_out, controls, target );
                    append_pauli( circ_out, control, pauli_axis::Z, 4u, tag.adjoint );
                    append_pauli( circ_out, target, pauli_axis::Z, 4u, tag.adjoint );
                    append_hadamard( circ_out, target );
                }
                else
                {
                    controls.clear();
                    controls.push_back( target );   
                    append_hadamard( circ_out, target );
                    append_toffoli( circ_out, controls, control );
                    append_pauli( circ_out, control, pauli_axis::Z, 4u, !tag.adjoint );
                    append_toffoli( circ_out, controls, control );
                    append_pauli( circ_out, control, pauli_axis::Z, 4u, tag.adjoint );
                    append_pauli( circ_out, target, pauli_axis::Z, 4u, tag.adjoint );
                    append_hadamard( circ_out, target );
                }
            }
            else
                circ_out.append_gate() = gate;
        }
        return circ_out;
    }

    // expand the cnot gates that are not supported by the architecture
    // assume that the corresponding matricies have been set up correctly
    void expand_cnots( circuit& circ_out, const circuit& circ_in ){

        circuit circ_aux;
        copy_metadata( circ_in, circ_aux );
        circ_aux = transform_v_clif(circ_in);
        circ_aux = transform_tof_clif(circ_aux);
        unsigned target, control, moreCnot3 = 0, aux = 0;
        std::vector<unsigned int> new_controls, control2, old_controls;
        
        copy_metadata( circ_in, circ_out );
        
        for ( const auto& gate : circ_aux )
        {
            target = gate.targets().front();
            new_controls.clear();
            new_controls.push_back( target );
            if( !gate.controls().empty() )
            {
                control = gate.controls().front().line();
                old_controls.clear();
                old_controls.push_back( control );
            }
            if ( is_toffoli( gate ) )
            {
                if( gate.controls().empty() ) // a NOT gate
                {
                    circ_out.append_gate() = gate;
                }
                else // CNOT gate
                {
                    moreCnot3 = 0;
                    aux = 0;
                    unsigned cnot3a, cnot3b;
                    for ( auto &p : trans_path[control][target].tpath )
                        if (p.getType() == cnot3 || p.getType() == cnot3i)
                        {
                            if(aux < 1)
                            {
                                cnot3a = p.getA();
                                cnot3b = p.getC();
                            }
                            else
                            {
                                cnot3a = p.getA();
                            }
                            ++aux;
                        }
                    for ( auto &p : trans_path[control][target].tpath ) {
                        switch ( p.getType() )
                        {
                            case cab : //std::cout << "cab" << std::endl;
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                break;
                            case cba : //std::cout << "cba" << std::endl;
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case tab : //std::cout << "tab" << std::endl;
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case tba : //std::cout << "tba" << std::endl;
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                break;
                            case cabi : //std::cout << "cabi" << std::endl;
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case cbai : //std::cout << "cbai" << std::endl;
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                break;
                            case tabi : //std::cout << "tabi" << std::endl;
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getA(), p.getB() );
                                break;
                            case tbai : //std::cout << "tbai" << std::endl;
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case nop : //std::cout << "nop" << std::endl;
                                append_cnot( circ_out, p.getA(), p.getB() );
                                break;
                            case flip : //std::cout << "flip" << std::endl;
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                append_cnot( circ_out, p.getB(), p.getA() );
                                append_hadamard( circ_out, p.getA() );
                                append_hadamard( circ_out, p.getB() );
                                break;
                            case cnot3 :
                                if(moreCnot3 == 0)  // if it is the first cnot3
                                {                   // just append the four cnot gates
                                    append_cnot( circ_out, p.getA(), p.getB() );
                                    append_cnot( circ_out, p.getB(), p.getC() );
                                    append_cnot( circ_out, p.getA(), p.getB() );
                                    append_cnot( circ_out, p.getB(), p.getC() );
                                    ++moreCnot3;    // update the number of cnot3
                                }
                                else    // if it is more than "two arrows"
                                {       // we have to calculate the number of cnots
                                    unsigned c = pow(2,moreCnot3) + pow(2,moreCnot3+1) - 2; // the calculation: 2^n + 2^(n+1) - 2 -> n=number of "arrows" - 1
                                    append_cnot( circ_out, p.getB(), p.getC() );                    // append the cnot   
                                    for (int i = 0, j = circ_out.num_gates()-(c+1); i < c; ++i, ++j)// and copy all the cnots placed before
                                        circ_out.append_gate() = circ_out[j];
                                    append_cnot( circ_out, p.getB(), p.getC() );                    // append again the cnot
                                    ++moreCnot3;
                                }
                                break;
                            case cnot3i :

                                if(moreCnot3 == 0)  // if it is the first cnot3
                                {                   // just append the four cnot gates
                                    append_hadamard( circ_out, cnot3a );
                                    append_hadamard( circ_out, cnot3b );
                                    append_cnot( circ_out, p.getC(), p.getB() );
                                    append_cnot( circ_out, p.getB(), p.getA() );
                                    append_cnot( circ_out, p.getC(), p.getB() );
                                    append_cnot( circ_out, p.getB(), p.getA() );
                                    ++moreCnot3;    // update the number of cnot3
                                }
                                else    // if it is more than "two arrows"
                                {       // we have to calculate the number of cnots
                                    unsigned c = pow(2,moreCnot3) + pow(2,moreCnot3+1) - 2; // the calculation: 2^n + 2^(n+1) - 2 -> n=number of "arrows" - 1
                                    append_cnot( circ_out, p.getB(), p.getA() );                    // append the cnot   
                                    for (int i = 0, j = circ_out.num_gates()-(c+1); i < c; ++i, ++j)// and copy all the cnots placed before
                                        circ_out.append_gate() = circ_out[j];
                                    append_cnot( circ_out, p.getB(), p.getA() );                    // append again the cnot
                                    ++moreCnot3;
                                }
                                if(moreCnot3 == aux)
                                {
                                    append_hadamard( circ_out, cnot3a );
                                    append_hadamard( circ_out, cnot3b );  
                                }
                                break;                                

                                // if(moreCnot3 == 0)  // if it is the first cnot3
                                // {                   // just append the four cnot gates
                                //     append_hadamard( circ_out, p.getA() );
                                //     append_hadamard( circ_out, p.getC() );
                                //     append_cnot( circ_out, p.getB(), p.getA() );
                                //     append_cnot( circ_out, p.getC(), p.getB() );
                                //     append_cnot( circ_out, p.getB(), p.getA() );
                                //     append_cnot( circ_out, p.getC(), p.getB() );
                                //     ++moreCnot3;    // update the number of cnot3
                                // }
                                // else    // if it is more than "two arrows"
                                // {       // we have to calculate the number of cnots
                                //     unsigned c = pow(2,moreCnot3) + pow(2,++moreCnot3) - 2; // the calculation: 2^n + 2^(n+1) - 2 -> n=number of "arrows" - 1
                                //     append_cnot( circ_out, p.getC(), p.getB() );                    // append the cnot   
                                //     for (int i = 0, j = circ_out.num_gates()-(c+1); i < c; ++i, ++j)// and copy all the cnots placed before
                                //         circ_out.append_gate() = circ_out[j];
                                //     append_cnot( circ_out, p.getC(), p.getB() );                    // append again the cnot
                                //     ++moreCnot3;
                                // }
                                // append_hadamard( circ_out, p.getA() );
                                // append_hadamard( circ_out, p.getC() );
                                break;
                            default : std::cout << "ERROR expand_cnots" << std::endl;
                        }
                    }
                }
            }
            else
            {
               circ_out.append_gate() = gate;
            }
        }
    }

 
    unsigned custo;
    matrix mapeamento;
    matrix matrix_custo;
    std::vector<int> vector_costs;

    void extract_matrix(std::vector<int>& m)
    {
        for (int i = 0; i < m.size(); ++i)
        {
            for (int j = 0; j < m.size(); ++j)
            {
                matrix_custo[i][j] = trans_cost[m[i]][m[j]];
                // std::cout << " " << matrix_q[i][j];
            }
            // std::cout << std::endl;
        }
    }

    void all_combinations(int offset, int k, std::vector<int>& qubits, std::vector<int>& combination, matrix& matrix_circuit) 
    {
        // long long int total = 0;
        if (k == 0) 
        {
            // pretty_print(combination);
            unsigned int aux;
            do
            {

                extract_matrix(combination);
                aux = 0;
                for (int i = 0; i < combination.size(); ++i)
                {
                    for (int j = 0; j < combination.size(); ++j)
                    {
                        aux += matrix_custo[i][j]*matrix_circuit[i][j];
                    }
                }
                
                if(aux < custo+10)
                {
                    if(aux < custo)
                        custo = aux;
                    vector_costs.push_back(aux);
                    // mapeamento.clear();
                    mapeamento.push_back(combination);
                    for (int i = 0; i < vector_costs.size(); ++i) 
                    {
                        if(vector_costs[i] > aux+10)
                        {
                            vector_costs.erase(vector_costs.begin()+i);
                            mapeamento.erase(mapeamento.begin()+i);
                            --i;
                        }
                    }
                }

           } while ( std::next_permutation(combination.begin(),combination.begin()+combination.size()) );
            return;
        }
        for (int i = offset; i <= qubits.size() - k; ++i) 
        {
            combination.push_back(qubits[i]);
            all_combinations(i+1, k-1, qubits, combination, matrix_circuit);
            combination.pop_back();
        }
    }

    void clear_matrix(matrix& m)
    {
        for (int i = 0; i < m.size(); ++i)
            m[i].clear();
        m.clear();
    }

    void initialize_matrix(matrix& m, unsigned size)
    {
        std::vector<int> aux;
        for (int i = 0; i < size; ++i)
            aux.push_back(0);
        for (int i = 0; i < size; ++i)
            m.push_back(aux);
    }

    circuit try_all( const circuit& circ_in, bool verbose, bool rm_dup )
    {
        circuit aux, circ_out;
        copy_circuit(circ_in, aux);
        unsigned target, control;
        std::vector<int> qubits;
        std::vector<int> combination;
        matrix matrix_circuit;
        custo = INT_MAX;

        vector_costs.clear(); // clear the vector with the mapping costs
        clear_matrix(mapeamento); // clear the matrix with the mappings
        clear_matrix(matrix_custo); // clear the matrix with the mapping cost

        initialize_matrix(matrix_circuit, circ_in.lines()); // initialize the matrix with zeros
        initialize_matrix(matrix_custo, circ_in.lines()); // initialize the matrix with zeros

        // std::cout << "Number of initial gates: " << circ_in.num_gates() << std::endl;
        // std::cout << "\t" << circ_in.num_gates();

        for ( const auto& gate : circ_in )
        {
            if( !gate.controls().empty() ) // if is not a NOT gate
            {
                target = gate.targets().front();
                control = gate.controls().front().line();
                ++matrix_circuit[control][target];
            }
        }

        for (int i = 0; i < graph_size; ++i) // initialize a vector with the qubits
            qubits.push_back(i);

        // test all the combinations
        all_combinations(0, circ_in.lines(), qubits, combination, matrix_circuit);
        
        // printing all the mappings found in the interval
        if(verbose)
        {
            for (int i = 0; i < mapeamento.size(); ++i)
            {
                std::cout << "Need to add " << vector_costs[i] << " gates ->";
                for (int j = 0; j < mapeamento[i].size(); ++j)
                    std::cout << " " << mapeamento[i][j];
                std::cout << std::endl;
            }    
            std::cout << "Best mapping without optimization: " << custo << " total: " << custo + circ_in.num_gates() << std::endl;
        }
        
        circuit minimo;
        std::vector<int> map_minimo;

        for (int i = 0; i < mapeamento.size(); ++i)
        {
            clear_circuit(aux);
            copy_circuit(circ_in, aux);
            for(unsigned i = circ_in.lines() ; i < graph_size; i++)
                add_line_to_circuit( aux, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
            
            permute_lines(aux, &mapeamento[i][0]);
            expand_cnots( circ_out, aux );
            if(circ_out.num_gates() == custo)
            {
                for (int j = 0; j < mapeamento[i].size(); ++j)
                    std::cout << " " << mapeamento[i][j];             
                std::cout << std::endl;
            }
            if(rm_dup)
                circ_out = remove_dup_gates( circ_out);

            if(i == 0 || circ_out.num_gates() < minimo.num_gates())
            {
                map_minimo.clear();
                for (int j = 0; j < mapeamento[i].size(); ++j) 
                    map_minimo.push_back(mapeamento[i][j]);
                clear_circuit(minimo);
                copy_circuit(circ_out, minimo);
            }
            clear_circuit(circ_out);
        }

        // std::cout << "\t" << minimo.num_gates();

        // bool d = false;
        // for (int i = 0; i < 5; ++i)
        // {
        //     d = false;
        //     for (int j = 0; j < minimo.num_gates(); ++j)
        //     {
        //         if(minimo[j].targets().front() == i)
        //         {
        //             if(minimo[j].controls().empty())
        //             {
        //                 if(d)
        //                 {
        //                     // std::cout << i << " posicao: " << j << std::endl;
        //                     minimo.remove_gate_at(j);
        //                     --j;
        //                     // getchar();
        //                 }
        //                 else
        //                     d = true;
        //             }
        //             else
        //             {
        //                 d = false;
        //             }
        //         }
        //         else if(!minimo[j].controls().empty())
        //             d = false;
        //     }
        // }
        // std::cout << minimo << std::endl;
        // std::cout << "\t" << circ_in.num_gates() << "\t" << minimo.num_gates();        
        std::cout << "Best mapping: " << minimo.num_gates() << " ->";
        for (int i = 0; i < map_minimo.size(); ++i)
            std::cout << " " << map_minimo[i];             
        std::cout << std::endl;
       
        return minimo;
    }

    int get_max_element(matrix& m, unsigned& l, unsigned& c)
    {
        int h = -1;
        for (int i = 0; i < m.size(); ++i)
        {
            for (int j = 0; j < m[i].size(); ++j)
            {
                if(i != j && m[i][j] > h)
                {
                    h = m[i][j];
                    l = i;
                    c = j;
                }
            }
        }
        return h;
    }

    void mark_element(matrix& m, unsigned& l, unsigned& c)
    {
        for (int i = 0; i < m.size(); ++i)
        {
            m[l][i] = -1;
            m[i][l] = -1;
            m[c][i] = -1;
            m[i][c] = -1;
        }
    }

    void print_matrix_circuit(matrix& m)
    {
        std::cout << std::endl;
        for (int i = 0; i < m.size(); ++i)
        {
            for (int j = 0; j < m[i].size(); ++j)
            {
                std::cout << " " << m[i][j];    
            }
            std::cout << std::endl;
        }
    }

    void get_min_element(std::vector<unsigned> a, unsigned& l, unsigned& c)
    {
        int min = INT_MAX;
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
            {
                std::cout << trans_cost[v][w] << " ";
            }
            std::cout << std::endl;
        }
    }

    void mapping( const circuit& circ_in )
    {
        unsigned target, control, max;
        matrix matrix_circuit;
        std::vector< std::pair< unsigned, std::pair< unsigned,unsigned > > > allocation;

        initialize_matrix(matrix_circuit, circ_in.lines());
        
        // Count the number of CNOTs 
        for ( const auto& gate : circ_in )
        {
            if( !gate.controls().empty() ) // if is not a NOT gate
            {
                target = gate.targets().front();
                control = gate.controls().front().line();
                ++matrix_circuit[control][target];
            }
        }

        // Print the matrix of CNOTs
        print_matrix_circuit(matrix_circuit);

        // print_matrix();

        // do the mapping
        for (int i = 0; i < circ_in.lines(); ++i)
        while(1)
        {
            max = get_max_element(matrix_circuit, control, target);
            if(max == -1)
                break;
            matrix_circuit[control][target] = -1;
            allocation.push_back( std::make_pair( max, std::make_pair(control,target) ) );
            // mark_element(matrix_circuit, control, target);
            // std::cout << "=========================" << std::endl;
            // print_matrix(matrix_circuit);
            // allocation.push_back(control);
            // allocation.push_back(target);
            // std::cout << "maior: " << max << " [" << control << "][" << target << "]" << std::endl; 
        }

        for (int i = 0; i < allocation.size(); ++i)
        {
            std::cout << "valor: " << allocation[i].first << " (" << allocation[i].second.first << ", " << allocation[i].second.second << ")" << std::endl;
        }

        for (int i = 0; i < allocation.size(); ++i)
        {
                        
        }
    }
    
}
