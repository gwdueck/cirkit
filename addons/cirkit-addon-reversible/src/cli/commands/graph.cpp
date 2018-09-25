/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
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

#include "graph.hpp"
#include <alice/rules.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/circuit.hpp>

#include <fstream>
#include <core/utils/string_utils.hpp>
#include <reversible/IBMgraph.hpp>



namespace cirkit
{

using boost::program_options::value;

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

graph_command::graph_command( const environment::ptr& env )
    : cirkit_command( env, "Manipulate the graph for IBMs architectures" )
{
    opts.add_options()
    ( "read,r", value( &filename ), "read graph from a file" )
    ( "create,c",  "create the transformation matrix" )
    ( "verbose,v", "verbose mode")
    ( "print,p",  "print current graph" )
    ( "matrix_cost,q",  "Print matrix cost and transformations" )
    ( "transform,t", "transform non supported cnot gates")
    ( "rm_dup,m",  "Remove duplicate gates" )
    ( "delete,d",  "delete current graph" )
    ( "file,w", value( &filename ), "Write graph, matrix and transformations to a file" )
    ( "from_file,f", value( &filename ), "Read matrix and transformations from a file" )
    ( "mapping,x", "Realize the mapping for a current circuit")
    ;
    add_new_option();
}

bool graph_command::execute()
{
    bool verbose = false;
    

    if( is_set( "verbose" ) ){
        verbose = true;
    }
    if( is_set( "read" ) ){
        read_graph( filename );
    }
    if( is_set( "print" ) ){
        print_graph( );
    }
    if( is_set( "matrix_cost" ) )
    {
        print_matrix( );
    }
    if( is_set( "file" ) )
    {
        write_to_file( filename );
    }
    if( is_set( "from_file" ) )
    {
        read_from_file( filename );
    }
    if( is_set( "create" ) ){
        create_trans( verbose );
    }
    if( is_set( "delete" ) ){
        delete_graph( );
    }
    if( is_set( "mapping" ) ){
        if( env->store<circuit>().current_index() < 0 ){
            std::cout << "no current circuit available" << std::endl;
            return true;
        }
        auto& circuits = env->store<circuit>();
        circuit circ_out, circ_in = circuits.current();
        the_mapping( circ_out, circ_in);
        if ( is_set( "new" ) )
        {
            circuits.extend();
        }
        circuits.current() = circ_out;
    }

    if( is_set( "transform" ) ){
        if( env->store<circuit>().current_index() < 0 ){
            std::cout << "no current circuit available" << std::endl;
            return true;
        }
        auto& circuits = env->store<circuit>();
        circuit circ_working = circuits.current();
        circuit circ_result;
        expand_cnots( circ_result, circ_working );
        if ( is_set( "new" ) )
        {
            circuits.extend();
        }
        if ( is_set( "rm_dup" ) )
        {
            circ_result = remove_dup_gates( circ_result );
        }
        circuits.current() = circ_result;
    }
    return true;
}

command::log_opt_t graph_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}

// command::rules_t graph_command::validity_rules() const
// {
//   return {has_store_element<circuit>( env )};
// }


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
