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
    ( "delete,d",  "delete current graph" )
    ( "print,p",  "print current graph" )
    ( "create,c",  "create the transformation matrix" )
    ( "verbose,v", "verbose mode")
    ( "transform,t", "transform non supported cnot gates")
    ( "rm_dup,r",  "Remove duplicate gates" )
    ;
    add_new_option();
}

    
bool graph_command::execute()
{
    bool verbose = false;
    
    
    if( is_set( "verbose" ) )
    {
        verbose = true;
    }

    if( is_set( "read" ) )
    {
        read_graph( filename );
    }
    if( is_set( "print" ) ){
        print_graph( );
    }
    if( is_set( "create" ) ){
        create_trans( verbose );
    }
    if( is_set( "delete" ) ){
        delete_graph( );
    }
    if( is_set( "transform" ) ){
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



}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
