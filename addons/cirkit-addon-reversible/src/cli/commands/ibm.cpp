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

#include "ibm.hpp"

#include <fstream>
#include <algorithm>
#include <limits>

#include <alice/rules.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/program_options.hpp>
#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/optimization/simplify.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/rotation_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/functions/ibm_helper.hpp>
#include <reversible/io/write_qc.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/variable.hpp>
#include <cli/commands/ibm.hpp>
#include <reversible/functions/ibm_helper.hpp>



namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

ibm_command::ibm_command( const environment::ptr& env )
    : cirkit_command( env, "Translate Clifford+T circuits to IBM Q\nArchitecture: qx2 (default) or qx4" )
{
    opts.add_options()
    ( "all_perm,a",  "Try all permutations" )
    ( "rm_dup,r",  "Remove duplicate gates" )
    ( "ibm_qx4,4", "The IBM Qx4 is the target")
    ( "verbose,v",  "verbose" )
    ( "template,t", "use template transformations--instead of swap based")
    ;
  add_new_option();
}


command::rules_t ibm_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

    
bool ibm_command::execute()
{

    auto& circuits = env->store<circuit>();
    circuit circ_working = circuits.current();
    circuit circ_IBM;
    unsigned start = circ_working.lines()+1;
    for(unsigned i = start ; i <= 5u; i++)
    {
        add_line_to_circuit( circ_working, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
    }
    std::cout << " " << circ_working.num_gates();
    if( !is_set( "all_perm" ) )
    {
        if ( is_set( "ibm_qx4" ) )
        {
            circ_IBM = transform_to_IBMQ( circ_working, map_method_qx4, is_set( "template" ) );
        }
        else
        {
            circ_IBM = transform_to_IBMQ( circ_working, map_method_qx2, is_set( "template" ) );
        }
        
        if ( is_set( "new" ) )
        {
            circuits.extend();
        }
        if ( is_set( "rm_dup" ) )
        {
            circ_IBM = remove_dup_gates( circ_IBM);
        }
        circuits.current() = circ_IBM;
    }
    else // go through all permutations
    {
        int perm[5] = {0, 1, 2, 3, 4}, inv_perm[5], best_perm[5] = {0, 1, 2, 3, 4};
        unsigned best_cost = UINT_MAX;
        circuit circ_best;
        do
        {
            permute_lines( circ_working , perm );
            if ( is_set( "ibm_qx4" ) )
            {
                circ_IBM = transform_to_IBMQ( circ_working, map_method_qx4, is_set( "template" ) );
            }
            else
            {
                circ_IBM = transform_to_IBMQ( circ_working, map_method_qx2, is_set( "template" ) );
            }
            if ( is_set( "new" ) )
            {
                circuits.extend();
            }
            
            if ( is_set( "rm_dup" ) )
            {
                circ_IBM = remove_dup_gates( circ_IBM );
            }
            circuits.current() = circ_IBM;
            if( is_set( "verbose" ) )
            {
                for( int i = 0; i < 5; i++ )
                {
                    std::cout << perm[i] << " ";
                }
                std::cout << "gates = " << circ_IBM.num_gates() << std::endl;
            }
            
            if( best_cost > circ_IBM.num_gates() )
            {
  //              std::cout << "new best_cost = " << circ_IBM.num_gates() << "\n";
                best_cost = circ_IBM.num_gates();
                circ_best = circ_IBM;
                for( int i = 0; i < 5; i++ )
                {
                    best_perm[i] = perm[i];
                }
            }

            // undo the permutation
            for( int i = 0; i < 5; i++ )
            {
                inv_perm[perm[i]] = i;
            }
            permute_lines( circ_working , inv_perm );
        } while ( std::next_permutation(perm,perm+5) );
        if ( is_set( "new" ) )
        {
            circuits.extend();
        }
        circuits.current() = circ_best;
        // std::cout << "best permutation = ";
        // for( int i = 0; i < 5; i++ )
        // {
        //     std::cout << best_perm[i] << " ";
        // }
        // std::cout << "gates = " << best_cost << std::endl;
        std::cout << " " << best_cost;
    }
    return true;
}

command::log_opt_t ibm_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}



}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
