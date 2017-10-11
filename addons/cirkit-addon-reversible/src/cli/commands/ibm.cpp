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

#include <alice/rules.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/program_options.hpp>
#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <cli/reversible_stores.hpp>
#include <cli/commands/permute_lines.hpp>
#include <reversible/optimization/simplify.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/io/write_qc.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/variable.hpp>

/* map methods for CNOT gates are as follows:
 0 - no papping possible (eg CNOT(1,1)
 1 - CNOT gate Exists
 2 - target and controls must be interchanged
 3 - map target to qubit 2
 4 - map control to qubit 2
 5 - map target to qubit 2 and interchange control and qubit 2
 */
int static const map_method_qx2[5][5] = {{0,1,1,3,3}, {2,0,1,3,3}, {2,2,0,2,2}, {3,3,1,0,1}, {3,3,1,2,0}};
int static const map_method_qx4[5][5] = {{0,2,2,5,4}, {1,0,2,5,4}, {1,1,0,2,1}, {4,4,1,0,1}, {4,4,2,2,0}};
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
    
    if( !is_set( "all_perm" ) )
    {
        if ( is_set( "ibm_qx4" ) )
        {
            circ_IBM = transform_to_IBMQ( circ_working, map_method_qx4 );
        }
        else
        {
            circ_IBM = transform_to_IBMQ( circ_working, map_method_qx2 );
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
    else
    {
        int perm[5] = {0, 1, 2, 3, 4}, inv_perm[5];
        do
        {
            permute_lines( circ_working , perm );
            if ( is_set( "ibm_qx4" ) )
            {
                circ_IBM = transform_to_IBMQ( circ_working, map_method_qx4 );
            }
            else
            {
                circ_IBM = transform_to_IBMQ( circ_working, map_method_qx2 );
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
            if( is_set( "verbose" ) )
            {
                for( int i = 0; i < 5; i++ )
                {
                    std::cout << perm[i] << " ";
                }
                std::cout << "gates = " << circ_IBM.num_gates() << std::endl;
            }

            // undo the permutation
            for( int i = 0; i < 5; i++ )
            {
                inv_perm[perm[i]] = i;
            }
            permute_lines( circ_working , inv_perm );
        } while ( std::next_permutation(perm,perm+5) );
    }
    return true;
}

command::log_opt_t ibm_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}

// transform a Clifford+T circuit to be IBM compliant
circuit transform_to_IBMQ( const circuit& circ, const int map_method[5][5] )
{
    circuit circ_IBM;
    unsigned target, control;
    std::vector<unsigned int> new_controls, control2;
    control2.push_back( 2u ); /* use line 2 as control */
    
    copy_metadata( circ, circ_IBM );
    
    // all IBM circuits have exactly 5 lines
    for(unsigned i = circ.lines()+1 ; i <= 5u; i++)
    {
        add_line_to_circuit( circ_IBM, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
    }

    
    // iterate through the gates
    for ( const auto& gate : circ )
    {
        
        target = gate.targets().front();
        new_controls.clear();
        new_controls.push_back( target );
        if( !gate.controls().empty() )
        {
            control = gate.controls().front().line();
        }
        
        if ( is_toffoli( gate ) )
        {
            if( gate.controls().empty() ) // a NOT gate
            {
                append_toffoli( circ_IBM, gate.controls(), target );
            }
            else // CNOT gate
            {
                
                
                switch ( map_method[control][target] )
                {
                    case 1:
                        append_toffoli( circ_IBM, gate.controls(), target );
                        break;
                        
                    case 2 : // invert CNOT
                        
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, target );
                        append_toffoli( circ_IBM, new_controls, control );
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, target );
                        break;
                    case 3 : // swap target with 2
                        append_toffoli( circ_IBM, new_controls, 2u );
                        append_hadamard( circ_IBM, 2u );
                        append_hadamard( circ_IBM, target );
                        append_toffoli( circ_IBM, new_controls, 2u );
                        append_hadamard( circ_IBM, 2u );
                        
                        append_toffoli( circ_IBM, gate.controls(), 2u );
                        
                        append_hadamard( circ_IBM, 2u );
                        append_toffoli( circ_IBM, new_controls, 2u );
                        append_hadamard( circ_IBM, 2u );
                        append_hadamard( circ_IBM, target );
                        append_toffoli( circ_IBM, new_controls, 2u );
                        break;
                    case 4 : // swap control with 2
                        append_toffoli( circ_IBM, control2, control );
                        append_hadamard( circ_IBM, 2u );
                        append_hadamard( circ_IBM, control );
                        append_toffoli( circ_IBM, control2, control );
                        append_hadamard( circ_IBM, 2u );
                        
                        append_toffoli( circ_IBM, control2, target );
                        
                        append_hadamard( circ_IBM, 2u );
                        append_toffoli( circ_IBM, control2, control );
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, 2u );
                        append_toffoli( circ_IBM, control2, control );
                        break;
                    case 5: // swap target with qubit 2 and interchange control and qubit 2
                        std::cout << "case 5\n";
                        append_toffoli( circ_IBM, new_controls, 2u );
                        append_hadamard( circ_IBM, 2u );
                        append_hadamard( circ_IBM, target );
                        append_toffoli( circ_IBM, new_controls, 2u );
                        
                        append_hadamard( circ_IBM, control );
                        append_toffoli( circ_IBM, control2, control );
                        append_hadamard( circ_IBM, control );
                        
                        append_toffoli( circ_IBM, new_controls, 2u );
                        append_hadamard( circ_IBM, 2u );
                        append_hadamard( circ_IBM, target );
                        append_toffoli( circ_IBM, new_controls, 2u );
                        break;
                        
                }
            }
        }
        else if ( is_pauli( gate ) )
        {
            const auto& tag = boost::any_cast<pauli_tag>( gate.type() );
            append_pauli( circ_IBM, target, tag.axis, tag.root, tag.adjoint );
        }
        else if ( is_hadamard( gate ) )
        {
            append_hadamard( circ_IBM, target );
        }
        else
        {
            assert( false );
        }
    }

    
    return circ_IBM;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
