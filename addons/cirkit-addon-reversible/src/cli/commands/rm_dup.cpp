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

#include "rm_dup.hpp"

#include <fstream>
#include <algorithm>

#include <alice/rules.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/program_options.hpp>
#include <reversible/circuit.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/gate.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/pauli_tags.hpp>

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

rm_dup_command::rm_dup_command( const environment::ptr& env )
    : cirkit_command( env, "rm_dup circuit" )
{
    opts.add_options()
    ;
  add_new_option();
}


command::rules_t rm_dup_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

bool rm_dup_command::execute()
{
	auto& circuits = env->store<circuit>();
    circuit circ = circuits.current();
    circuit circ_rm;
    circ_rm = remove_dup_gates(circ);
    gate g;
    unsigned i = 0, j;
    while(i < circ_rm.num_gates() - 1 )
    {
        j = i + 1;
        bool done = false;
        bool incr_i = true;
        while( ( !done ) && ( j < circ_rm.num_gates() ) )
        {
            // if ( is_toffoli( gate ) )
            // {
            //     prepend_toffoli( circ_invert, gate.controls(), gate.targets().front() );
            // }
            if ( is_hadamard( circ_rm[i] ) && is_hadamard( circ_rm[i+1] ) &&
                 is_toffoli( circ_rm[i+2] ) &&  
                 is_hadamard( circ_rm[i+3] ) && is_hadamard( circ_rm[i+4]))
            {
                std::cout << "ASdsaDsaD" << std::endl;
                // if( circ_rm[i+2].controls() )
                // {

                // }
                //prepend_hadamard( circ_invert,  gate.targets().front() );
            }
            // else if ( is_pauli( gate ) )
            // {
            //     const auto& tag = boost::any_cast<pauli_tag>( gate.type() );
            //     prepend_pauli( circ_invert,  gate.targets().front(), tag.axis, tag.root, !tag.adjoint );
            // }

            // if( is_inverse( circ_rm[i], circ_rm[j] ) )
            // {
            //     circ_rm.remove_gate_at(j);
            //     circ_rm.remove_gate_at(i);
            //     done = true;
            //     i = 0; // overkill, but to be safe
            //     incr_i = false;
            // }
            // if ( !done && gates_can_merge( circ_rm[i], circ_rm[j], g) )
            // {
            //     circ_rm.remove_gate_at(j);
            //     circ_rm[i] = g;
            //     done = true;
            //     i = 0; // overkill, but to be safe
            //     incr_i = false;
            // }
            // if(!done && gates_can_move( circ_rm[i], circ_rm[j]) )
            // {
            //     j++;
            // }
            // else{
            //     done = true;
            // }
            j++;
        }
        i++;
        if ( incr_i )
        {
            i++;
        }
    }
    if ( is_set( "new" ) )
    {
        circuits.extend();    
    }
    circuits.current() = circ_rm;
    return true;
}

command::log_opt_t rm_dup_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
