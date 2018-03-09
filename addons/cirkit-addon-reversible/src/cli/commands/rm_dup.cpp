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
#include <reversible/functions/add_gates.hpp>
#include <reversible/io/print_circuit.hpp>

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
    circuit circ_rm = circuits.current();
    //circuit circ_rm;
    //circ_rm = remove_dup_gates(circ);
    unsigned i = 0, j, target;
    std::vector<unsigned int> controls;
    while(i < circ_rm.num_gates())
    {
        if ( i+5 < circ_rm.num_gates() && 
            is_hadamard( circ_rm[i] ) &&  
            is_hadamard( circ_rm[i+1] ) &&
            is_toffoli( circ_rm[i+2] ) &&
            is_hadamard( circ_rm[i+3] ) && 
            is_hadamard( circ_rm[i+4]))
        {
            if((circ_rm[i+2].targets().front() == circ_rm[i].targets().front()) ||
                circ_rm[i+2].targets().front() == circ_rm[i+1].targets().front())
            {
                target = circ_rm[i+2].controls().front().line();
                controls.push_back(circ_rm[i+2].targets().front());
                insert_toffoli( circ_rm, i,controls, target );
                for (int j = 0; j < 5; ++j)
                    circ_rm.remove_gate_at(i);              
            }
        }
        else if ( i+10 < circ_rm.num_gates() &&
                is_toffoli( circ_rm[i] ) &&
                is_hadamard( circ_rm[i+1] ) &&
                is_hadamard( circ_rm[i+2]) &&
                is_toffoli( circ_rm[i+3]) &&
                is_hadamard( circ_rm[i+4]) &&
                is_toffoli( circ_rm[i+5] ) &&
                is_hadamard( circ_rm[i+6] ) &&
                is_toffoli( circ_rm[i+7] ) &&
                is_hadamard( circ_rm[i+8] ) &&
                is_hadamard( circ_rm[i+9] ) &&
                is_toffoli( circ_rm[i+10] ))
        {
            if(circ_rm[i+5].targets().front() == circ_rm[i+4].targets().front())
            {
                target = circ_rm[i].controls().front().line();
                controls.push_back(circ_rm[i+5].controls().front().line());
            }
            else
            {
                target = circ_rm[i+5].targets().front();
                controls.push_back(circ_rm[i].targets().front());
            }
            insert_toffoli( circ_rm, i,controls, target );
            for (int j = 0; j < 11; ++j)
                circ_rm.remove_gate_at(i);                  
        }
        controls.clear();
        ++i;
    }
    circ_rm = remove_dup_gates(circ_rm);
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
