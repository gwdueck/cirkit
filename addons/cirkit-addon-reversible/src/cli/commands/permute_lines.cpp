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

#include "permute_lines.hpp"

#include <cmath>
#include <vector>

#include <alice/rules.hpp>

#include <core/utils/string_utils.hpp>
#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <cli/reversible_stores.hpp>

using namespace boost::program_options;

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

permute_lines_command::permute_lines_command( const environment::ptr& env )
  : cirkit_command( env, "Permute the lines of a circuit" )
{
  opts.add_options()
    ( "permutation,p", value( &permutation ), "Create permute_lines from permutation (starts with 0, space separated)" )
    ( "new,n",                                "Add a new entry to the store; if not set, the current entry is overriden" )
    ;
}

command::rules_t permute_lines_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

bool permute_lines_command::execute()
{
    auto& circuits = env->store<circuit>();
    std::vector<int> perm;
    parse_string_list( perm, permutation );
    circuit circ_working = circuits.current();
    permute_lines( circ_working , &perm[0] );
    if ( is_set( "new" ) )
    {
        circuits.extend();
    }
    circuits.current() = circ_working;
    return true;
}

    
    // permute the line in the circuit
    // assume it is a qc circuit
    void permute_lines( circuit& circ , int perm[])
    {
        unsigned target, control;
        for ( auto& gate : circ )
        {
            assert( gate.targets().size() == 1 );
            target = gate.targets().front();
            gate.remove_target(target);
            gate.add_target(perm[target]);
            if( !gate.controls().empty() )
            {
                assert( gate.controls().size() == 1 );
                control = gate.controls().front().line();
                gate.remove_control( make_var(control) );
                gate.add_control( make_var(perm[control]) );
                
            }
            
        }
        
    }
    
}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
