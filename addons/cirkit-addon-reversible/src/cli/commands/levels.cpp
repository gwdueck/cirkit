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

#include "levels.hpp"

#include <vector>

#include <boost/format.hpp>

#include <alice/rules.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/remove_dup_gates.hpp>

namespace cirkit
{

levels_command::levels_command( const environment::ptr& env )
  : cirkit_command( env, "Prints the number of levels/nRearranges the gates accordingly." )
{
}

command::rules_t levels_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

bool levels_command::execute()
{
  const auto& circuits = env->store<circuit>();
  const auto& circ = circuits.current();

    circuit result = circ;
    unsigned i = 0, j;
    while(i < result.num_gates() - 1 )
    {
       if( gates_can_move( result[i], result[i+1]))
       {
           std::cout << i << " can move\n";
       }
       else{
           std::cout << i << " can NOT move\n";
       }
        i++;
    }
/*
  for ( const auto& g : circ )
  {
      
    auto key = gate_controls.size() - 1u;


    if ( is_toffoli( g ) )
    {
      key = 0u;
    }
    else if ( is_fredkin( g ) )
    {
      key = 1u;
    }
    else if ( is_stg( g ) )
    {
      key = 2u;
    }
    else if ( is_pauli( g ) )
    {
      key = 3u;
    }
    else if ( is_hadamard( g ) )
    {
      key = 4u;
    }

    auto& v = gate_controls[key];
    const auto num_controls = g.controls().size();
    if ( v.size() <= num_controls )
    {
      v.resize( num_controls + 1u );
    }
    v[num_controls]++;

    if ( max_controls < num_controls )
    {
      max_controls = num_controls;
    }

  }
*/
  return true;
}

command::log_opt_t levels_command::log() const
{
  return boost::none;
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
