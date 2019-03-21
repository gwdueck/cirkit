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
#include <algorithm>

#include <boost/format.hpp>

#include <alice/rules.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/gate.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/io/print_circuit.hpp>

namespace cirkit
{

levels_command::levels_command( const environment::ptr& env )
  : cirkit_command( env, "Prints the number of levels\nRearranges the gates accordingly." )
{
    add_new_option();
}

command::rules_t levels_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}
/* This function will calculate the number of levels in the circuit.
 The algorithm is as follows:
 for each gate
     do
         if the gate can join the level to the left
             remember the position
         move gate to the left
     until gate cannot be moved
     put gate in pos (added to that level)
 
 */
    
bool levels_command::execute()
{
  auto& circuits = env->store<circuit>();
    
  std::vector<int> glevel; /* level of the gates */
    
    circuit result = circuits.current();
    unsigned i = 1, max_lev = 1;
    int pos, k, j;
    glevel.push_back( max_lev );
    bool blocked = false;
    while (i < result.num_gates() )
    {
        //std::cout << " i = " << i << " ";
        pos = -1;
        j = i - 1; /* index of the end of a level */
        do
        {
            //std::cout << " j = " << j << std::endl;
            bool can_join = true; /* can i join the previous level? */
            k = j; /* iterate over the level */
            while( ( k >= 0 ) && ( glevel[k] == glevel[j] ) && can_join )
            {
                can_join = gates_do_not_intersect( result[k], result[i] );
                k--;
            }
            if( can_join )
            {
                pos = j + 1;
            }
            k = j; /* iterate over the level */
            while ( ( k >= 0 ) && ( glevel[k] == glevel[j] ) )
            {
                blocked = !gates_can_move( result[k], result[i] );
                k--;
            }
            j = k;
            
        } while ( !blocked && ( j >= 0 ) );
    
        if ( pos >= 0 )
        {
            if ( pos < (int) i )
            {
                //std::cout << "move to " << pos << "\n";
                result.insert_gate( pos ) = result[i];
                result.remove_gate_at( i + 1 );
                glevel.insert( glevel.begin()+pos , glevel[pos-1] );
            }
            else
            {
                glevel.push_back( max_lev );
            }
        }
        else
        {
            max_lev++;
            glevel.push_back( max_lev );
        }
/*
        for(int tt = 0; tt < (int) glevel.size(); tt++ )
            std::cout << glevel[tt] << " ";
        std::cout << std::endl;
        print_circuit( subcircuit( result, 0, i + 1 ) );
  */
        i++;
    }
    // std::cout << " " << max_lev;
    
    std::cout << "The circuit has " << max_lev << " levels" << std::endl;
    // std::cout << " levels: " << max_lev;
    if ( is_set( "new" ) )
    {
        circuits.extend();
    }
    circuits.current() = result;

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
