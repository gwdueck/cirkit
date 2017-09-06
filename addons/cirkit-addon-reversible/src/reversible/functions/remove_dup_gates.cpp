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

#include "remove_dup_gates.hpp"

#include <reversible/gate.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/pauli_tags.hpp>

namespace cirkit
{

circuit remove_dup_gates( const circuit& circ )
{
    circuit result = circ;
    unsigned i = 0, j;
    while(i < result.num_gates() - 1 )
    {
        j = i + 1;
        bool done = false;
        bool incr_i = true;
        while( ( !done ) && ( j < result.num_gates() ) )
        {
            if( equal( result[i], result[j] ) )
            {
                result.remove_gate_at(j);
                result.remove_gate_at(i);
                done = true;
                i = 0; // overkill, but to safe 
                incr_i = false;
            }
            if(!done && gates_can_move( result[i], result[j]) )
            {
                j++;
            }
            else{
                done = true;
            }
        }
        if ( incr_i )
        {
            i++;
        }
    }
    return result;
}

// check if two gates are equal
// only consderred a few self-inverse gates
bool equal(const gate& g1, const gate& g2 )
{
    if( !same_type( g1, g2 )){
        return false;
    }
    if( !(is_toffoli( g1 ) || is_hadamard( g1 ) ) )
    {
        return false;
    }
    if(( g1.targets().size() != g2.targets().size() ) ||
       ( g1.controls().size() != g2.controls().size() ))
    {
        return false;
    }
    if( ( g1.targets().size() != 1 ) ||
       ( g1.targets().front() != g2.targets().front()) )
    {
        return false;
    }
    if( ( g1.controls().size() == 1 ) &&
       ( g1.controls().front().line() != g2.controls().front().line() ) )
    {
        return false;
    }
    
    return true;
}
 
// check if the two gates may be inter changed
// WARNING: only written for Clifford+T gates
bool gates_can_move( const gate& g1, const gate& g2 )
{
    unsigned target_g2 = g2.targets().front();
    // only move hadamard if it does not intersect with control nor target
    if ( is_hadamard( g1 ) )
    {
        if ( target_g2 == g1.targets().front() ){
            return false;
        }
        if ( g2.controls().empty() )
        {
            return true;
        }
        else{
            unsigned control_g2 = g2.controls().front().line();
            if ( g1.targets().front() == control_g2 )
            {
                return false;
            }
        }
        return true;
    }
    return false;
}
    
}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
