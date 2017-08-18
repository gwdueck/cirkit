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
    unsigned i = 1;
    while(i < result.num_gates())
    {
        if( equal(result[i-1], result[i]) )
        {
            result.remove_gate_at(i);
            result.remove_gate_at(i-1);
            i--;
        }
        else
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
    
}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
