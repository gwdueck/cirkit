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

/**
 * @file remove_dup_gates.hpp
 *
 * @brief remove dupplicat adjacent gates that are self inverse
 *
 * @author Gerhard Dueck
 * @since  2.1
 */

#ifndef REMOVE_DUP_GATES_HPP
#define REMOVE_DUP_GATES_HPP

#include <reversible/circuit.hpp>

namespace cirkit
{

circuit remove_dup_gates( const circuit& circ );
bool equal(const gate& g1, const gate& g2 );
bool gates_can_move( const gate& g1, const gate& g2 );
bool gates_do_not_intersect( const gate& g1, const gate& g2 );
bool gates_can_merge( const gate& g1, const gate& g2, gate& res);
bool is_T_gate( const gate& g );
bool is_T_star_gate( const gate& g );
    
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
