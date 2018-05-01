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
 * @file ibm_helper.hpp
 *
 * @brief some help functions for IBM QX
 *
 * @author Gerhard Dueck
 * @since  2.3
 */

#ifndef IBM_HELPER_HPP
#define IBM_HELPER_HPP

#include <reversible/circuit.hpp>

/* NOTE (Apr 22, 2018): added new notation
 map methods for CNOT gates are as follows:
 0 - no mapping possible (eg CNOT(1,1))
 1 - CNOT gate Exists
 2 - target and controls must be interchanged -- FLIP(c,t)
 3 - map target to qubit 2 given CNOT(c,2) and CNOT(t,2) -- TBA(t,2)
 4 - map control to qubit 2, given CNOT(2,c) and CNOT(2,t) -- CBA(c,2)
 5 - map target to qubit 2 and interchange control and qubit 2, given CNOT(c,2) and CNOT(2,t)  -- TAB(t,2) FLIP(c,2)
 6 - map control to qubit 2, given CNOT(c,2) and CNOT(2,t) -- CAB(c,2)
 */

int static const map_method_qx2[5][5] ={{0,1,1,3,3},
                                        {2,0,1,3,3},
                                        {2,2,0,2,2},
                                        {3,3,1,0,1},
                                        {3,3,1,2,0}};
int static const map_method_qx4[5][5] ={{0,2,2,5,4},
                                        {1,0,2,5,4},
                                        {1,1,0,2,1},
                                        {6,6,1,0,1},
                                        {4,4,2,2,0}};

namespace cirkit
{

void permute_lines( circuit& circ , int perm[]);

circuit transform_to_IBMQ( const circuit& circ, const int map_method[5][5], bool templ );

int levels(const circuit& circ, circuit& result );
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
