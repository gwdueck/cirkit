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
 * @file clifford_templates.hpp
 *
 * @brief structures for Clifford+T templates
 *
 * @author Gerhard Dueck
 * @since  2.1
 */

#ifndef CLIFFORD_TEMPLATES_HPP
#define CLIFFORD_TEMPLATES_HPP

#include <reversible/target_tags.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/circuit.hpp>
#include <reversible/functions/remove_dup_gates.hpp>


namespace cirkit
{

bool is_CNOT_gate( const gate& g );

enum Cliff_Gate_Type { H, T, Ts, S, Ss, Z, Y, RZ, V, Vs, X, CNOT };

const static std::string gate_name[] = { "h", "t", "T", "s", "S", "z", "y", "r", "v", "V", "x", "c" };

static bool ( *is_Gate[] )( const gate& g ) = {
	&is_hadamard,
	&is_T_gate,
	&is_T_star_gate,
	&is_S_gate,
	&is_S_star_gate,
	&is_Z_gate,
	&is_Y_gate,
	&is_RZ_gate,
	&is_V_gate,
	&is_V_star_gate,
	&is_X_gate,
	&is_CNOT_gate };

static std::map<char, Cliff_Gate_Type> cliff_map = {
	{'h', H},
	{'t', T}, 
	{'T', Ts},
	{'s', S},
	{'S', Ss},
	{'z', Z},
	{'y', Y},
	{'r', RZ},
	{'v', V},
	{'V', Vs},
	{'x', X},
	{'c', CNOT}
};

class Cliff_Gate{
public:
	Cliff_Gate_Type gtype;
	int target, control;
	void convert_gate( const gate &g );
};

class Clifford_Template{
public:
	std::vector<Cliff_Gate> gates_matched;
	std::vector<Cliff_Gate> gates_replaced;
	int num_qubits;
	void print();
	void read( std::ifstream &infile );
	void convert_circ( circuit &circ );
	void clear();
};

static std::vector<Clifford_Template> cliff_templates;

}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
