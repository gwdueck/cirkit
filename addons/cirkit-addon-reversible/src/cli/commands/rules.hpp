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
 * @file alex.hpp
 *
 * @brief Test
 *
 * @author A.G.A.
 * @2
 */

#ifndef RULES_HPP
#define RULES_HPP

#include <reversible/circuit.hpp>
#include <reversible/functions/add_gates.hpp>

#include <cli/cirkit_command.hpp>

namespace cirkit
{
	
	bool verify_rule_Done ( const gate& ga, const gate& gb );
	void apply_rule_Done ( circuit& circ, unsigned itGateIndex, unsigned nextGateIndex );
	bool verify_rule_Dthree ( const gate& ga, const gate& gb );
	void apply_rule_Dthree ( circuit& circ, circuit::const_iterator& itGate, circuit::const_iterator& nextGate );
	bool verify_rule_Dfour ( const gate& ga, const gate& gb );
	void apply_rule_Dfour ( circuit& circ, circuit::const_iterator& itGate, circuit::const_iterator& nextGate );
	void apply_rule_Rfive( gate& ga, gate& gb );
	void remove_line_control_target( gate& control, const gate& target );
	bool line_control_target( const gate& control, const gate& target );
	bool different_polarity_controls( const gate& ga, const gate& gb );
	bool single_control( const gate& ga, const gate& gb );
	bool controls_same_line( const gate& ga, const gate& gb );
	bool targets_same_line( const gate& ga, const gate& gb );
	void swap_gates( circuit::const_iterator itGate, circuit::const_iterator nextGate );
	bool verify_rule_Dtwo (  gate ga, gate gb );
	void apply_rule_Dtwo ( circuit::const_iterator itGate, circuit::const_iterator nextGate );
	bool verify_rule_Rfour ( const gate& ga, const gate& gb );
	bool verify_rule_Dsix ( const gate& ga, const gate& gb );
	bool verify_rule_Dseven ( const gate& ga, const gate& gb );
	void apply_rule_Dseven ( circuit& circ, circuit::const_iterator& itGate, circuit::const_iterator& nextGate );
}

#endif

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
