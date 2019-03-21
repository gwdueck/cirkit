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

#include "ibm_helper.hpp"

#include <boost/lexical_cast.hpp>

#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/rotation_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/functions/remove_dup_gates.hpp>

namespace cirkit
{
// permute the line in the circuit
// assume it is a qc circuit
// void permute_lines( circuit& circ , int perm[])
// {
//     unsigned target, control;
//     for ( auto& gate : circ )
//     {
//         assert( gate.targets().size() == 1 );
//         target = gate.targets().front();
//         gate.remove_target(target);
//         gate.add_target(perm[target]);
//         if( !gate.controls().empty() )
//         {
//             assert( gate.controls().size() == 1 );
//             control = gate.controls().front().line();
//             gate.remove_control( make_var(control) );
//             gate.add_control( make_var(perm[control]) );
//         }
//     }
// }

void permute_lines( circuit& circ , int perm[])
{
    unsigned target;
    std::vector<unsigned> control;
    std::vector<bool> polarity;

    for ( auto& gate : circ )
    {
        assert( gate.targets().size() == 1 );
        target = gate.targets().front();
        gate.remove_target(target);
        gate.add_target(perm[target]);
        if( !gate.controls().empty() )
        {
            // assert( gate.controls().size() == 1 );
            control.clear();
            polarity.clear();
            for ( auto& c : gate.controls() )
            {
                control.push_back( c.line() );
                polarity.push_back( c.polarity() );
            }
            for (int i = 0; i < control.size(); ++i)
                gate.remove_control( make_var( control[i], polarity[i] ) );
            for (int i = 0; i < control.size(); ++i)
                gate.add_control( make_var( perm[control[i]],  polarity[i]) );
        }
    }
}

// transform a Clifford+T circuit to be IBM compliant
circuit transform_to_IBMQ( const circuit& circ, const int map_method[5][5], bool templ )
{
    circuit circ_IBM;
    unsigned target, control;
    std::vector<unsigned int> new_controls, control2, old_controls;
    control2.push_back( 2u ); /* use line 2 as control */
    
    copy_metadata( circ, circ_IBM );
    
    // all IBM circuits have exactly 5 lines
    for(unsigned i = circ.lines()+1 ; i <= 5u; i++)
    {
        add_line_to_circuit( circ_IBM, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
    }
    
    
    // iterate through the gates
    for ( const auto& gate : circ )
    {
        
        target = gate.targets().front();
        new_controls.clear();
        new_controls.push_back( target );
        if( !gate.controls().empty() )
        {
            control = gate.controls().front().line();
            old_controls.clear();
            old_controls.push_back( control );
        }
        
        if ( is_toffoli( gate ) )
        {
            if( gate.controls().empty() ) // a NOT gate
            {
                append_toffoli( circ_IBM, gate.controls(), target );
            }
            else if ( gate.controls().size() == 1 ) // CNOT gate
            {
                //std::cout << "CNOT case " << map_method[control][target] << "\n";
                
                switch ( map_method[control][target] )
                {
                    case 1:
                        append_toffoli( circ_IBM, gate.controls(), target );
                        break;
                        
                    case 2 : // invert CNOT
                        
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, target );
                        append_toffoli( circ_IBM, new_controls, control );
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, target );
                        break;
                    case 3 : // swap target with 2
                        if( !templ )
                        {
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            
                            append_toffoli( circ_IBM, gate.controls(), 2u );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                        }
                        else // use the "template" transformation
                        {
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, 2u );
                        }
                        break;
                    case 4 : // map control to qubit 2, given CNOT(2,c) -- CBA(c,2)
                        if( !templ )
                        {
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            
                            append_toffoli( circ_IBM, control2, target );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, control );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                        }
                        else
                        {
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, target );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, target );
                            
                            
                        }
                        break;
                    case 5: // map control to qubit 2, given CNOT(c,2) -- CAB(c,2)
                         if( !templ )
                        {
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            
                            append_toffoli( circ_IBM, control2, target );
                            
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            
                        }
                        else
                        {
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_toffoli( circ_IBM, control2, target );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_toffoli( circ_IBM, control2, target );
                        }
                        break;
                    case 6 : // swap target with qubit 2 and interchange control and qubit 2
                        if( !templ )
                        {
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                        }
                        else
                        {
                           // append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, control2, control );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            
                           // append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, control );
                        }
                        break;
                }
            }
            else
                assert( false );
        }
        else if ( is_v( gate ) )
        {
            const auto& tag = boost::any_cast<v_tag>( gate.type() );

            if( gate.controls().empty() ) // a V gate without controls
            {
                append_v( circ_IBM,  gate.controls(), target, tag.adjoint );
            }
            else // V gate with controls
            {
                //std::cout << "CNOT case " << map_method[control][target] << "\n";
                
                unsigned method;
                bool tag1;
                if (map_method[control][target] < map_method[target][control])
                {
                    method = map_method[control][target];
                    control = gate.controls().front().line();
                    target = gate.targets().front();
                    new_controls.clear();
                    new_controls.push_back( target );
                    old_controls.clear();
                    old_controls.push_back( control );
                }
                else
                {
                    method = map_method[target][control];
                    target  = gate.controls().front().line();
                    control = gate.targets().front();
                    new_controls.clear();
                    new_controls.push_back( target );
                    old_controls.clear();
                    old_controls.push_back( control );
                }
                // std::cout << method << std::endl;
                if( tag.adjoint )
                    tag1 = false; 
                else 
                    tag1 = true;
                switch ( method )
                {
                    case 1:
                        append_hadamard( circ_IBM, gate.targets().front() );
                        append_toffoli( circ_IBM, old_controls, target );
                        append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                        append_toffoli( circ_IBM, old_controls, target );
                        append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                        append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                        append_hadamard( circ_IBM, gate.targets().front() );
                        break;
                        
                    case 2 : // invert CNOT
                        append_hadamard( circ_IBM, gate.targets().front() );
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, target );
                        append_toffoli( circ_IBM, new_controls, control );
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, target );
                        append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, target );
                        append_toffoli( circ_IBM, new_controls, control );
                        append_hadamard( circ_IBM, control );
                        append_hadamard( circ_IBM, target );
                        append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                        append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                        append_hadamard( circ_IBM, gate.targets().front() );
                        break;
                    case 3 : // swap target with 2
                        if( !templ )
                        {
                            append_hadamard( circ_IBM, gate.targets().front() );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            
                            append_toffoli( circ_IBM, gate.controls(), 2u );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            
                            append_toffoli( circ_IBM, gate.controls(), 2u );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                            append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                            append_hadamard( circ_IBM, gate.targets().front() );
                        }
                        else // use the "template" transformation
                        {
                            append_hadamard( circ_IBM, gate.targets().front() );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, 2u );
                            append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, 2u );
                            append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                            append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                            append_hadamard( circ_IBM, gate.targets().front() );
                            
                        }
                        break;
                    case 4 : // map control to qubit 2, given CNOT(2,c) -- CBA(c,2)
                        if( !templ )
                        {
                            append_hadamard( circ_IBM, gate.targets().front() );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            
                            append_toffoli( circ_IBM, control2, target );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, control );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            
                            append_toffoli( circ_IBM, control2, target );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, control );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                            append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                            append_hadamard( circ_IBM, gate.targets().front() );
                        }
                        else
                        {
                            append_hadamard( circ_IBM, gate.targets().front() );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, target );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, target );
                            append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, target );
                            append_hadamard( circ_IBM, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, target );
                            append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                            append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                            append_hadamard( circ_IBM, gate.targets().front() );
                        }
                        break;
                    case 5 : // map control to qubit 2, given CNOT(c,2) -- CAB(c,2)
                        if( !templ )
                        {
                            append_hadamard( circ_IBM, gate.targets().front() );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            
                            append_toffoli( circ_IBM, control2, target );
                            
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            
                            append_toffoli( circ_IBM, control2, target );
                            
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, control );
                            append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                            append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                            append_hadamard( circ_IBM, gate.targets().front() );
                        }
                        else
                        {
                            append_hadamard( circ_IBM, gate.targets().front() );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_toffoli( circ_IBM, control2, target );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_toffoli( circ_IBM, control2, target );
                            append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_toffoli( circ_IBM, control2, target );
                            append_toffoli( circ_IBM, old_controls, 2u );
                            append_toffoli( circ_IBM, control2, target );
                            append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                            append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                            append_hadamard( circ_IBM, gate.targets().front() );
                        }
                        break;
                    case 6: // swap target with qubit 2 and interchange control and qubit 2
                        if( !templ )
                        {
                            append_hadamard( circ_IBM, gate.targets().front() );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            
                            append_hadamard( circ_IBM, control );
                            append_toffoli( circ_IBM, control2, control );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                            append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                            append_hadamard( circ_IBM, gate.targets().front() );
                        }
                        else
                        {
                            append_hadamard( circ_IBM, gate.targets().front() );
                            // append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, control2, control );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            
                           // append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, control );
                            append_pauli( circ_IBM,  target, pauli_axis::Z, 4u, tag1 );
                            // append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, control );
                            
                            append_toffoli( circ_IBM, control2, control );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            append_toffoli( circ_IBM, control2, control );
                            append_toffoli( circ_IBM, new_controls, 2u );
                            
                           // append_hadamard( circ_IBM, 2u );
                            append_hadamard( circ_IBM, target );
                            append_hadamard( circ_IBM, control );
                            append_pauli( circ_IBM,  gate.controls()[0].line(), pauli_axis::Z, 4u, !tag1 );
                            append_pauli( circ_IBM,  gate.targets().front(), pauli_axis::Z, 4u, !tag1 );
                            append_hadamard( circ_IBM, gate.targets().front() );
                        }
                        break;
                }
            }
        }
        else if ( is_pauli( gate ) )
        {
            const auto& tag = boost::any_cast<pauli_tag>( gate.type() );
            append_pauli( circ_IBM, target, tag.axis, tag.root, tag.adjoint );
            //            std::cout << "Pauli " << target <<  std::endl;
            
        }
        else if ( is_hadamard( gate ) )
        {
            append_hadamard( circ_IBM, target );
        }
        else if ( is_rotation( gate ) )
        {
            const auto& tag = boost::any_cast<rotation_tag>( gate.type() );
            append_rotation( circ_IBM, target, rotation_axis::Z, tag.rotation );
        }
        else
        {
            assert( false );
        }
        //        std::cout << "size = " << circ_IBM.num_gates() << std::endl;
    }
    
    
    return circ_IBM;
}


int levels(const circuit& circ, circuit& result )
{
    std::vector<int> glevel; /* level of the gates */
    
    result = circ;
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
        i++;
    }
    
    return max_lev;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
