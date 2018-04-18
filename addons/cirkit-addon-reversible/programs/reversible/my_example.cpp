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
 * @author Gerhard Dueck
 */

#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <boost/lexical_cast.hpp>

#include <reversible/circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/io/read_qc.hpp>
#include <reversible/optimization/simplify.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/functions/ibm_helper.hpp>

using namespace cirkit;

void all_stats(const circuit& circ, std::string fname );

inline void cpy_perm(int v1[], int v2[]){
    for(int i = 0; i < 5; i++)
        v1[i] = v2[i];
}

int main( int argc, char ** argv )
{
    std::string fname, fname_qc;
    while(std::cin >> fname){
        fname_qc = fname + ".qc";
        std::cout << fname_qc << "\n";
        circuit circ = read_qc( fname_qc );
        all_stats( circ, fname );
        
       // std::cout << circ << std::endl;
    }
    return 0;
}

void all_stats(const circuit& circ, std::string fname ){
    circuit circ_working = circ, tempc;
    circuit circ_IBM_QX[4], // stored in this order
    circ_IBM_QX_best[4];    // [0] best QX2 circuit using swap transformations
                            // [1] best QX4 circuit using swap transformations
                            // [2] best QX2 circuit using template transformations
                            // [3] best QX4 circuit using template transformations
    bool template_flag[4] = {false, false, true, true};
    bool QX2_flag[4] = {true, false, true, false}, first_time = true;
    
    int perm[5] = {0, 1, 2, 3, 4}, inv_perm[5],
    best_perm[4][5] = {{0, 1, 2, 3, 4}, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4}};
    int best[4] = {INT_MAX, INT_MAX, INT_MAX, INT_MAX};
    
    std::ofstream outfile;
    outfile.open("stats.tex", std::ios_base::app);
    outfile << fname << " & " << circ_working.lines();
    
    unsigned start = circ_working.lines()+1;
    for(unsigned i = start ; i <= 5u; i++)
    {
        add_line_to_circuit( circ_working, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
    }

    do
    {
        permute_lines( circ_working , perm );
        for(int i = 0; i < 4; i++)
        {
            if( QX2_flag[i] )
            {
                circ_IBM_QX[i] = transform_to_IBMQ( circ_working, map_method_qx2, template_flag[i] );
            }
            else
            {
                circ_IBM_QX[i] = transform_to_IBMQ( circ_working, map_method_qx4, template_flag[i] );
            }
            circ_IBM_QX[i] = remove_dup_gates( circ_IBM_QX[i] );
            if( best[i] > circ_IBM_QX[i].num_gates() )
            {
                best[i] = circ_IBM_QX[i].num_gates();
                circ_IBM_QX_best[i] = circ_IBM_QX[i];
                cpy_perm(best_perm[i], perm);
            }
        }
 
        // undo the permutation
        for( int i = 0; i < 5; i++ )
        {
            inv_perm[perm[i]] = i;
        }
        permute_lines( circ_working , inv_perm );
        if( first_time )
        {
            for( int i = 0; i < 4; i++ )
            {
                outfile << " & " << best[i] << " & " << levels(circ_IBM_QX_best[i], tempc);
                std::cout << "level = " << levels(circ_IBM_QX_best[i], tempc) << " gates = " << best[i] << std::endl;
            }
            outfile << " \\\\ \\hline\n & ";
            first_time = false;
        }
    } while ( std::next_permutation(perm,perm+5) );
    
    
    for( int i = 0; i < 4; i++ )
    {
        outfile << " & " << best[i] << " & " << levels(circ_IBM_QX_best[i], tempc);
        std::cout << "level = " << levels(circ_IBM_QX_best[i], tempc) << " gates = " << best[i] << std::endl;
    }
    outfile << " \\\\ \\hline\n";
    outfile.close();
    
}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
