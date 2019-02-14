/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
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

#include "qxg.hpp"

#include <fstream>
#include <algorithm>
#include <core/utils/timer.hpp>
#include <alice/rules.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/program_options.hpp>
#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/io/write_qc.hpp>
#include <reversible/io/print_circuit.hpp>
#include <cli/commands/ibm.hpp>
#include <cli/commands/permute_lines.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/ibm_helper.hpp>

using namespace boost::program_options;

using matrix = std::vector<std::vector<int>>;

//Matrix with the cost of each possible cnot (QX2)
matrix map_qx2 = { {0,0,0,10,10},
                                {4,0,0,10,10},
                                {4,4,0,4,4},
                                {10,10,0,0,0},
                                {10,10,0,4,0}};
//Matrix with the cost of each possible cnot (QX4)
matrix map_qx4 = { {0,4,4,10,10},
                                {0,0,4,10,10},
                                {0,0,0,4,0},
                                {12,12,0,0,0},
                                {10,10,4,4,0}};
//Matrix with the cost of each possible cnot (QX3)
matrix map_qx3 = { {0, 0, 10, 24, 38, 52, 74, 80, 94, 88, 66, 52, 46, 32, 10, 4},
                                {4, 0, 0, 10, 24, 38, 80, 94, 108, 94, 80, 66, 52, 38, 24, 18},
                                {18, 4, 0, 0, 10, 24, 66, 80, 94, 80, 66, 52, 38, 24, 10, 24},
                                {24, 18, 4, 0, 4, 10, 52, 66, 80, 66, 52, 38, 24, 10, 0, 10},
                                {38, 24, 10, 0, 0, 0, 38, 52, 66, 52, 38, 24, 10, 4, 10, 24},
                                {52, 46, 32, 10, 4, 0, 32, 38, 52, 46, 24, 10, 4, 10, 24, 46},
                                {66, 80, 66, 52, 38, 24, 0, 0, 10, 24, 10, 0, 10, 24, 38, 52},
                                {80, 94, 80, 66, 52, 38, 4, 0, 4, 10, 0, 10, 24, 38, 52, 66},
                                {94, 108, 94, 80, 66, 52, 10, 0, 0, 4, 10, 24, 38, 52, 66, 80},
                                {80, 94, 80, 66, 52, 38, 24, 10, 0, 0, 0, 10, 24, 38, 52, 66},
                                {66, 80, 74, 52, 38, 24, 18, 4, 10, 4, 0, 4, 18, 24, 38, 60},
                                {52, 66, 60, 38, 24, 10, 4, 10, 24, 10, 0, 0, 4, 10, 24, 46},
                                {38, 52, 38, 24, 10, 0, 10, 24, 38, 24, 10, 0, 0, 0, 10, 24},
                                {24, 38, 24, 10, 0, 10, 32, 38, 52, 46, 24, 10, 4, 0, 0, 10},
                                {10, 24, 18, 4, 10, 24, 46, 52, 66, 60, 38, 24, 18, 4, 0, 4},
                                {0, 10, 24, 10, 24, 38, 52, 66, 80, 66, 52, 38, 24, 10, 0, 0}};

matrix path_qx3={{0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1},
                            {-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,-1,0,-1,0,0,0,0,0,0,0,0,0,1,0},
                            {0,0,0,1,0,1,0,0,0,0,0,0,0,-1,0,0},
                            {0,0,0,0,-1,0,0,0,0,0,0,0,-1,0,0,0},
                            {0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0},
                            {0,0,0,0,0,0,-1,0,-1,0,1,0,0,0,0,0},
                            {0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0},
                            {0,0,0,0,0,0,0,-1,0,-1,0,-1,0,0,0,0},
                            {0,0,0,0,0,0,-1,0,0,0,1,0,-1,0,0,0},
                            {0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0},
                            {0,0,0,0,1,0,0,0,0,0,0,0,-1,0,1,0},
                            {0,0,0,-1,0,0,0,0,0,0,0,0,0,-1,0,-1},
                            {1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}};

matrix path_qx5={{0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1},
                            {1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,-1},
                            {0,0,-1,0,1,0,0,0,0,0,0,0,0,0,1,0},
                            {0,0,0,-1,0,-1,0,0,0,0,0,0,0,-1,0,0},
                            {0,0,0,0,1,0,-1,0,0,0,0,0,-1,0,0,0},
                            {0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0},
                            {0,0,0,0,0,0,-1,0,-1,0,1,0,0,0,0,0},
                            {0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0},
                            {0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0},
                            {0,0,0,0,0,0,0,-1,0,-1,0,-1,0,0,0,0},
                            {0,0,0,0,0,0,-1,0,0,0,1,0,-1,0,0,0},
                            {0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0},
                            {0,0,0,0,1,0,0,0,0,0,0,0,-1,0,1,0},
                            {0,0,0,-1,0,0,0,0,0,0,0,0,0,-1,0,-1},
                            {1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0}};

matrix map_qx5={{0,4,10,20,32,42,52,64,74,64,54,42,30,20,10,4},
                            {0,0,0,12,24,34,44,56,66,76,66,54,42,32,22,10},
                            {10,4,0,0,12,22,32,44,54,64,54,42,30,20,10,4},
                            {20,14,4,0,0,10,20,32,42,52,44,32,20,10,0,10},
                            {30,24,14,4,0,4,14,20,30,40,32,20,14,4,10,20},
                            {42,30,20,10,0,0,4,10,20,30,22,10,4,10,22,30},
                            {54,42,32,22,12,0,0,0,10,20,12,0,10,22,34,42},
                            {64,52,42,32,22,10,4,0,4,10,0,10,20,32,44,52},
                            {76,64,54,44,34,22,10,0,0,4,10,20,30,42,54,64},
                            {66,74,64,54,44,32,20,10,0,0,0,10,20,32,44,54},
                            {54,62,52,42,32,20,14,4,10,4,0,4,14,20,32,42},
                            {44,52,42,32,22,10,4,10,20,10,0,0,4,10,22,32},
                            {34,42,32,22,12,0,10,22,32,22,12,0,0,0,12,22},
                            {22,30,20,10,0,10,20,32,42,32,22,10,4,0,0,10},
                            {10,20,10,4,10,20,30,42,52,42,32,20,14,4,0,4},
                            {0,10,0,10,22,32,42,54,64,54,44,32,20,10,0,0}};

matrix map_qx20= {{0, 0, 10, 66, 52, 0, 10, 10, 38, 52, 10, 10, 24, 24, 38, 24, 24, 24, 38, 38},
                                {0, 0, 0, 52, 38, 10, 0, 0, 24, 38, 10, 10, 10, 10, 24, 24, 24, 24, 24, 24},
                                {10, 0, 0, 52, 38, 24, 10, 0, 24, 38, 24, 24, 10, 10, 24, 38, 24, 24, 24, 24},
                                {66, 52, 52, 0, 0, 52, 52, 38, 10, 0, 52, 38, 24, 24, 10, 52, 38, 38, 24, 24},
                                {52, 38, 38, 0, 0, 38, 38, 24, 0, 0, 38, 24, 10, 10, 10, 38, 24, 24, 24, 24},
                                {0, 10, 24, 52, 38, 0, 0, 10, 24, 38, 0, 0, 10, 24, 38, 10, 10, 10, 24, 38},
                                {10, 0, 10, 52, 38, 0, 0, 0, 24, 38, 0, 0, 10, 10, 24, 10, 10, 10, 24, 24},
                                {10, 0, 0, 38, 24, 10, 0, 0, 10, 24, 10, 10, 0, 0, 10, 24, 10, 10, 10, 10},
                                {38, 24, 24, 10, 0, 24, 24, 10, 0, 0, 24, 10, 0, 0, 10, 24, 10, 10, 10, 10},
                                {52, 38, 38, 0, 0, 38, 38, 24, 0, 0, 38, 24, 10, 10, 0, 38, 24, 24, 10, 10},
                                {10, 10, 24, 52, 38, 0, 0, 10, 24, 38, 0, 0, 10, 24, 38, 0, 10, 10, 24, 38},
                                {10, 10, 24, 38, 24, 0, 0, 10, 10, 24, 0, 0, 0, 10, 24, 10, 0, 0, 10, 24},
                                {24, 10, 10, 24, 10, 10, 10, 0, 0, 10, 10, 0, 0, 0, 10, 10, 0, 0, 10, 10},
                                {24, 10, 10, 24, 10, 24, 10, 0, 0, 10, 24, 10, 0, 0, 0, 24, 10, 10, 0, 0},
                                {38, 24, 24, 10, 10, 38, 24, 10, 10, 0, 38, 24, 10, 0, 0, 38, 24, 10, 0, 0},
                                {24, 24, 38, 52, 38, 10, 10, 24, 24, 38, 0, 10, 10, 24, 38, 0, 0, 10, 24, 38},
                                {24, 24, 24, 38, 24, 10, 10, 10, 10, 24, 10, 0, 0, 10, 24, 0, 0, 0, 10, 24},
                                {24, 24, 24, 38, 24, 10, 10, 10, 10, 24, 10, 0, 0, 10, 10, 10, 0, 0, 0, 10},
                                {38, 24, 24, 24, 24, 24, 24, 10, 10, 10, 24, 10, 10, 0, 0, 24, 10, 0, 0, 0},
                                {38, 24, 24, 24, 24, 38, 24, 10, 10, 10, 38, 24, 10, 0, 0, 38, 24, 10, 0, 0}};

matrix path_qx20 = {{0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    {1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
                                    {0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
                                    {0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
                                    {0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
                                    {1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0},
                                    {0,1,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0},
                                    {0,1,1,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0},
                                    {0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0},
                                    {0,0,0,1,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0},
                                    {0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0},
                                    {0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,1,1,0,0},
                                    {0,0,0,0,0,0,0,1,1,0,0,1,0,1,0,0,1,1,0,0},
                                    {0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,1,1},
                                    {0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1},
                                    {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0},
                                    {0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,1,0,0},
                                    {0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0},
                                    {0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,1},
                                    {0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0}};

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

qxg_command::qxg_command( const environment::ptr& env )
    : cirkit_command( env, "IBM QX mapping algorithm" )
{
    opts.add_options()
     ( "qx4,4", "IBM QX4 matrix")
     ( "qx3,3", "IBM QX3 matrix")
     ( "qx5,5", "IBM QX5 matrix")
     ( "QS1_1,20", "QS1_1 matrix")
    ;
  add_new_option();
}


command::rules_t qxg_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

void get_cnots_matrix(circuit& circ, matrix& cnots)
{
    for ( const auto& gate : circ )
        if ( is_toffoli( gate ) && !gate.controls().empty() )
            ++cnots[gate.controls().front().line()][gate.targets().front()];
}

void print_matrix(matrix& m)
{
    std::cout << std::endl;
    for (int i = 0; i < m.size(); ++i){
        for (int j = 0; j < m.size(); ++j){
            std::cout << " " << m[i][j];   
        }
        std::cout << std::endl;
    }
}

std::pair<int,int> get_most_used_cnot(matrix& m)
{
    int h = 0, x, y;
    for (int i = 0; i < m.size(); ++i)
        for (int j = 0; j < m.size(); ++j)
            if( m[i][j] > h ) { h = m[i][j]; x = i; y = j; }
    return std::make_pair(x,y);
}

std::pair<int,int> get_cheaper_cnot(matrix& m)
{
    int h = 0, x, y;
    for (int i = 0; i < m.size(); ++i)
        for (int j = 0; j < m.size(); ++j)
            if( m[i][j] > h ) { h = m[i][j]; x = i; y = j; }
    return std::make_pair(x,y);
}


void create_matrix(matrix& m, unsigned size)
{
    std::vector <int> p;
    for (int i = 0; i < size; ++i)
        p.push_back(0);
    for (int i = 0; i < size; ++i)
        m.push_back(p);
}

void get_best_qubit( matrix& m)
{
    unsigned control = 0, target = 0, qc, qt;
    for (int i = 0; i < m.size(); ++i)
    {
        unsigned c = 0, t= 0;
        for (int j = 0; j < m.size(); ++j)
        {
            if (m[i][j] == -1)
                ++t;
            if (m[i][j] == 1)
                ++c;
        }
        if( c > control)
        {
            qc = i;
            control = c;
        }
        if( t > target)
        {
            qt = i;
            target = t;
        }
    }
    std::cout << "control: " << qc << " - " << control << std::endl;
    std::cout << "target: " << qt << " - " << target << std::endl;
}


void qxg(circuit& circ, matrix& map, matrix& path )
{
    std::pair<int,int> cnot;
    std::map<int, int> permutation;
    matrix cnots;
    
    create_matrix(cnots, map.size());
    get_cnots_matrix(circ, cnots);
    // print_matrix(cnots);
    cnot = get_most_used_cnot(cnots);
    std::cout << "most used: CNOT(" << cnot.first << "," << cnot.second << ")" << std::endl;
    get_best_qubit(path);
}

bool qxg_command::execute()
{
    auto& circuits = env->store<circuit>();
    circuit circ = circuits.current();
    //auto settings = make_settings();
    
    if ( is_set( "QS1_1" ) )
    {
        if(circ.lines() > 20)
        {
            std::cout << "Only up to 20 variables!" << std::endl;
            return true;
        }
        qxg(circ, map_qx20, path_qx20);
    }
    else if ( is_set( "qx3" ) )
    {
        if(circ.lines() > 16)
        {
            std::cout << "Only up to 16 variables!" << std::endl;
            return true;
        }
        qxg(circ, map_qx3, map_qx3);
    }
    else if ( is_set( "qx5" ) )
    {
        if(circ.lines() > 16)
        {
            std::cout << "Only up to 16 variables!" << std::endl;
            return true;
        }
        qxg(circ, map_qx5, path_qx5);
    }
    else
    {
        if(circ.lines() > 5)
        {
            std::cout << "Only up to 5 variables! Try another option." << std::endl;
            return true;
        }
        if ( is_set( "qx4" ) )
            qxg(circ, map_qx4, map_qx4);
        else
            qxg(circ, map_qx2, map_qx2);
    }

    return true;
}

command::log_opt_t qxg_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}


}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
