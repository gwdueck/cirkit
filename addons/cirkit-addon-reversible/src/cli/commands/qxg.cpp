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

using namespace boost::program_options;

typedef std::vector<std::vector<int>> matrix;

//Matrix with the cost of each possible cnot (QX2)
static const matrix map_qx2 = { {0,0,0,10,10},{4,0,0,10,10},{4,4,0,4,4},{10,10,0,0,0},{10,10,0,4,0}};
//Matrix with the cost of each possible cnot (QX4)
static const matrix map_qx4 = { {0,4,4,18,10},{0,0,4,18,10},{0,0,0,4,0},{10,10,0,0,0},{10,10,4,4,0}};
//Matrix with the cost of each possible cnot (QX3)
static const matrix map_qx3 = { {0, 0, 10, 24, 38, 52, 74, 80, 94, 88, 66, 52, 46, 32, 10, 4},
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

static const matrix path_qx3={{0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1},
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

static const matrix map_qx20= {{0, 0, 10, 66, 52, 0, 10, 10, 38, 52, 10, 10, 24, 24, 38, 24, 24, 24, 38, 38},
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

static const matrix path_qx20 = {{0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
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
     ( "QS1_1,20", "QS1_1 matrix")
    ;
  add_new_option();
}


command::rules_t qxg_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

//Change two lines of the circuit
void manipulate_matrix( matrix& m1, int x, int y, matrix& m2, std::vector<int>& perm, const matrix& map )
{
    int aux;

    aux = perm[x];
    perm[x] = perm[y];
    perm[y] = aux;
    for(int j=0; j<m1.size(); ++j)
    {
        aux = m1[x][j];
        m1[x][j] = m1[y][j];
        m1[y][j] = aux;
        m2[x][j] = m1[x][j] * map[x][j]; 
        m2[y][j] = m1[y][j] * map[y][j];
    }
    for(int i=0; i<m1.size(); ++i)
    {
        aux = m1[i][x];
        m1[i][x] = m1[i][y];
        m1[i][y] = aux;
        m2[i][x] = m1[i][x] * map[i][x]; 
        m2[i][y] = m1[i][y] * map[i][y];
    }
}

void print_matrix( const matrix& m)
{
    for(int i=0; i<m.size(); ++i)
    {
        for(int j=0; j<m.size(); ++j)
        {
            std::cout << m[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// //Search matrix for the qubit with higher cost
// int higher_cost( const matrix& m1, const matrix& m2, const std::vector<int>& p)
// {
//     int cost, qtd_cnot; 
//     int higher_cost = 0, index = 0, higher_qtd_cnot = 0;
    
//     for(int i=0; i<m1.size(); ++i)
//     {
//         cost = 0;
//         qtd_cnot = 0;
//         if (std::find(p.begin(), p.end(), i) == p.end())
//         {
//             for(int j=0; j<m1.size(); ++j)
//             {
//                 cost += m1[i][j] + m1[j][i];
//                 qtd_cnot += m2[i][j] + m2[j][i];
//             }
//             if( cost > higher_cost)
//             {
//                 higher_cost = cost;
//                 higher_qtd_cnot = qtd_cnot;
//                 index = i;
//             }
//             else if(cost == higher_cost && qtd_cnot > higher_qtd_cnot)
//             {
//                 higher_qtd_cnot = qtd_cnot;
//                 index = i;
//             }
//         }
//     }
//     return index;
// }

// //Search matrix for the pair qubit with higher cost
int higher_cost( const matrix& m1, const matrix& m2, const std::vector<int>& p)
{
    int cost, qtd_cnot; 
    int higher_cost = 0, index = 0, higher_qtd_cnot = 0;
    
    for(int i=0; i<m1.size(); ++i)
    {
        cost = 0;
        qtd_cnot = 0;
        if (std::find(p.begin(), p.end(), i) == p.end())
        {
            for(int j=0; j<m1.size(); ++j)
            {
                cost = m1[i][j]; //+ m1[j][i];
                qtd_cnot = m2[i][j]; //+ m2[j][i];
                if( cost > higher_cost)
                {
                    higher_cost = cost;
                    index = i;
                }
                else if(cost == higher_cost && qtd_cnot > higher_qtd_cnot)
                {
                    higher_qtd_cnot = qtd_cnot;
                    index = i;
                }
            }
        }
    }
    return index;
}

//Create cnots matrix and cost matrix
unsigned int initial_matrix(circuit circ, matrix& cnots, matrix& map_cost, const matrix& map)
{
    unsigned int target, control;
    unsigned int cost = 0;
    for ( const auto& gate : circ )
    {
        if ( is_toffoli( gate ) && !gate.controls().empty() )
        {   
            target = gate.targets().front();
            control = gate.controls().front().line();
            cost += map[control][target];
            ++cnots[control][target];
            map_cost[control][target] = cnots[control][target] * map[control][target];
        }
    }
    return cost;
}

//Update the circuit cost
unsigned int matrix_cost(const matrix& m)
{
    unsigned int cost = 0;
    for(int i=0; i<m.size(); ++i)
        for(int j=0; j<m.size(); ++j)
            cost += m[i][j];
    return cost;
}

void print_permutation( const std::vector<int>& perm )
{
    for(int i=0; i<perm.size(); ++i)
        std::cout << " " << perm[i];    
    std::cout << std::endl;
}

void print_results( const matrix& cnots, const std::vector<int>& perm, const unsigned int cost)
{
    std::cout << "Best permutation found (gates =  " << cost << "):";
    //print_permutation(perm);
    
    //This print the permutation like the ibm command
    //std::cout << "ibm representation: ";
    for(int i=0; i<cnots.size(); ++i)
    {
        for(int j=0; j<cnots.size(); ++j)
        {
            if(perm[j] == i)
            {
                std::cout << " " << j;
                break;
            }
        }
    }
    std::cout << std::endl;  
}

int search_qubit_column(const matrix& mapping, const unsigned int target, const unsigned int value)
{
    for (int i = 0; i < mapping.size(); ++i)
    {
        if(mapping[i][target] == value && i != target)
            return i;
    }
    return (-1);
}

int search_qubit_row(const matrix& mapping, const unsigned int control, const unsigned int value)
{
    for (int i = 0; i < mapping.size(); ++i)
    {
        if(mapping[control][i] == value && i != control)
            return i;   
    }
    return (-1);
}

void find_path(const matrix& mapping, const unsigned int control, const unsigned int target, std::vector<int>& permute, const matrix& path, const unsigned int& path_size)
{
    unsigned int x;
    if(permute.size() >= path_size)
    {
        // std::cout << "Tamanho maximo" << std::endl;
        return;
    }
    else
    {
        for (int i = 0; i < path.size(); ++i)
        {
            if (path[control][i] != 0 && std::find(permute.begin(), permute.end(), i) == permute.end())
            {   
                permute.push_back(i);
                // for (int c = 0; c < permute.size(); ++c){std::cout << " " << permute[c];}std::cout << std::endl;
                // std::cout << "target: " << target << std::endl;
                if(i == target)
                {
                    // std::cout << "AAAAAAAAAAAA" << std::endl;
                    x = (permute.size() - 2) * 14;
                    if(path[control][i] == -1)
                        if(x + 4 == mapping[permute[0]][target])
                        {
                            // std::cout << "Achou1!" << std::endl;
                            return;
                        }
                    if(path[permute[0]][permute[1]] == -1)
                    {
                        if(x + 4 == mapping[permute[0]][target])
                        {
                            for (int y = 0; y < permute.size(); ++y)
                                permute[y] = permute[y] * -1;
                            // std::cout << "Achou2!" << std::endl;
                            return;  
                        }
                    }
                    if(path[control][i] == 1)
                        if(x - 4 == mapping[permute[0]][target])
                        {
                            // std::cout << "Achou3!" << std::endl;
                            return;
                        }
                    if(path[permute[0]][permute[1]] == 1)
                    {
                        if(x - 4 == mapping[permute[0]][target])
                        {
                            for (int y = 0; y < permute.size(); ++y)
                                permute[y] = permute[y] * -1;
                            // std::cout << "Achou4!" << std::endl;
                            return;
                        }
                    }
                    permute.pop_back();
                }
                else
                {
                    find_path(mapping, i, target, permute, path, path_size);
                    if(permute.back() != target && permute.back() != target*-1)
                    {
                        permute.pop_back();
                    }
                    else
                    {
                        // std::cout << "Achou5!" << std::endl;
                    }
                }
            }   
        }
    }
}

void invert_cnot(circuit& circ, const unsigned int control, const unsigned int target)
{
    std::vector<unsigned int> controls;
    controls.push_back( target );
    append_hadamard( circ, control );
    append_hadamard( circ, target );
    append_toffoli( circ, controls, control );
    append_hadamard( circ, control );
    append_hadamard( circ, target );
    return circ;
}

void swap_gates_seven(circuit& circ, const unsigned int control, const unsigned int target)
{
    std::vector<unsigned int> controls;
    controls.push_back( control );
    append_toffoli( circ, controls, target );
    append_hadamard( circ, control );
    append_hadamard( circ, target );
    append_toffoli( circ, controls, target );
    append_hadamard( circ, control );
    append_hadamard( circ, target );
    append_toffoli( circ, controls, target );
    return circ;
}

void swap_gates_five(circuit& circ, const unsigned int control, const unsigned int target, const unsigned hadamard)
{
    std::vector<unsigned int> controls;
    controls.push_back( control );
    append_toffoli( circ, controls, target );
    append_hadamard( circ, control );
    append_hadamard( circ, target );
    append_toffoli( circ, controls, target );
    append_hadamard( circ, hadamard );
    return circ;
}

void swap_gates_five_back(circuit& circ, const unsigned int control, const unsigned int target, const unsigned hadamard)
{
    std::vector<unsigned int> controls;
    controls.push_back( control );
    append_hadamard( circ, hadamard );
    append_toffoli( circ, controls, target );
    append_hadamard( circ, control );
    append_hadamard( circ, target );
    append_toffoli( circ, controls, target );
}

circuit matrix_to_circuit( circuit circ, const matrix& cnots, const std::vector<int>& perm, const matrix& mapping, const matrix& path)
{
    //piece of code from ibm.cpp
    std::vector<int> permute;
    unsigned int start = circ.lines() + 1;
    unsigned int target, control;
    unsigned int qubit;
    unsigned int comeco, fim;

    std::vector<unsigned int> new_controls;
    circuit circ_qx;
   
    for(unsigned int i = start ; i <= cnots.size(); i++)
    {
        add_line_to_circuit( circ, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
    }   
    
    copy_metadata(circ, circ_qx);
    
    // std::cout << "tamanho do circuito no comeco: " << circ_qx.num_gates() << std::endl;
    for(int i=0; i<cnots.size(); ++i)
        for(int j=0; j<cnots.size(); ++j)
            if(perm[j] == i)
                permute.push_back(j);
    
    permute_lines( circ , &permute[0] );
    permute.clear();
    // iterate through the gates
    for ( const auto& gate : circ )
    {
        target = gate.targets().front();
        new_controls.clear();
        new_controls.push_back( target );
        if( !gate.controls().empty() )
        {
            control = gate.controls().front().line();
        }
        
        if ( is_toffoli( gate ) )
        {
            if( gate.controls().empty() ) // a NOT gate
            {
                append_toffoli( circ_qx, gate.controls(), target );
            }
            else // CNOT gate
            {
                // std::cout << "controle: " << control << " target: " << target << std::endl;
                // std::cout << "esperado: " << mapping[control][target] << std::endl;
                // comeco = circ_qx.num_gates();
                // std::cout << "comeco: " << comeco << std::endl;
                if ( mapping[control][target] == 0 )
                {
                    append_toffoli( circ_qx, gate.controls(), target );
                }
                else if( mapping[control][target] == 4 ) // invert CNOT
                {
                    invert_cnot(circ_qx, control, target);
                }
                else
                {
                    unsigned int path_size;
                    permute.clear();
                    if((mapping[control][target] - 4) % 7 == 0)
                    {
                        // std::cout << "-4" << std::endl;
                        path_size = ((mapping[control][target] - 4) / 14) + 2;
                        permute.clear();
                        permute.push_back(control);
                        find_path(mapping, control, target, permute, path, path_size);
                        // for (int i = 0; i < permute.size(); ++i){std::cout << " " << permute[i];}std::cout << std::endl;
                        if(permute[0] < 0 || permute[1] < 0)
                        {
                            // std::cout << "-4 negativo" << std::endl;
                            for (int y = 0; y < permute.size(); ++y)
                                permute[y] = permute[y] * -1;
                            for (int i = permute.size() - 1; i > 1; --i)
                            {
                                if(path[permute[i]][permute[i-1]] == -1)
                                    swap_gates_seven(circ_qx, permute[i-1], permute[i]);
                                else
                                    swap_gates_seven(circ_qx, permute[i], permute[i-1]);
                            }                

                            invert_cnot(circ_qx, permute[0], permute[1]);                           
                            
                            for (int i = 1; i < permute.size() - 1; ++i)
                            {
                                if(path[permute[i]][permute[i+1]] == -1)
                                    swap_gates_seven(circ_qx, permute[i+1], permute[i]);
                                else
                                    swap_gates_seven(circ_qx, permute[i], permute[i+1]);
                            } 
                        }
                        else
                        {
                            // std::cout << "-4 positivo" << std::endl;
                            for (int i = 0; i < permute.size() - 1; ++i)
                            {
                                if(i == permute.size() - 2)
                                {
                                    invert_cnot(circ_qx, permute[i], permute[i+1]); 
                                }
                                else
                                {
                                    if(path[permute[i]][permute[i+1]] == -1)
                                        swap_gates_seven(circ_qx, permute[i+1], permute[i]);
                                    else
                                        swap_gates_seven(circ_qx, permute[i], permute[i+1]);
                                }
                                
                            }                 
                                         
                            for (int i = permute.size() - 1; i > 1; --i)
                            {
                                if(path[permute[i]][permute[i-1]] == -1)
                                    swap_gates_seven(circ_qx, permute[i-1], permute[i]);
                                else
                                    swap_gates_seven(circ_qx, permute[i], permute[i-1]);
                            }
                        }
                    }
                    else if((mapping[control][target] + 4) % 7 == 0)
                    {
                        // std::cout << "+4" << std::endl;
                        path_size = ((mapping[control][target] + 4) / 14) + 2;
                        permute.clear();
                        permute.push_back(control);
                        find_path(mapping, control, target, permute, path, path_size);
                        // for (int i = 0; i < permute.size(); ++i){std::cout << " " << permute[i];}std::cout << std::endl;
                        if(permute[0] < 0 || permute[1] < 0)
                        {
                            // std::cout << "+4 negativo" << std::endl;
                            for (int y = 0; y < permute.size(); ++y)
                                permute[y] = permute[y] * -1;
                            for (int i = permute.size() - 1; i > 2; --i)
                            {
                                if(path[permute[i]][permute[i-1]] == -1)
                                    swap_gates_seven(circ_qx, permute[i-1], permute[i]);
                                else
                                    swap_gates_seven(circ_qx, permute[i], permute[i-1]);
                            }                
                            
                            if(path[permute[2]][permute[1]] == -1)
                                swap_gates_five(circ_qx, permute[1], permute[2], permute[2]);
                            else
                                swap_gates_five(circ_qx, permute[2], permute[1], permute[1]);
                            
                            new_controls.clear();
                            new_controls.push_back(permute[0]);
                            append_toffoli(circ_qx, new_controls, permute[1]);

                            if(path[permute[2]][permute[1]] == -1)
                                swap_gates_five_back(circ_qx, permute[1], permute[2], permute[2]);
                            else
                                swap_gates_five_back(circ_qx, permute[2], permute[1], permute[1]);

                            for (int i = 2; i < permute.size() - 1; ++i)
                            {
                                if(path[permute[i]][permute[i+1]] == -1)
                                    swap_gates_seven(circ_qx, permute[i+1], permute[i]);
                                else
                                    swap_gates_seven(circ_qx, permute[i], permute[i+1]);
                            } 
                        }
                        else
                        {
                            // std::cout << "+4 positivo" << std::endl;
                            for (int i = 0; i <= permute.size() - 3; ++i)
                            {
                                // std::cout << "primeira fase: " << permute[i] << std::endl;
                                if(i == permute.size() - 3)
                                {
                                    // std::cout << "igual" << std::endl;
                                    if(path[permute[i]][permute[i+1]] == -1)
                                        swap_gates_five(circ_qx, permute[i+1], permute[i], permute[i+1]);
                                    else
                                        swap_gates_five(circ_qx, permute[i], permute[i+1], permute[i]);

                                    new_controls.clear();
                                    new_controls.push_back(permute[i+1]);
                                    append_toffoli(circ_qx, new_controls, permute[i+2]);

                                    if(path[permute[i+1]][permute[i]] == -1)
                                        swap_gates_five_back(circ_qx, permute[i], permute[i+1], permute[i]);
                                    else
                                        swap_gates_five_back(circ_qx, permute[i+1], permute[i], permute[i+1]);
                                }
                                else
                                {
                                    // std::cout << "nao igual" << std::endl;
                                    if(path[permute[i]][permute[i+1]] == -1)
                                        swap_gates_seven(circ_qx, permute[i+1], permute[i]);
                                    else
                                        swap_gates_seven(circ_qx, permute[i], permute[i+1]);
                                }
                            }

                            for (int i = permute.size() - 2; i > 1; --i)
                            {
                                // std::cout << "segunda fase: " << permute[i] << std::endl;
                                if(path[permute[i]][permute[i-1]] == -1)
                                    swap_gates_seven(circ_qx, permute[i-1], permute[i]);
                                else
                                    swap_gates_seven(circ_qx, permute[i], permute[i-1]);
                            } 

                        }
                    }                 
                }
            }
            // fim = circ_qx.num_gates();
            // std::cout << "fim: " << fim << std::endl;
            // std::cout << "diferenca: " << fim - comeco << std::endl;

            // for (int i = 0; i < permute.size(); ++i){std::cout << " " << permute[i];}std::cout << std::endl;
        }
        else if ( is_pauli( gate ) )
        {
            const auto& tag = boost::any_cast<pauli_tag>( gate.type() );
            append_pauli( circ_qx, target, tag.axis, tag.root, tag.adjoint );
            
        }
        else if ( is_hadamard( gate ) )
        {
            append_hadamard( circ_qx, target );
        }
        else
        {
            assert( false );
        }
    }
    //std::cout << "TAMANHO: " << circ_qx.num_gates() << std::endl;
    return circ_qx;
}

circuit qxg(circuit& circ, const matrix& map, const matrix& path, properties::ptr& statistics )
{
    properties_timer t( statistics );
    circuit circ_qx;
    unsigned int cost, lower_cost, h, c, q;
    std::vector <int> p;
    std::vector<int> perm;
    std::vector<int> best_perm;
    matrix cnots;
    matrix map_cost;
    matrix aux;

    for (int i = 0; i < map.size(); ++i)
        p.push_back(0);
    
    for (int i = 0; i < map.size(); ++i)
    {
        cnots.push_back(p);
        map_cost.push_back(p);
        perm.push_back(i);
        best_perm.push_back(i);
        aux.push_back(p);
    }
    p.clear();
    cost = initial_matrix(circ, cnots, map_cost, map);
    cost = cost + circ.num_gates();
    //std::cout << circ.num_gates() << std::endl;
    // std::cout << "initial cost: " << cost << std::endl;
    std::cout << "initial gates: " << circ.num_gates() << std::endl;
    lower_cost = cost;
    // std::cout << "cnots matrix: " << std::endl;
    // print_matrix(cnots);
    // std::cout << "cost matrix: " << std::endl;
    // print_matrix(map_cost);
    
    // srand (time(NULL));
    unsigned int it = 0;
    do
    {
        // h = rand() % cnots.size();
        // i = rand() % cnots.size();
        h = higher_cost(map_cost, cnots, p);
        //std::cout << "Higher cost: " << h << std::endl;
        q = h;
        for(unsigned int i=0; i<cnots.size(); ++i)
        {
            manipulate_matrix( cnots, h, i, map_cost, perm, map );
            cost = matrix_cost(map_cost) + circ.num_gates();
            if(cost < lower_cost)
            {
                q = i;
                lower_cost = cost;
                //std::cout << "changed [" << h << "] [" << i << "] CUSTO: " << lower_cost << std::endl;
                it = 0;
                for(int j=0; j<cnots.size(); ++j)
                    best_perm[j] = perm[j];
            }
            manipulate_matrix( cnots, i, h, map_cost, perm, map );
        }
        if(h == q)
        {
            //std::cout << "Not better permutation found! [" << h << "] in the list" << std::endl;
            if (std::find(p.begin(), p.end(), h) == p.end())
                p.push_back(h);
            else
                p.clear();
            // for (int i = 0; i < p.size(); ++i)
            // {
            //     std::cout << " " << p[i];
            // }
            // std::cout << std::endl;
        }
        //std::cout << "Swapping [" << h << "] [" << q << "]" << std::endl;  
        manipulate_matrix( cnots, h, q, map_cost, perm, map );
        it++;
        // if(it % cnots.size() == 0)
        // {
        //     std::cout << "limpando p..." << std::endl;
        //     p.clear();
        // }
    } while (it < 3 * cnots.size());
    circ_qx = matrix_to_circuit(circ, cnots, best_perm, map, path);
    //circ_qx = remove_dup_gates( circ_qx );
    //print_results(cnots, best_perm, circ_qx.num_gates());
    return circ_qx;
}

circuit optimize_circuit(circuit& circ,  properties::ptr& statistics)
{
    properties_timer t( statistics );
    return remove_dup_gates( circ );
}

bool qxg_command::execute()
{
    auto& circuits = env->store<circuit>();
    circuit aux = circuits.current();
    circuit circ_qx, circ;
    copy_circuit(aux, circ);
    //auto settings = make_settings();

    // unsigned int path_size;
    // if((map_qx3[12][0] - 4) % 7 == 0)
    //     path_size = ((map_qx3[12][0] - 4) / 14) + 2;
    // else if((map_qx3[12][0] + 4) % 7 == 0)
    //     path_size = ((map_qx3[12][0] + 4) / 14) + 2;

    // std::cout << path_size << std::endl;
    // std::vector<int> permute;
    // permute.push_back(12);
    // find_path(map_qx3, 12, 0, permute, path_qx3, path_size);
    // for (int i = 0; i < permute.size(); ++i){std::cout << " " << permute[i];}std::cout << std::endl;
    // return true;

    if ( is_set( "QS1_1" ) )
    {
        if(circ.lines() > 20)
        {
            std::cout << "Only up to 20 variables!" << std::endl;
            return true;
        }
        circ_qx = qxg(circ, map_qx20, path_qx20, statistics);
        print_runtime();
        circ_qx = optimize_circuit(circ_qx, statistics);
        print_runtime();
        std::cout << "After: " << circ_qx.num_gates() << std::endl;
    }
    else if ( is_set( "qx3" ) )
    {
        if(circ.lines() > 16)
        {
            std::cout << "Only up to 16 variables!" << std::endl;
            return true;
        }
        circ_qx = qxg(circ, map_qx3, path_qx3, statistics);
        std::cout << "Before: " << circ_qx.num_gates() << std::endl;
        print_runtime();
        circ_qx = optimize_circuit(circ_qx, statistics);
        std::cout << "After: " << circ_qx.num_gates() << std::endl;
        print_runtime();
    }
    else
    {
        if(circ.lines() > 5)
        {
            std::cout << "Only up to 5 variables! Try another option." << std::endl;
            return true;
        }
        if ( is_set( "qx4" ) )
            circ_qx = qxg(circ, map_qx4, map_qx4, statistics);
        else
            circ_qx = qxg(circ, map_qx2, map_qx2, statistics);
        print_runtime();
        circ_qx = optimize_circuit(circ_qx, statistics);
    }
    if ( is_set( "new" ) )
        circuits.extend();    
    circuits.current() = circ_qx;

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
