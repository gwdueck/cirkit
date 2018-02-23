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
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/functions/remove_dup_gates.hpp>

typedef std::vector<std::vector<int>> matrix;

//Matrix with the cost of each possible cnot (QX2)
static const matrix map_qx2 = {{0,0,0,10,10}, {4,0,0,10,10}, {4,4,0,4,4}, {10,10,0,0,0}, {10,10,0,4,0}};
//Matrix with the cost of each possible cnot (QX4)
static const matrix map_qx4 = {{0,4,4,10,10}, {0,0,4,10,10}, {0,0,0,4,0}, {10,10,0,0,0}, {10,10,4,4,0}};
//Matrix with the cost of each possible cnot (QX3)
static const matrix map_qx3 = {{0,0,10,10,10,10,10,10,10,10,10,10,10,10,10,4},
                    {4,0,0,10,10,10,10,10,10,10,10,10,10,10,10,10},
                    {10,4,0,0,10,10,10,10,10,10,10,10,10,10,10,10},
                    {10,10,4,0,4,10,10,10,10,10,10,10,10,10,0,10},
                    {10,10,10,0,0,0,10,10,10,10,10,10,10,4,10,10},
                    {10,10,10,10,4,0,10,10,10,10,10,10,4,10,10,10},
                    {10,10,10,10,10,10,0,0,10,10,10,0,10,10,10,10},
                    {10,10,10,10,10,10,4,0,4,10,0,10,10,10,10,10},
                    {10,10,10,10,10,10,10,0,0,4,10,10,10,10,10,10},
                    {10,10,10,10,10,10,10,10,0,0,0,10,10,10,10,10},
                    {10,10,10,10,10,10,10,4,10,4,0,4,10,10,10,10},
                    {10,10,10,10,10,10,4,10,10,10,0,0,4,10,10,10},
                    {10,10,10,10,10,0,10,10,10,10,10,0,0,0,10,10},
                    {10,10,10,10,0,10,10,10,10,10,10,10,4,0,0,10},
                    {10,10,4,10,10,10,10,10,10,10,10,10,10,4,0,4},
                    {0,10,10,10,10,10,10,10,10,10,10,10,10,10,0,0}};



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
    ;
  add_new_option();
}


command::rules_t qxg_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

//Change two lines of the circuit
void manipulate_matrix( matrix& m1, int x, int y, matrix& m2, std::vector<int>& p, const matrix& map )
{
    int aux;

    aux = p[x];
    p[x] = p[y];
    p[y] = aux;
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

//Search matrix for the qubit with higher cost
int higher_cost( const matrix& m1, const matrix& m2, const std::vector<int>& p)
{
    int cost, higher_cost = 0, index = 0, higher_qtd_cnot = 0, qtd_cnot;
    
    for(int i=0; i<m1.size(); ++i)
    {
        cost = 0;
        qtd_cnot = 0;
        if (std::find(p.begin(), p.end(), i) == p.end())
        {
            for(int j=0; j<m1.size(); ++j)
            {
                cost += m1[i][j] + m1[j][i];
                qtd_cnot += m2[i][j] + m2[j][i];
            }
            if(cost > higher_cost)
            {
                higher_cost = cost;
                higher_qtd_cnot = qtd_cnot;
                index = i;
            }
            else if(cost == higher_cost && qtd_cnot > higher_qtd_cnot)
            {
                higher_qtd_cnot = qtd_cnot;
                index = i;
            }
        }
            
    }
    return index;
}

//Create cnots matrix and cost matrix
unsigned initial_matrix(circuit circ, matrix& cnots, matrix& map_cost, const matrix& map)
{
    unsigned target, control;
    unsigned cost = 0;
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
unsigned matrix_cost(const matrix& m)
{
    unsigned cost = 0;
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

void print_results( const matrix& cnots, const std::vector<int>& perm, const unsigned cost)
{
    std::cout << "Best permutation found (gates =  " << cost << "):";
    print_permutation(perm);
    
    //This print the permutation like the ibm command
    std::cout << "ibm representation: ";
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

bool qxg_command::execute()
{
    auto& circuits = env->store<circuit>();
    circuit circ = circuits.current();
    unsigned cost, lower_cost;
    int h, c, q;
    std::vector <int> p;
    std::vector<int> perm;
    std::vector<int> best_perm;

    if ( is_set( "qx3" ) )
    {
        matrix cnots(16, std::vector<int>(16));
        matrix map_cost(16, std::vector<int>(16));
        for (int i = 0; i < cnots.size(); ++i)
        {           
            perm.push_back(i);
            best_perm.push_back(i);
        }
        cost = initial_matrix(circ, cnots, map_cost, map_qx3);
        cost = cost + circ.num_gates();
        std::cout << circ.num_gates() << std::endl;
        std::cout << "initial cost: " << cost << std::endl;
        lower_cost = cost;
        std::cout << "cnots matrix: " << std::endl;
        print_matrix(cnots);
        std::cout << "cost matrix: " << std::endl;
        print_matrix(map_cost);

        unsigned int it = 0;
        do
        {
            h = higher_cost(map_cost, cnots, p);
            q = h;
            for(int i=0; i<cnots.size(); ++i)
            {
                manipulate_matrix( cnots, h, i, map_cost, perm, map_qx3 );
                cost = matrix_cost(map_cost) + circ.num_gates();
                if(cost <= lower_cost)
                {
                    q = i;
                    lower_cost = cost;
                    for(int j=0; j<cnots.size(); ++j)
                        best_perm[j] = perm[j];
                }
                manipulate_matrix( cnots, i, h, map_cost, perm, map_qx3 );
            }
            p.push_back(q);
            manipulate_matrix( cnots, h, q, map_cost, perm, map_qx3 );
            it++;
        } while (it < cnots.size());
        print_results(cnots, best_perm, lower_cost);
    }
    else
    {
        matrix cnots(5, std::vector<int>(5));
        matrix map_cost(5, std::vector<int>(5));
        for (int i = 0; i < cnots.size(); ++i)
        {           
            perm.push_back(i);
            best_perm.push_back(i);
        }

        if ( is_set( "qx4" ) )
            cost = initial_matrix(circ, cnots, map_cost, map_qx4);
        else
            cost = initial_matrix(circ, cnots, map_cost, map_qx2);

        cost = cost + circ.num_gates();
        //std::cout << circ.num_gates() << std::endl;
        //std::cout << "initial cost: " << cost << std::endl;
        lower_cost = cost;
        //std::cout << "cnots matrix: " << std::endl;
        //print_matrix(cnots);
        //std::cout << "cost matrix: " << std::endl;
        //print_matrix(map_cost);

        unsigned int it = 0;
        do
        {
            //std::cout << "Perm: ";
            //print_permutation(perm);
            h = higher_cost(map_cost, cnots, p);
            //std::cout << "qubit higher cost: " << h << std::endl;
            if ( is_set( "qx4" ) )
            {
                q = h;
                //print_permutation(perm);
                for(int i=0; i<cnots.size(); ++i)
                {
                    manipulate_matrix( cnots, h, i, map_cost, perm, map_qx4);
                    //print_permutation(perm);
                    cost = matrix_cost(map_cost) + circ.num_gates();
                    //std::cout << " i: " << i << " " << cost << std::endl;
                    if(cost <= lower_cost)
                    {
                        q = i;
                        lower_cost = cost;
                        for(int j=0; j<cnots.size(); ++j)
                            best_perm[j] = perm[j];
                    }
                    manipulate_matrix( cnots, i, h, map_cost, perm, map_qx4 );
                }
                p.push_back(q);
                manipulate_matrix( cnots, h, q, map_cost, perm, map_qx4 );
            }
            else
            {
                q = h;
                //print_permutation(perm);
                for(int i=0; i<cnots.size(); ++i)
                {
                    manipulate_matrix( cnots, h, i, map_cost, perm, map_qx2 );
                    //print_permutation(perm);
                    cost = matrix_cost(map_cost) + circ.num_gates();
                    //std::cout << " i: " << i << " " << cost << std::endl;
                    if(cost <= lower_cost)
                    {
                        q = i;
                        lower_cost = cost;
                        for(int j=0; j<cnots.size(); ++j)
                            best_perm[j] = perm[j];
                    }
                    manipulate_matrix( cnots, i, h, map_cost, perm, map_qx2 );
                }
                p.push_back(q);
                manipulate_matrix( cnots, h, q, map_cost, perm, map_qx2 );              
            }
            
            it++;
        } while (it < cnots.size()); 
        
        print_results(cnots, best_perm, lower_cost);
    }

    //piece of code from ibm.cpp
    unsigned start = circ.lines()+1;
    circuit circ_qx;
    for(unsigned i = start ; i <= 5u; i++)
    {
        add_line_to_circuit( circ, "i" + boost::lexical_cast<std::string>(i) , "o" + boost::lexical_cast<std::string>(i));
    }
    circuits.extend();
    circuits.current() = circ;

    if ( is_set( "qx4" ) )
    {

    }
    else if ( is_set( "qx3" ) )
    {

    }
    else 
    { 
        
        //circ_qx = remove_dup_gates( circ_qx );
        //circuits.extend();
        //circuits.current() = circ_qx;
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
