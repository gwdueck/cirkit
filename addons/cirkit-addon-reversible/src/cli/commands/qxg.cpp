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
static const matrix map_qx2 = { {0,0,0,10,10},
                                {4,0,0,10,10},
                                {4,4,0,4,4},
                                {10,10,0,0,0},
                                {10,10,0,4,0}};
//Matrix with the cost of each possible cnot (QX4)
static const matrix map_qx4 = { {0,4,4,10,10},
                                {0,0,4,10,10},
                                {0,0,0,4,0},
                                {12,12,0,0,0},
                                {10,10,4,4,0}};
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

static const matrix path_qx5={{0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1},
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

static const matrix map_qx5={{0,4,10,20,32,38,48,60,66,60,54,42,30,20,10,4},
                            {0,0,0,12,24,34,44,56,66,72,66,54,42,32,22,10},
                            {10,4,0,0,12,22,32,44,54,64,54,42,30,20,10,4},
                            {20,10,4,0,0,10,20,32,42,52,40,32,20,10,0,10},
                            {30,20,10,4,0,4,10,20,30,40,28,20,10,4,10,20},
                            {42,30,20,10,0,0,4,10,20,30,18,10,4,10,22,30},
                            {50,38,28,18,12,0,0,0,10,20,12,0,10,18,30,38},
                            {60,48,38,28,22,10,4,0,4,10,0,10,20,28,40,48},
                            {72,60,50,40,34,22,10,0,0,4,10,20,30,38,46,60},
                            {62,70,60,50,44,32,20,10,0,0,0,10,20,32,40,50},
                            {50,58,48,38,32,20,10,4,10,4,0,4,10,20,28,38},
                            {40,48,38,28,22,10,4,10,20,10,0,0,4,10,18,28},
                            {34,38,28,18,12,0,10,22,32,18,12,0,0,0,12,22},
                            {22,30,20,10,0,10,20,32,42,28,22,10,4,0,0,10},
                            {10,20,10,4,10,20,30,38,48,38,32,20,10,4,0,4},
                            {0,10,0,10,22,28,38,50,56,50,44,32,20,10,0,0}};

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
     ( "qx5,5", "IBM QX5 matrix")
     ( "QS1_1,20", "QS1_1 matrix")
    ;
  add_new_option();
}


command::rules_t qxg_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}


// Create cnots matrix and cost matrix
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

// Permute the matrix
matrix permute_matrix(matrix& cnots, const std::map<int, int>& permutation, matrix& aux)
{
    for(auto it : permutation)
        for (int j = 0; j < cnots.size(); ++j)
            aux[it.second][j] = cnots[it.first][j];
    for(auto it : permutation)
        for (int i = 0; i < cnots.size(); ++i)
            cnots[i][it.second] = aux[i][it.first];
    return cnots;
}

// Get the permutation cost
unsigned permute_cost(const matrix& cnots, const matrix& map)
{
	unsigned int cost = 0;
	for (int i = 0; i < cnots.size(); ++i)
		for (int j = 0; j < cnots.size(); ++j)
			cost += cnots[i][j] * map[i][j];
	return cost;
}

// Find the higher value in the matrix
std::pair<int,int> get_position_higher_value_matrix(matrix& m)
{
    unsigned int h = 0, x, y;
    for (int i = 0; i < m.size(); ++i)
    {
        for (int j = 0; j < m.size(); ++j)
        {
            if( m[i][j] >= h )
            {
                h = m[i][j];
                x = i;
                y = j;
            }
        }
    }
    return std::make_pair (x,y);
}

// Get a mapping from the matrix
std::pair<int,int> get_mapping(const matrix& map, std::pair<int,int>& qubit, std::map<int, int>& permutation, std::vector<int> p)
{
    std::map<int, int>::iterator p1, p2;
    p1 = permutation.find(qubit.first);
    p2 = permutation.find(qubit.second);
    if(p1 != permutation.end() && p2 != permutation.end())
    {
        //Do nothing
        // std::cout << "Do nothing" << std::endl;
        return std::make_pair(p1->second,p2->second);
    }
    else if(p1 == permutation.end() && p2 == permutation.end())
    {
        unsigned int l = 1000, x, y;
        for (int i = 0; i < map.size(); ++i)
        {
            for (int j = 0; j < map.size(); ++j)
            {
                if(i != j && map[i][j] < l && std::find(p.begin(), p.end(), i) == p.end() && std::find(p.begin(), p.end(), j) == p.end())
                {
                    l = map[i][j];
                    x = i;
                    y = j;
                }
            }
        }
        // std::cout << "nenhum " << x << " " << y << std::endl;
        return std::make_pair(x,y);
    }
    else if(p1 == permutation.end() && p2 != permutation.end())
    {
        unsigned int l = 1000, x, y;
        y = p2->second;
        for(int i = 0; i < map.size(); ++i)
        {
            if(map[i][y] < l && std::find(p.begin(), p.end(), i) == p.end())
            {
                l = map[i][y];
                x = i;
            }    
        }
        // std::cout << "sem controle " << x << " " << y << std::endl;
        return std::make_pair(x,y);
    }
    else if(p1 != permutation.end() && p2 == permutation.end())
    {
        unsigned int l = 1000, x, y;
        x = p1->second;
        for(int j = 0; j < map.size(); ++j)
        {
            if(map[x][j] < l && std::find(p.begin(), p.end(), j) == p.end())
            {
                l = map[x][j];
                y = j;
            }    
        }
        // std::cout << "sem alvo " << x << " " << y << std::endl;
        return std::make_pair(x,y);
    }
    std::cout << "====> should not happen" << std::endl;
    return std::make_pair(0,0); // added -- not reachable???
}

// Map the last qubit
std::map<int, int> complete_permutation(std::map<int, int>& permutation, unsigned int s)
{
    std::vector<int> x, y;
    for (int i = 0; i < s; ++i)
    {
        x.push_back(i);
        y.push_back(i);
    }
    for(auto it : permutation)
    {
        x.erase(std::remove(x.begin(), x.end(), it.first), x.end());
        y.erase(std::remove(y.begin(), y.end(), it.second), y.end());
    }
    for (int i = 0; i < x.size(); ++i)
        permutation.insert(std::pair<int, int>(x[i], y[i]));
    return permutation;
}

// To duplicate a matrix
matrix copy_matrix(matrix& aux, const matrix& m)
{
	for (int i = 0; i < m.size(); ++i)
        for (int j = 0; j < m.size(); ++j)
            aux[i][j] = m[i][j];
    return aux;
}


// Main function -- Find all the mappings
circuit qxg(circuit& circ, const matrix& map, const matrix& path, properties::ptr& statistics )
{
    properties_timer t( statistics );
    circuit circ_qx;
    unsigned int cost;
    std::vector <int> p;
    std::pair<int,int> qubit1, qubit2;
    std::map<int, int> permutation;
    matrix cnots, map_cost, aux;

    // Set the matrices
    for (int i = 0; i < map.size(); ++i)
        p.push_back(0);
    
    for (int i = 0; i < map.size(); ++i)
    {
        cnots.push_back(p);
        map_cost.push_back(p);
        aux.push_back(p);
    }
    p.clear();

    // Get the default permutation cost
    cost = initial_matrix(circ, cnots, map_cost, map);
    std::cout << "number of gates: " << circ.num_gates() << std::endl;
    std::cout << "default permutation: " << cost << std::endl;
    
    // Copy matrix to aux
    aux = copy_matrix(aux, cnots);
    
    for(int i = 0; i < 16; i++){
        for(int j = 0; j < 16; j++) std::cout << aux[i][j] << " ";
        std::cout << std::endl;
    }
    
    // Map the qubits until 1 left    
    while(permutation.size() < cnots.size() - 1)
    {   
        // Get higher CNOT count 
        qubit1 = get_position_higher_value_matrix(aux);
        
        // If the higher qubit is zero... (maybe this can be changed)
        if(aux[qubit1.first][qubit1.second] == 0) break;
        
        std::cout << "mapping, cost " << qubit1.first << " " << qubit1.second<< " " << aux[qubit1.first][qubit1.second] << std::endl;
        
        aux[qubit1.first][qubit1.second] = 0;
        
        // Get mapping for the second qubit
        qubit2 = get_mapping(map, qubit1, permutation, p);
        
        std::cout << qubit1.first << "  " << qubit1.second << " map to " << qubit2.first << "  " << qubit2.second << std::endl;
        
        // Saving the mapping
        permutation.insert(std::pair<int, int>(qubit1.first, qubit2.first));
        permutation.insert(std::pair<int, int>(qubit1.second, qubit2.second));
        
        // Vector to save the qubits already mapped
        p.push_back(qubit2.first);
        p.push_back(qubit2.second);
    }

    // Finish the mapping
    permutation = complete_permutation(permutation, cnots.size());
    
    // Permute the matrix with the CNOTs
    cnots = permute_matrix(cnots, permutation, aux);
    
    // Get the cost of the permutation
    cost = permute_cost(cnots, map);
    
    // Print permutation
    std::cout << "final permutation: " << cost << std::endl;
    std::cout << "total gates: " << cost + circ.num_gates() << std::endl;
    for(auto it : permutation)
    	std::cout << " " << it.second;
    std::cout << std::endl;
    return circ_qx;
}

bool qxg_command::execute()
{
    auto& circuits = env->store<circuit>();
    circuit aux = circuits.current();
    circuit circ_qx, circ;
    copy_circuit(aux, circ);    
    //auto settings = make_settings();

    
    if ( is_set( "QS1_1" ) )
    {
        if(circ.lines() > 20)
        {
            std::cout << "Only up to 20 variables!" << std::endl;
            return true;
        }
        circ_qx = qxg(circ, map_qx20, path_qx20, statistics);
        // print_runtime();
    }
    else if ( is_set( "qx3" ) )
    {
        if(circ.lines() > 16)
        {
            std::cout << "Only up to 16 variables!" << std::endl;
            return true;
        }
        circ_qx = qxg(circ, map_qx3, map_qx3, statistics);
        // print_runtime();

    }
    else if ( is_set( "qx5" ) )
    {
        if(circ.lines() > 16)
        {
            std::cout << "Only up to 16 variables!" << std::endl;
            return true;
        }
        circ_qx = qxg(circ, map_qx5, map_qx5, statistics);
        print_runtime();
        // std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
        // std::cout.precision(9);
        // std::cout << "\t" << statistics->get<double>( "runtime" ) << std::endl;
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
        // print_runtime();
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
