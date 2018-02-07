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

#include "stg_map_luts.hpp"

#include <core/utils/range_utils.hpp>
#include <core/utils/timer.hpp>
#include <classical/utils/truth_table_utils.hpp>
#include <classical/xmg/xmg_cover.hpp>
#include <classical/xmg/xmg_extract.hpp>
#include <classical/xmg/xmg_flow_map.hpp>
#include <classical/xmg/xmg_simulate.hpp>
#include <reversible/functions/add_gates.hpp>

namespace cirkit
{

/******************************************************************************
 * Functions                                                                  *
 ******************************************************************************/

tt xmg_truth_table_from_lut_rec( const xmg_graph& xmg, xmg_node node, std::unordered_map<xmg_node, tt>& node_to_tt )
{
  const auto it = node_to_tt.find( node );
  if ( it != node_to_tt.end() )
  {
    return it->second;
  }

  assert( !xmg.is_input( node ) );

  tt f;

  if ( xmg.is_maj( node ) )
  {
    auto children = xmg.children( node );
    tt tt0;

    if ( children[0].node == 0 ) /* AND or OR */
    {
      tt0 = tt_const0();
    }
    else
    {
      tt0 = xmg_truth_table_from_lut_rec( xmg, children[0].node, node_to_tt );
    }

    const auto tt1 = xmg_truth_table_from_lut_rec( xmg, children[1].node, node_to_tt );
    const auto tt2 = xmg_truth_table_from_lut_rec( xmg, children[2].node, node_to_tt );

    f = tt_maj( children[0].complemented ? ~tt0 : tt0, children[1].complemented ? ~tt1 : tt1, children[2].complemented ? ~tt2 : tt2 );
  }
  else if ( xmg.is_xor( node ) )
  {
    auto children = xmg.children( node );
    const auto tt0 = xmg_truth_table_from_lut_rec( xmg, children[0].node, node_to_tt );
    const auto tt1 = xmg_truth_table_from_lut_rec( xmg, children[1].node, node_to_tt );

    f = ( children[0].complemented ? ~tt0 : tt0 ) ^ ( children[1].complemented ? ~tt1 : tt1 );
  }
  else
  {
    assert( false );
  }

  node_to_tt[node] = f;

  return f;
}

tt xmg_truth_table_from_lut( const xmg_graph& xmg, xmg_node root )
{
  assert( xmg.has_cover() && xmg.cover().has_cut( root ) && xmg.cover().num_leafs( root ) <= 6u );

  std::unordered_map<xmg_node, tt> node_to_tt;

  auto i = 0u;
  for ( const auto& leaf : xmg.cover().cut( root ) )
  {
    node_to_tt[leaf] = tt_nth_var( i++ );
  }

  return xmg_truth_table_from_lut_rec( xmg, root, node_to_tt );
}

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

class stg_map_luts_impl
{
public:
  stg_map_luts_impl( circuit& circ, const xmg_graph& function,
                     const std::vector<unsigned>& line_map,
                     const std::vector<unsigned>& ancillas,
                     const stg_map_luts_params& params,
                     stg_map_luts_stats& stats )
    : circ( circ ),
      function( function ),
      line_map( line_map ),
      ancillas( ancillas ),
      params( params ),
      stats( stats )
  {
  }

  void run()
  {
    assert( !function.has_cover() );

    const int num_inputs = function.inputs().size();

    /* if very small fall back to precomputed database */
    if ( num_inputs <= params.max_cut_size )
    {
      xmg_tt_simulator sim;
      auto tt = simulate_xmg_function( function, function.outputs().front().first, sim );
      stg_map_precomp( circ, tt.to_ulong(), num_inputs, line_map, *params.map_precomp_params, *stats.map_precomp_stats );
      return;
    }

    const auto mapping = compute_mapping();

    if ( !mapping.has_cover() )
    {
      stg_map_esop( circ, mapping, line_map, *params.map_esop_params, *stats.map_esop_stats );
      return;
    }

    /* LUT based mapping */
    auto pi_index = 0u;
    auto anc_index = 0u;
    auto ins_index = 0u;
    std::vector<unsigned> lut_to_line( mapping.size() );
    std::vector<unsigned> synth_order( 2 * ( mapping.cover().lut_count() - 1 ) + 1, 0 );

    for ( const auto& input : mapping.inputs() )
    {
      lut_to_line[input.first] = line_map[pi_index++];
    }

    auto root = mapping.outputs().front().first.node;
    for ( const auto& index : mapping.topological_nodes() )
    {
      if ( !mapping.cover().has_cut( index ) ) continue;

      if ( index == root )
      {
        lut_to_line[index] = line_map[pi_index++];
        synth_order[ins_index] = index;
      }
      else
      {
        lut_to_line[index] = ancillas[anc_index++];
        synth_order[ins_index] = synth_order[synth_order.size() - 1 - ins_index] = index;
        ++ins_index;
      }
    }

    std::unordered_map<unsigned, std::pair<unsigned, unsigned>> esop_circ_cache;
    std::unordered_map<unsigned, uint64_t> truth_table_cache;

    for ( auto index : synth_order )
    {
      const int num_inputs = mapping.cover().num_leafs( index );
      std::vector<unsigned> local_line_map;
      local_line_map.reserve( num_inputs + 1u );

      for ( auto fanin : mapping.cover().cut( index ) )
      {
        local_line_map.push_back( lut_to_line[fanin] );
      }
      local_line_map.push_back( lut_to_line[index] );

      if ( num_inputs == 0 )
      {
        assert( false );
      }
      else if ( num_inputs == 1 )
      {
        assert( false );
        //assert( mapping.lut_truth_table( index ) == 1 );
        append_cnot( circ, make_var( local_line_map[0], false ), local_line_map[1] );
      }
      else if ( num_inputs <= params.max_cut_size )
      {
        uint64_t func;
        const auto it = truth_table_cache.find( index );
        if ( it == truth_table_cache.end() )
        {
          func = xmg_truth_table_from_lut( mapping, index ).to_ulong();
          truth_table_cache[index] = func;
        }
        else
        {
          func = it->second;
        }
        stg_map_precomp( circ, func, num_inputs, local_line_map, *params.map_precomp_params, *stats.map_precomp_stats );
      }
      else
      {
        const auto it = esop_circ_cache.find( index );
        if ( it == esop_circ_cache.end() )
        {
          if ( params.map_esop_params->progress )
          {
            std::cout << "\n";
          }
          const auto lut = xmg_extract_lut( mapping, index );

          const auto begin = circ.num_gates();
          stg_map_esop( circ, lut, local_line_map, *params.map_esop_params, *stats.map_esop_stats );
          esop_circ_cache.insert( {index, {begin, circ.num_gates()}} );

          if ( params.map_esop_params->progress )
          {
            std::cout << "\e[A";
          }
        }
        else
        {
          for ( auto i = it->second.first; i < it->second.second; ++i )
          {
            circ.append_gate() = circ[i];
          }
        }
      }
    }
  }

private:
  xmg_graph compute_mapping()
  {
    switch ( params.strategy )
    {
    case stg_map_luts_params::mapping_strategy::mindb:
      return compute_mapping_mindb();
    case stg_map_luts_params::mapping_strategy::bestfit:
      return compute_mapping_bestfit();
    default:
      assert( false );
    }
  }

  xmg_graph compute_mapping_mindb()
  {
    const auto num_inputs = function.inputs().size();

    for ( auto k = 3; k <= params.max_cut_size; ++k )
    {
      xmg_graph mapping = function;
      xmg_flow_map( mapping, make_settings_from( std::make_pair( "cut_size", static_cast<unsigned>( k ) ) ) );

      if ( mapping.cover().lut_count() - 1 <= ancillas.size() )
      {
        /* mapping is small enough, we're lucky */
        return mapping;
      }
      else if ( k == params.max_cut_size )
      {
        /* last round and no fit, we need to merge some LUTs */
        auto num_ancilla = mapping.cover().lut_count() - 1; /* no need to store the output LUT */

        if ( ancillas.empty() )
        {
          /* mapping won't help, so we return our function */
          return function;
        }

        while ( num_ancilla > ancillas.size() )
        {
          //abc::Gia_ManMergeTopLuts( mapping );
          assert( false );
          --num_ancilla;
        }

        /* now we merged top LUTs, but the big LUT can now have more inputs than our function */
        auto root = mapping.outputs().front().first.node;
        if ( mapping.cover().num_leafs( root ) > num_inputs )
        {
          /* mapping is not better than function */
          return function;
        }

        return mapping;
      }
      /* else we continue */
    }

    assert( false );
  }

  xmg_graph compute_mapping_bestfit()
  {
    const auto num_inputs = function.inputs().size();

    for ( auto k = 4u; k < num_inputs; ++k )
    {
      xmg_graph mapping = function;
      xmg_flow_map( mapping, make_settings_from( std::make_pair( "cut_size", k ) ) );

      /* Do we have enough ancillas? */
      if ( mapping.cover().lut_count() - 1 <= ancillas.size() )
      {
        return mapping;
      }
    }

    /* if no mapping fits, we return the function */
    return function;
  }

private:
  circuit& circ;
  const xmg_graph& function;
  const std::vector<unsigned>& line_map;
  const std::vector<unsigned>& ancillas;
  const stg_map_luts_params& params;
  stg_map_luts_stats& stats;
};

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

void stg_map_luts( circuit& circ, const xmg_graph& function,
                   const std::vector<unsigned>& line_map,
                   const std::vector<unsigned>& ancillas,
                   const stg_map_luts_params& params,
                   stg_map_luts_stats& stats )
{
  stg_map_luts_impl impl( circ, function, line_map, ancillas, params, stats );
  impl.run();
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
