/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2016  EPFL
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "feathering.hpp"

#include <iostream>

#include <boost/any.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>

#include <core/utils/graph_utils.hpp>
#include <core/utils/timer.hpp>
#include <classical/functions/compute_levels.hpp>
#include <classical/utils/aig_utils.hpp>

namespace cirkit
{

/******************************************************************************
 * Types                                                                      *
 ******************************************************************************/

/******************************************************************************
 * Private functions                                                          *
 ******************************************************************************/

properties::ptr make_properties( const std::pair<std::string, boost::any>& pair )
{
  auto prop = std::make_shared<properties>();
  prop->set( pair.first, pair.second );
  return prop;
}

boost::dynamic_bitset<> complements_from_edges( const std::vector<aig_edge>& edges, const aig_graph& aig )
{
  const auto& complement = boost::get( boost::edge_complement, aig );
  boost::dynamic_bitset<> b( 2u );

  for ( const auto& e : edges )
  {
    b.set( complement[e] ? 1u : 0u );
  }

  return b;
}

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

aig_graph output_feathering( const aig_graph& aig, unsigned levels,
                             const properties::ptr& settings, const properties::ptr& statistics )
{
  using boost::format;
  using boost::str;

  /* settings */
  const auto respect_edges = get( settings, "respect_edges", true );
                            /**** description: ********************************************
                             * If true, then for each node within the feathering level,   *
                             * an output is added according to the polarities of outgoing *
                             * edges.  Otherwise, always two outputs are created for each *
                             * node, if not existing.                                     *
                             **************************************************************/
  const auto output_name   = get( settings, "output_name", std::string( "FO_%d_%d" ) );
  const auto verbose       = get( settings, "verbose", false );

  /* timing */
  properties_timer t( statistics );

  auto new_aig     = aig;
  const auto& info = aig_info( new_aig );

  const auto cl_statistics = std::make_shared<properties>();
  const auto vertex_levels = compute_levels( new_aig, make_properties( { "push_to_outputs", true } ), cl_statistics );
  const auto max_level     = cl_statistics->get<unsigned>( "max_level" );
  const auto in_edges      = precompute_ingoing_edges( new_aig );

  if ( verbose )
  {
    std::cout << "[i] output levels" << std::endl;
    for ( const auto& o : info.outputs )
    {
      std::cout << format( "[i] %s : %d" ) % o.second % vertex_levels.at( o.first.node ) << std::endl;
    }
    std::cout << format( "[i] max_level : %d" ) % max_level << std::endl;
  }

  for ( const auto& p : vertex_levels )
  {
    /* valid node? */
    if ( p.second + levels >= max_level && ( !respect_edges || p.second != max_level ) )
    {
      if ( verbose )
      {
        std::cout << "[i] add outputs to " << p.first << std::endl;
      }

      /* check for inverted output */
      auto complements = respect_edges ? complements_from_edges( in_edges.at( p.first ), new_aig ) : boost::dynamic_bitset<>( 2u, 3u );

      auto pos = complements.find_first();

      while ( pos != boost::dynamic_bitset<>::npos )
      {
        auto complemented = ( pos != 0u );

        /* check whether there is already an output? */
        auto exists = false;
        for ( const auto& output : info.outputs )
        {
          if ( output.first.node == p.first && output.first.complemented == complemented )
          {
            if ( verbose )
            {
              std::cout << "[i] output already exists for polarity " << !complemented << std::endl;
            }
            exists = true;
            break;
          }
        }

        /* add only if there is no output */
        if ( !exists )
        {
          aig_create_po( new_aig, {p.first, complemented}, str( format( output_name ) % p.first % p.second ) );
        }

        pos = complements.find_next( pos );
      }
    }
  }

  return new_aig;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
