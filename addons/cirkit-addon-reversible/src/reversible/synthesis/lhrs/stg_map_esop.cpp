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

#include "stg_map_esop.hpp"

#include <boost/format.hpp>

#include <core/properties.hpp>
#include <core/utils/timer.hpp>
#include <classical/abc/gia/gia.hpp>
#include <classical/optimization/esop_minimization.hpp>
#include <classical/xmg/xmg_io.hpp>
#include <classical/xmg/xmg_simulate.hpp>
#include <reversible/functions/add_circuit.hpp>
#include <reversible/optimization/esop_post_optimization.hpp>
#include <reversible/synthesis/esop_synthesis.hpp>

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

void stg_map_esop( circuit& circ, const xmg_graph& function,
                   const std::vector<unsigned>& line_map,
                   const stg_map_esop_params& params,
                   stg_map_esop_stats& stats )
{
  /* using the `dumpfile' variable, we have the chance to write internal data structures into a file */
  if ( !params.dumpfile.empty() )
  {
    write_verilog( function, boost::str( boost::format( "%s/function-%d.aig" ) % params.dumpfile % stats.dumpfile_counter ) );
  }

  /* if we don't collapse, we don't synthesize */
  if ( params.nocollapse )
  {
    stats.dumpfile_counter++;
    return;
  }

  /* collapse AIG into initial ESOP
   *
   * note: we are using a "lambda trick" to initialize the esop directly and do some timing
   */
  auto esop = [&function, &params, &stats]() {
    increment_timer t( &stats.cover_runtime );

    Cudd mgr;
    xmg_bdd_simulator sim( mgr );
    const auto bdd = simulate_xmg_function( function, function.outputs().front().first, sim );

    /* get initial cover using exact PSDKRO optimization */
    exp_cache_t exp_cache;
    count_cubes_in_exact_psdkro( mgr.getManager(), bdd.getNode(), exp_cache );

    char * var_values = new char[mgr.ReadSize()];
    std::fill( var_values, var_values + mgr.ReadSize(), 2 );

    abc::Vec_Wec_t *esop = abc::Vec_WecAlloc( 0u );
    generate_exact_psdkro( mgr.getManager(), bdd.getNode(), var_values, -1, exp_cache, [&mgr, &bdd, &esop, &var_values]() {
        auto * level = abc::Vec_WecPushLevel( esop );
        for ( auto i = 0; i < mgr.ReadSize(); ++i )
        {
          if ( var_values[i] == 2 ) continue;
          abc::Vec_IntPush( level, ( i << 1u ) | !var_values[i] );
        }
        abc::Vec_IntPush( level, -1 );
      } );

    delete[] var_values;

    return gia_graph::esop_ptr( esop, &abc::Vec_WecFree );
  }();

  /* perform ESOP minimization, if we enabled it */
  if ( params.script != exorcism_script::none )
  {
    esop = [&esop, &function, &params, &stats]() {
      increment_timer t( &stats.exorcism_runtime );
      const auto em_settings = make_settings_from( std::make_pair( "progress", params.progress ), std::make_pair( "script", params.script ) );
      return exorcism_minimization( esop, function.inputs().size(), function.outputs().size(), em_settings );
    }();
  }

  /* we also use the `dumpfile' variable to write the resulting ESOP */
  if ( !params.dumpfile.empty() )
  {
    write_esop( esop, function.inputs().size(), function.outputs().size(),
                boost::str( boost::format( "%s/esop-%d.esop" ) % params.dumpfile % stats.dumpfile_counter++ ) );
  }

  /* finally, we perform ESOP based synthesis and possibly apply post optimization */
  if ( params.optimize_postesop )
  {
    circuit circ_local;
    esop_synthesis( circ_local, esop, function.inputs().size(), function.outputs().size() );
    auto circ_opt = esop_post_optimization( circ_local );
    append_circuit( circ, circ_opt, gate::control_container(), line_map );
  }
  else
  {
    const auto es_settings = make_settings_from( std::make_pair( "line_map", line_map ) );
    esop_synthesis( circ, esop, function.inputs().size(), function.outputs().size(), es_settings );
  }
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
