//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************


/** \class
 *  \brief
 */


#ifndef __AE_ALIGN_H__
#define __AE_ALIGN_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_vis_a_vis.h>
#include <ae_dna.h>




// =================================================================
//                          Class declarations
// =================================================================







class ae_align : public ae_object
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================

    // =================================================================
    //                             Destructors
    // =================================================================

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    static ae_vis_a_vis* search_alignement_direct( const ae_dna* chrom1, const int32_t seed1,
                                                   const ae_dna* chrom2, const int32_t seed2, int16_t needed_score );
    static ae_vis_a_vis* search_alignement_indirect( const ae_dna* chrom1, const int32_t seed1,
                                                     const ae_dna* chrom2, const int32_t seed2, int16_t needed_score );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_align( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_align( const ae_align &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_ALIGN_H__
