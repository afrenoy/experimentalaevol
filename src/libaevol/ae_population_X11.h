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
 
 
 #ifndef __AE_POPULATION_X11_H__
#define  __AE_POPULATION_X11_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_population.h>
#include <ae_X11_window.h>




// =================================================================
//                          Class declarations
// =================================================================





 
class ae_population_X11 : public ae_population
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_population_X11( ae_exp_manager* exp_m );
    //~ ae_population_X11( gzFile backup_file );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_population_X11( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    virtual void display( ae_X11_window* win );
    virtual void display_grid( ae_X11_window* win, double** cell_grid );

    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_population_X11( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_population_X11( const ae_population &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  
    // =================================================================
    //                           Protected Methods
    // =================================================================
    void compute_colormap(void);  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    char ** _col_map;

};


// =====================================================================
//                          Accessors definitions
// =====================================================================


#endif // __AE_POPULATION_X11_H__