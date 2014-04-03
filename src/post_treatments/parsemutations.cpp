//
//  parsemutations.cpp
//  Aevol4.3local
//
//  Created by Antoine Frénoy on 02/04/2014.
//
//

#include "parsemutations.h"

void print_help( void );


int main(int argc, char** argv)
{
  // =====================
  //  Parse command line
  // =====================
  
  // Default values
  bool        verbose           = false;
  int32_t     begin_gener       = 0;
  int32_t     end_gener         = -1;
  char tree_file_name[50];
  
  const char * short_options = "hvncb:i:r:e:";
  static struct option long_options[] = {
    {"help",      no_argument,       NULL,  'h'},
    {"verbose",   no_argument,       NULL,  'v'},
    {"begin",     required_argument, NULL,  'b'},
    {"end",       required_argument,  NULL, 'e' },
    {0, 0, 0, 0}
  };
  
  int option;
  while( (option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1 )
  {
    switch( option )
    {
      case 'h' : print_help(); exit(EXIT_SUCCESS);  break;
      case 'v' : verbose = true;                    break;
      case 'b' : begin_gener  = atol(optarg);       break;
      case 'e' :
      {
        if ( strcmp( optarg, "" ) == 0 )
        {
          printf( "%s: error: Option -e or --end : missing argument.\n", argv[0] );
          exit( EXIT_FAILURE );
        }
        
        end_gener = atol( optarg );
        
        break;
      }
    }
  }
  
  if ( end_gener == -1 )
  {
    printf( "%s: error: You must provide a generation number.\n", argv[0] );
    exit( EXIT_FAILURE );
  }
  
  // Load the simulation
#ifndef __NO_X
  ae_exp_manager* exp_manager = new ae_exp_manager_X11();
#else
  ae_exp_manager* exp_manager = new ae_exp_manager();
#endif
  exp_manager->load( end_gener, false, true, false );
  
  int32_t tree_step = exp_manager->get_tree_step();
  int32_t nb_indivs = exp_manager->get_nb_indivs();
  int32_t nb_geners = end_gener - begin_gener;
  
  // Load the tree files and copy replication reports into big table reports
  ae_tree * tree = NULL;
  
  reports = new ae_replication_report**[nb_geners];
  int32_t i;
  for(i=0;i<end_gener - begin_gener;i++)
  {
    reports[i] = new ae_replication_report*[nb_indivs];
  }
  
  int32_t loading_step;
  for (loading_step = end_gener; loading_step >= begin_gener+tree_step; loading_step-=tree_step){
    sprintf( tree_file_name,"tree/tree_%06"PRId32".ae", loading_step );
    printf("  Loading tree file %s\n",tree_file_name);
    tree = new ae_tree( exp_manager, tree_file_name );
    int32_t generation,individual;
    for (generation=0;generation<tree_step;generation++) {
      for (individual=0;individual<nb_indivs;individual++){
        //printf("about to treat %"PRId32" %"PRId32" %"PRId32"\n",loading_step,generation,individual);
        reports[generation+loading_step-tree_step][individual]=new ae_replication_report(*(tree->get_report_by_index(generation+1,individual)));
      }
    }
    delete tree;
  }
  
  printf("Done with loading\n");

  // Open the output file
  char output_file_name[101];
  snprintf( output_file_name, 100, "mutations-b%06"PRId32"-e%06"PRId32".txt", begin_gener, end_gener);
  
  FILE* output_file = fopen(output_file_name, "w");
  if ( output_file == NULL )
  {
    fprintf(stderr, "File %s could not be created, exiting.\n", output_file_name);
    fprintf(stderr, "Please check your permissions in this directory.\n");
    exit(EXIT_FAILURE);
  }
  
  fprintf(output_file,"Small test\n");
  
  
  fclose(output_file);
  for(i=0;i<end_gener - begin_gener;i++)
  {
    delete [] reports[i];
  }
  delete [] reports;
  delete exp_manager;
  
  exit(EXIT_SUCCESS);
  
}



void print_help( void )
{
  printf( "Analysis and lineage of mutations post-treatment \n" );
  printf( "Usage : parsemutations [-b gen0 -e gen1] [-h]\n");
  printf( "\n" );
  printf( "\t-b gener1 or --begin gener1 : \n" );
  printf( "\t                  Retrieve the lineage up to generation gener1.\n" );
  printf( "\t                  There must be a genome backup file for this\n" );
  printf( "\t                  generation. If not specified, the program \n" );
  printf( "\t                  retrieves the lineage up to generation 0.\n");
  printf( "\n" );
  printf( "\t-e end_gener or --end end_gener : \n" );
  printf( "\t                  Retrieve the lineage of the individual of end_gener \n" );
  printf( "\n" );
  
}
