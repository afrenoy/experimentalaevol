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
  // Parse the parameters
  int32_t     begin_gener       = 0;
  int32_t     end_gener         = -1;
  char tree_file_name[50];
  
  const char * short_options = "hb:e:";
  static struct option long_options[] = {
    {"help",      no_argument,       NULL,  'h'},
    {"begin",     required_argument, NULL,  'b'},
    {"end",       required_argument,  NULL, 'e' },
    {0, 0, 0, 0}
  };
  
  int option;
  while( (option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1 )
  {
    switch( option )
    {
      case 'h' : print_help(); exit(EXIT_SUCCESS); break;
      case 'b' : begin_gener = atol(optarg); break;
      case 'e' : end_gener = atol(optarg); break;
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
  //printf("%"PRId32"\n",tree_step);
  for (loading_step = end_gener; loading_step >= begin_gener+tree_step; loading_step-=tree_step){
    sprintf( tree_file_name,"tree/tree_%06"PRId32".ae", loading_step );
    printf("  Loading tree file %s\n",tree_file_name);
    tree = new ae_tree( exp_manager, tree_file_name );
    int32_t generation,individual;
    for (generation=0;generation<tree_step;generation++) {
      for (individual=0;individual<nb_indivs;individual++){
        //printf("about to treat %"PRId32" %"PRId32" %"PRId32"\n",loading_step,generation,individual);
        reports[generation+loading_step-tree_step-begin_gener][individual]=new ae_replication_report(*(tree->get_report_by_index(generation+1,individual)));
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
  
  // Calculate reproductive success of each (set of) mutation
  reproductive_success=new int32_t*[nb_geners];
  for(i=0;i<nb_geners;i++)
  {
    reproductive_success[i] = new int32_t[nb_indivs];
  }
  
  int32_t generation,individual;
  
  // Dynamic programming algorithm
  // initialize the table with 1s
  for (generation=0;generation<nb_geners;generation++) {
    for (individual=0;individual<nb_indivs;individual++){
      reproductive_success[generation][individual]=0;
    }
  }
  /*
  // fill the last generation with 0s
  for (individual=0;individual<nb_indivs;individual++){
    reproductive_success[nb_geners-1][individual]=0;
  }
  */
  // then go backward from last generation
  for (generation=nb_geners-1;generation>0;generation--) {
    for (individual=0;individual<nb_indivs;individual++){
      // Find parent and update its reproductive success
      int32_t idparent = reports[generation][individual]->get_parent_id(); // the report that describes creation of individual at generation + 1 from parent_id at generation.
      reproductive_success[generation-1][idparent]+=reproductive_success[generation][individual]+1;
    }
  }
  // Note that at this stage we do not use the "first" reports (reports[0][x]) because they tell us where generation 1 come from, so what of the individuals of generation 0 are successful, but we do not know what mutations permitted to obtain generation 0 from generation -1.
  // We will however use these reports when we will analyse reproductive_success[0][x] (individuals at generation 1) because these reports explain what created this generation of individuals.
  
  // Assert the sum of reproductive success is correct at each generation
  for (generation=0;generation<nb_geners;generation++) {
    int32_t sum=0;
    for (individual=0;individual<nb_indivs;individual++){
      sum+=reproductive_success[generation][individual];
    }
    //printf("%"PRId32" %"PRId32" %"PRId32"\n",generation,sum,(nb_geners-1-generation)*nb_indivs);
    assert(sum==(nb_geners-1-generation)*nb_indivs);
  }
  
  // output all the analyzed mutations
  for (generation=0;generation<nb_geners;generation++){
    for (individual=0;individual<nb_indivs;individual++){
      double metabolic_effect=reports[generation][individual]->get_metabolic_error() - reports[generation][individual]->get_parent_metabolic_error();
      double secretion_effect=reports[generation][individual]->get_secretion_error() - reports[generation][individual]->get_parent_secretion_error();
      // negative value = smaller gap = beneficial mutation
      fprintf(output_file,"%"PRId32" %"PRId32" %"PRId32" %+.15f %+.15f\n", generation+begin_gener, individual, reproductive_success[generation][individual], metabolic_effect, secretion_effect);
    }
  }
  
  
  // Close the files and clean memory
  fclose(output_file);
  for(i=0;i<end_gener - begin_gener;i++)
  {
    delete [] reports[i];
    delete [] reproductive_success[i];
  }
  delete [] reports;
  delete [] reproductive_success;
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
