// This program adds a unique integer to each individual. It will be inherited by descendants.

 
// =================================================================
//                              Libraries
// =================================================================
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>


// =================================================================
//                            Project Files
// =================================================================
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_list.h>
#include <ae_exp_manager.h>




void print_help( char* prog_name );
void analyse_indiv( ae_individual* indiv, FILE* triangles_file, FILE* sequence_file, int16_t gu, ae_environment* env );
void analyse_gu( ae_genetic_unit* gen_unit, int32_t gen_unit_number, FILE* triangles_file, ae_environment* env );



int main( int argc, char* argv[] )
{
  // Initialize command-line option variables with default values  
  char* pop_file_name  = NULL;
  char* triangles_file_name  = NULL;
  char* sequence_file_name  = NULL;
  bool best_only = false;
  int16_t gu = -1;
  int32_t num_gener = -1;
  
  // Define allowed options
  const char * options_list = "hr:";
  static struct option long_options_list[] = {
    { "help", 1, NULL, 'h' },
    { "resume", 1, NULL, 'r' },
    { 0, 0, 0, 0 }
  };

  // Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 ) 
  {
    switch ( option )
    {
      case 'h' :
        print_help( argv[0] );
        exit( EXIT_SUCCESS );
        break;
      case 'r':
        num_gener = atol( optarg );
        break;  
    }
  }
  
  // Open the files

  ae_population* pop = NULL;
  ae_exp_manager* exp_manager = new ae_exp_manager();
  
  // Two possible sources: either the user provided a "full" simulation via a generation number (option '-r'), either he just provided a population file (option '-p').
  if ( num_gener == -1 )
  {
    printf("You must specify a generation number");
    exit(EXIT_FAILURE);
  }
  exp_manager->load( num_gener, false, false, false );
  pop = exp_manager->get_pop();

  ae_list_node<ae_individual*>* indiv_node = pop->get_indivs()->get_first();
  ae_individual* indiv      = NULL;
  int32_t nb=0;
  while( indiv_node != NULL )
  {
    indiv = (ae_individual *) indiv_node->get_obj();
    int *probes = new int32_t[5];
    probes[0]=nb;
    indiv->set_int_probes(probes);
    nb++;
    indiv_node = indiv_node->get_next();
  }

  exp_manager->save();
  delete exp_manager;

  return EXIT_SUCCESS;
}

void print_help( char* prog_name ) 
{
  printf( "\n\
Usage : addintprobe -h\n\
or :    addintprobe -r num_generation \n\
\t-h : display this screen\n\
\t-r num_generation  : modify num_generation to add a unique int probe to each individual, \n\
");
}
