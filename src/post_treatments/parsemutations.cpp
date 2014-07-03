//
//  parsemutations.cpp
//  Aevol4.3local
//
//  Created by Antoine Frénoy on 02/04/2014.
//
//

#include "parsemutations.h"

/*
 
 This post-treatment analyses all the reproductive events (being caracterised by a set of mutations) and outputs the reproductive success of each of these reproductive events.
 
 ngen: Number of generations to take into account to calculate reproductive success. Should be a few times higher than population radius.

 stepgen: Discretisation of the algorithm: we analyse successively (for all i) from i*stepgen to (i+1)*stepgen, so using data from i*stepgen to (i+1)*stepgen + ngen.
 Does not change the output but changes time and space complexity. If low, we use more CPU but less memory (because there are more independent and smaller steps)
 
*/


int main(int argc, char** argv)
{
  // Parse the parameters
  begin_gener = 0;
  end_gener = -1;
  ngen = 1000;
  stepgen = 500;
  
  const char * short_options = "hb:e:g:f:";
  static struct option long_options[] = {
    {"help",        no_argument,       NULL,  'h'},
    {"begin",       required_argument, NULL,  'b'},
    {"end",         required_argument, NULL,  'e' },
    {"granularity", required_argument, NULL,  'g' },
    {"future",      required_argument, NULL,  'f' },
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
      case 'g' : stepgen = atol(optarg); break;
      case 'f' : ngen = atol(optarg); break;
    }
  }
  
  // Load the simulation
  loadreports();
  
  // Compute reproductive success of the individuals
  computereproductivesuccess();
  
  // Close the files and clean memory
  clean();
  
  exit(EXIT_SUCCESS);
  
}

void loadreports( void )
{
  if ( end_gener == -1 )
  {
    printf( "error: You must provide a generation number.\n");
    exit( EXIT_FAILURE );
  }

#ifndef __NO_X
  exp_manager = new ae_exp_manager_X11();
#else
  exp_manager = new ae_exp_manager();
#endif
  exp_manager->load( end_gener, false, true, false );
  
  int32_t tree_step = exp_manager->get_tree_step();
  nb_indivs = exp_manager->get_nb_indivs();
  nb_geners = end_gener - begin_gener;
  
  // Load the tree files and copy replication reports into big table reports
  ae_tree * tree = NULL;
  
  reports = new ae_replication_report**[nb_geners];
  int32_t i;
  for(i=0;i<end_gener - begin_gener;i++){
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
}

void computereproductivesuccess( void )
{
  // Open the output file
  char output_file_name[101];
  snprintf( output_file_name, 100, "mutations-b%06"PRId32"-e%06"PRId32".txt", begin_gener, end_gener);
  
  output_file = fopen(output_file_name, "w");
  if ( output_file == NULL )
  {
    fprintf(stderr, "Failed to create file %s, exiting.\n", output_file_name);
    exit(EXIT_FAILURE);
  }

  // Write header to output file
  fprintf(output_file,"#Syntax:\n #generation\n #individual\n #total reproductive success\n #highest reproductive success among all generations\n #generation at which highest reproductive success is reached\n #metabolic effect of the mutation set\n #secretion effect of the mutation set\n\n");

  // Allocate memory
  reproductive_success_bygen=new int32_t**[stepgen+ngen];
  bigger_reproductive_success=new int32_t*[stepgen];
  gen_bigger_reproductive_success=new int32_t*[stepgen];
  reproductive_success=new int32_t*[stepgen];
  
  int32_t i,j;
  for(i=0;i<stepgen+ngen;i++){
    reproductive_success_bygen[i] = new int32_t*[nb_indivs];
    for(j=0;j<nb_indivs;j++){
      reproductive_success_bygen[i][j] = new int32_t[stepgen+ngen];
    }
  }
  
  for (i=0;i<stepgen;i++){
    bigger_reproductive_success[i]=new int32_t[nb_indivs];
    gen_bigger_reproductive_success[i]=new int32_t[nb_indivs];
    reproductive_success[i]=new int32_t[nb_indivs];
  }
  
  // For each stepgen generations
  int32_t nbstep=(nb_geners-ngen)/stepgen; // We do not analyse the last ngen generations to avoid a creating a bias
  int32_t step;
  for (step=0;step<nbstep;step++){
    int32_t gena=step*stepgen; // First generation we analyse
    // int32_t genb=(step+1)*stepgen; // Last generation we analyse
    // int32_t genc=genb+ngen; // Last generation we need to consider to be able to analyse until genb
    
    
    // Initialize memory with 0s
    int32_t generation,individual,relgeneration,reltargetgen;
    
    for (relgeneration=0;relgeneration<stepgen+ngen;relgeneration++) {
      for (individual=0;individual<nb_indivs;individual++){
        for (reltargetgen=0;reltargetgen<stepgen+ngen;reltargetgen++){
          reproductive_success_bygen[relgeneration][individual][reltargetgen]=0;
        }
      }
    }
    
    for (relgeneration=0;relgeneration<stepgen;relgeneration++) {
      for (individual=0;individual<nb_indivs;individual++){
        bigger_reproductive_success[relgeneration][individual]=0;
        gen_bigger_reproductive_success[relgeneration][individual]=0;
        reproductive_success[relgeneration][individual]=0;
      }
    }
    
    // Dynamic programming algorithm
    
    for (relgeneration=stepgen+ngen-1;relgeneration>0;relgeneration--) {
      generation=relgeneration+gena;
      for (individual=0;individual<nb_indivs;individual++){
        // Find parent and update its reproductive success
        int32_t idparent = reports[generation][individual]->get_parent_id(); // the report that describes creation of individual at generation + 1 from parent_id at generation.
        reproductive_success_bygen[relgeneration-1][idparent][relgeneration]+=1;
        for (reltargetgen=relgeneration+1;reltargetgen<stepgen+ngen;reltargetgen++){
          reproductive_success_bygen[relgeneration-1][idparent][reltargetgen]+=reproductive_success_bygen[relgeneration][individual][reltargetgen];
        }
      }
    }
    
    // Assert consistency: for all relgeneration, for all reltargetgen, the sum of reproductive success of all individuals should be equal to the number of individual
    for (relgeneration=0;relgeneration<stepgen+ngen;relgeneration++) {
      for (reltargetgen=relgeneration+1;reltargetgen<stepgen+ngen;reltargetgen++){
        int32_t sum_indivs=0;
        for (individual=0;individual<nb_indivs;individual++){
          sum_indivs+=reproductive_success_bygen[relgeneration][individual][reltargetgen];
        }
        assert(sum_indivs==nb_indivs);
      }
    }
    
    // Find total reproductive success and best reproductive success (among all generations) during ngen
    for (relgeneration=0;relgeneration<stepgen;relgeneration++) {
      for (individual=0;individual<nb_indivs;individual++){
        int32_t valuebest=0;
        int32_t indexbest=-1;
        int32_t totsuccess=0; // Start with zero because we do not count the focal individual in its offsprings
        for (reltargetgen=relgeneration+1;reltargetgen<=relgeneration+ngen;reltargetgen++){
          if (reproductive_success_bygen[relgeneration][individual][reltargetgen]>valuebest){
            valuebest=reproductive_success_bygen[relgeneration][individual][reltargetgen];
            indexbest=reltargetgen;
          }
          totsuccess+=reproductive_success_bygen[relgeneration][individual][reltargetgen];
        }
        bigger_reproductive_success[relgeneration][individual]=valuebest;
        gen_bigger_reproductive_success[relgeneration][individual]=indexbest;
        reproductive_success[relgeneration][individual]=totsuccess;
      }
    }
    
    // Assert consistency: for all generations, for all individuals, the sum of total reproductive success is nb_indivs*ngen
    for (relgeneration=0;relgeneration<stepgen;relgeneration++) {
      int32_t sum=0;
      for (individual=0;individual<nb_indivs;individual++){
        sum+=reproductive_success[relgeneration][individual];
      }
      assert(sum==nb_indivs*ngen);
    }
    
    // output all the analyzed mutations
    for (relgeneration=0;relgeneration<stepgen;relgeneration++){
      for (individual=0;individual<nb_indivs;individual++){
        double metabolic_effect=reports[gena+relgeneration][individual]->get_metabolic_error() - reports[gena+relgeneration][individual]->get_parent_metabolic_error();
        double secretion_effect=reports[gena+relgeneration][individual]->get_secretion_error() - reports[gena+relgeneration][individual]->get_parent_secretion_error();
        // negative value = smaller gap = beneficial mutation
        fprintf(output_file,"%"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %+.15f %+.15f\n", gena+relgeneration+begin_gener, individual, reproductive_success[relgeneration][individual], bigger_reproductive_success[relgeneration][individual], (gen_bigger_reproductive_success[relgeneration][individual]!=-1)?(gen_bigger_reproductive_success[relgeneration][individual]+gena+begin_gener):-1, metabolic_effect, secretion_effect);
      }
    }
    
    // We do not clean memory now because we will use it in next iteration
    
  }
}

void clean( void )
{
  fclose(output_file);
  
  int32_t i,j;
  for(i=0;i<stepgen+ngen;i++){
    for(j=0;j<nb_indivs;j++){
      delete [] reproductive_success_bygen[i][j];
    }
    delete [] reproductive_success_bygen[i];
  }
  delete [] reproductive_success_bygen;
  
  for (i=0;i<stepgen;i++){
    delete [] bigger_reproductive_success[i];
    delete [] gen_bigger_reproductive_success[i];
    delete [] reproductive_success[i];
  }
  delete [] bigger_reproductive_success;
  delete [] gen_bigger_reproductive_success;
  delete [] reproductive_success;
  
  for(i=0;i<end_gener - begin_gener;i++)
  {
    delete [] reports[i];
  }
  delete [] reports;
  
  delete exp_manager;
}

void print_help( void )
{
  printf( "Analysis and lineage of mutations post-treatment \n" );
  printf( "Usage : parsemutations -b gen0 -e gen1 -g gra -f fut [-h]\n");
  printf( "\n" );
  printf( "\t-b gen0 or --begin gen0 : \n" );
  printf( "\t                  First generatin to analyze\n" );
  printf( "\n" );
  printf( "\t-e gen1 or --end gen1 : \n" );
  printf( "\t                  Last generation to analyze \n" );
  printf( "\n" );
  printf( "\t-g stepgen or --granularity stepgen : \n" );
  printf( "\t                  How many generations to analyze at a time \n" );
  printf( "\n" );
  printf( "\t-f ngen or --futur ngen : \n" );
  printf( "\t                  How many future generations to consider when calculating reproductive success \n" );
  printf( "\n" );
}
