//
//  parsemutations.h
//  Aevol4.3local
//
//  Created by Antoine Frénoy on 02/04/2014.
//
//

#include <errno.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <sys/stat.h>


#include <ae_macros.h>
#include <ae_utils.h>
#ifndef __NO_X
#include <ae_exp_manager_X11.h>
#else
#include <ae_exp_manager.h>
#endif
#include <ae_individual.h>
#include <ae_genetic_unit.h>
#include <ae_list.h>
#include <ae_tree.h>
#include <ae_replication_report.h>
#include <ae_dna_replic_report.h>
#include <ae_mutation.h>

int32_t nb_indivs;
int32_t nb_geners;

ae_replication_report***  reports;
int32_t** reproductive_success;
int32_t*** reproductive_success_bygen;
int32_t** bigger_reproductive_success;
int32_t** bigger_targetgen;

/*
int32_t** snapshot2gen(int32_t gen0, int32_t gen1)
{
  int32_t
  
}
*/
