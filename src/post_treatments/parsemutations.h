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

ae_replication_report***  reports;
int32_t** reproductive_success;

int32_t get_nb_descendant_per_generation(int32_t generation, int32_t individual, int32_t nbgen){
  ae_replication_report* report = reports[generation][individual];
  return 0;
}

int32_t get_total_nb_descendant(int32_t generation, int32_t individual, ae_replication_report* rep){
  ae_replication_report* report = reports[generation][individual];
  return 0;
}