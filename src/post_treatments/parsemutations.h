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

void print_help( void );
void loadreports( void );
void computereproductivesuccess( void );
void clean ( void );
double snapshot2gen( int32_t gen0, int32_t gen1, int32_t* results);
void computerelatedness( void );

inline int32_t gety(int32_t individual);
inline int32_t getx(int32_t individual);

ae_exp_manager* exp_manager;
char tree_file_name[50];
FILE* output_file;
FILE* relatedness_file;

int32_t nb_indivs;
int32_t nb_geners;

int32_t begin_gener;
int32_t end_gener;
int32_t stepgen;
int32_t ngen;
int32_t rwindow;

int32_t popx;
int32_t popy;

ae_replication_report***  reports;
int32_t** reproductive_success;
int32_t*** reproductive_success_bygen;
int32_t** bigger_reproductive_success;
int32_t** gen_bigger_reproductive_success;

