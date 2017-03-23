/*
    Tests process in parallel.
    -----
    Return values:
    0      result is correct
    1      result is wrong
    2      problem size does not match
    253    no "data_output" file
    254    no "data_input" file
    -----
    Reference:
    http://mpi.deino.net/mpi_functions/MPI_Allgather.html
*/
#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab4_IO.h"
#include "mpi.h"
#include "timer.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

#define THRESHOLD 0.0001

int main (int argc, char* argv[]){
    struct node *nodehead;
    int nodecount;
    int *num_in_links, *num_out_links;
    double *r, *r_pre, *r_local;
    int i, j;
    double damp_const;
    int iterationcount = 0;
    int collected_nodecount;
    double *collected_r;
    double cst_addapted_threshold;
    double error;
    double start, end;
    int chunksize;

    int rank, numProcs;

    FILE *fp;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // Adjust the threshold according to the problem size
    cst_addapted_threshold = THRESHOLD;
    
    // Calculate the result
    if (get_node_stat(&nodecount, &num_in_links, &num_out_links)) return 254;
    if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;
    
    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    r_local = malloc(nodecount / numProcs * sizeof(double));

    chunksize = nodecount / numProcs;
    int begin_row = rank * chunksize;
    int end_row = (rank+1) * chunksize;

    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;
    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;

    // CORE CALCULATION 
    GET_TIME(start);
    // do broadcast r here
    do { 
        ++iterationcount; 
        vec_cp(r, r_pre, nodecount);
        // use allgather
        
        for ( i = begin_row; i < end_row; ++i) {
            r_local[i] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j) {
                r_local[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            }
            r_local[i] *= DAMPING_FACTOR;
            r_local[i] += damp_const;
        }
        MPI_Allgather(r_local, chunksize, MPI_DOUBLE, r, chunksize, MPI_DOUBLE, MPI_COMM_WORLD);
    } while(rel_error(r, r_pre, nodecount) >= EPSILON);
    GET_TIME(end);

    // post processing
    Lab4_saveoutput(r, nodecount, end - start);
    node_destroy(nodehead, nodecount);
    free(num_in_links); free(num_out_links);
    free(r); free(r_local); free(r_pre);
    MPI_Finalize();

}
