#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "drone_system.h"

// 1. Divide and Conquer
void merge_sort_parcels(Parcel arr[], int l, int r);
void strassen_matrix_mult(int n, double A[MAX_NODES][MAX_NODES], double B[MAX_NODES][MAX_NODES], double C[MAX_NODES][MAX_NODES]);

// 2. Greedy Methods
double fractional_knapsack(double capacity, Parcel parcels[], int n, long long *ops);
void huffman_coding(char data[], int freq[], int size);
void dijkstra(WorldState *world, int src, int dist[], int prev[]);
void prim_mst(WorldState *world);
void kruskal_mst(WorldState *world);

// 3. Dynamic Programming
double knapsack_01_dp(double W, double wt[], double val[], int n, long long *ops);
double tsp_dp(int mask, int pos, int n, double cost[MAX_NODES][MAX_NODES], double **memo);

// 4. Flow and Matching
int ford_fulkerson(int graph[MAX_NODES][MAX_NODES], int s, int t, int n);
int max_bipartite_matching(int bp_graph[MAX_NODES][MAX_NODES], int m, int n, int match[]);

// 5. String Matching
void naive_string_match(char *text, char *pattern);
void rabin_karp_match(char *text, char *pattern, int q);
void kmp_match(char *text, char *pattern);

#endif