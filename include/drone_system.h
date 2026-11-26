#ifndef DRONE_SYSTEM_H
#define DRONE_SYSTEM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#define MAX_NODES 100
#define MAX_DRONES 50
#define MAX_PARCELS 500
#define MAX_STRING 256

typedef struct {
    int drone_id;
    double max_battery;
    double current_battery;
    double payload_capacity;
    int current_node;
    bool is_active;
} Drone;

typedef struct {
    int parcel_id;
    int source_node;
    int dest_node;
    double weight;
    int priority;
    int deadline_minutes;
    bool delivered;
} Parcel;

typedef struct {
    int node_id;
    char name[MAX_STRING];
    double x, y;
    bool is_no_fly_zone;
} Node;

typedef struct {
    int u, v;
    double weight;
    int capacity;
    bool blocked;
} Edge;

typedef struct {
    Node nodes[MAX_NODES];
    int node_count;
    Edge edges[MAX_NODES * MAX_NODES];
    int edge_count;
    Drone drones[MAX_DRONES];
    int drone_count;
    Parcel parcels[MAX_PARCELS];
    int parcel_count;
} WorldState;

void init_world(WorldState *world);
void load_sample_world(WorldState *world);
void run_menu(WorldState *world);
void simulate_normal_operation(WorldState *world);
void simulate_constraint_shock(WorldState *world);
void simulate_stress_test(WorldState *world);
void run_custom_live_test(void);
void compare_payload_strategies_and_export(WorldState *world);
#endif