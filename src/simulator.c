#include "../include/drone_system.h"
#include "../include/algorithms.h"

#ifdef _WIN32
#include <windows.h>
static double now_us(void) {
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart * 1e6 / (double)freq.QuadPart;
}
#else
static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec / 1e3;
}
#endif

void init_world(WorldState *world) {
    world->node_count = 0;
    world->edge_count = 0;
    world->drone_count = 0;
    world->parcel_count = 0;
}

void load_sample_world(WorldState *world) {
    init_world(world);

    world->nodes[0] = (Node){0, "Central Hub", 0.0, 0.0, false};
    world->nodes[1] = (Node){1, "Sector A Delivery", 10.0, 15.0, false};
    world->nodes[2] = (Node){2, "Sector B Relay", 20.0, 5.0, false};
    world->nodes[3] = (Node){3, "Sector C Outpost", 30.0, 25.0, false};
    world->node_count = 4;

    // Edges with weights AND corridor capacities (units/hr)
    world->edges[0] = (Edge){0, 1, 18.0, 12, false};
    world->edges[1] = (Edge){1, 2, 14.0, 8,  true};  // Blocked corridor
    world->edges[2] = (Edge){0, 2, 22.0, 10, false};
    world->edges[3] = (Edge){2, 3, 20.0, 9,  false};
    world->edges[4] = (Edge){1, 3, 25.0, 10, false};
    world->edge_count = 5;

    world->drones[0] = (Drone){101, 100.0, 100.0, 5.0, 0, true};
    world->drones[1] = (Drone){102, 150.0, 150.0, 8.0, 0, true};
    world->drone_count = 2;

    world->parcels[0] = (Parcel){1, 0, 1, 2.5, 9, 30, false};
    world->parcels[1] = (Parcel){2, 0, 3, 4.0, 7, 45, false};
    world->parcels[2] = (Parcel){3, 1, 3, 1.5, 5, 60, false};
    world->parcel_count = 3;

    printf("[System] Loaded Autonomous Drone Logistics World (%d nodes, %d edges, %d drones, %d parcels).\n",
           world->node_count, world->edge_count, world->drone_count, world->parcel_count);
}

void simulate_normal_operation(WorldState *world) {
    printf("\n=== SCENARIO 1: NORMAL OPERATION (Dijkstra + MergeSort + MST + TSP) ===\n");
    int dist[MAX_NODES], prev[MAX_NODES];
    dijkstra(world, 0, dist, prev);
    for (int i = 0; i < world->node_count; i++) {
        printf("Shortest path Hub -> Node %d (%s): %d units\n",
               i, world->nodes[i].name, dist[i]);
    }

    merge_sort_parcels(world->parcels, 0, world->parcel_count - 1);
    printf("Top-priority parcel after Merge Sort: ID=%d (Priority=%d, Weight=%.2f kg)\n",
           world->parcels[0].parcel_id, world->parcels[0].priority, world->parcels[0].weight);

    prim_mst(world);
    kruskal_mst(world);

    // Multi-stop tour via bitmask DP with double** memo
    int n = world->node_count;
    if (n >= 2 && n <= 12) {
        double cost[MAX_NODES][MAX_NODES];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                cost[i][j] = (i == j) ? 0.0 : 1e6;
        for (int e = 0; e < world->edge_count; e++) {
            if (!world->edges[e].blocked) {
                int u = world->edges[e].u, v = world->edges[e].v;
                cost[u][v] = cost[v][u] = world->edges[e].weight;
            }
        }
        int states = 1 << n;
        double **memo = (double **)malloc((size_t)states * sizeof(double *));
        for (int m = 0; m < states; m++) {
            memo[m] = (double *)malloc((size_t)n * sizeof(double));
            for (int j = 0; j < n; j++) memo[m][j] = -1.0;
        }
        double tour = tsp_dp(1, 0, n, cost, memo);
        printf("[TSP Bitmask DP] Minimum closed delivery tour cost across %d nodes: %.2f\n", n, tour);
        for (int m = 0; m < states; m++) free(memo[m]);
        free(memo);
    }
}

void simulate_constraint_shock(WorldState *world) {
    printf("\n=== SCENARIO 2: CONSTRAINT SHOCK (Dynamic World Flow & Bipartite Matching) ===\n");
    int n = world->node_count;
    int sink = n - 1;
    int cap_graph[MAX_NODES][MAX_NODES] = {{0}};

    // Build capacity network directly from active WorldState edges
    for (int e = 0; e < world->edge_count; e++) {
        if (!world->edges[e].blocked) {
            cap_graph[world->edges[e].u][world->edges[e].v] += world->edges[e].capacity;
        }
    }
    int normal_flow = ford_fulkerson(cap_graph, 0, sink, n);

    // Inject shock: corridor 0 -> 2 suffers weather/no-fly closure
    int shock_graph[MAX_NODES][MAX_NODES] = {{0}};
    for (int e = 0; e < world->edge_count; e++) {
        if (!world->edges[e].blocked && !(world->edges[e].u == 0 && world->edges[e].v == 2)) {
            shock_graph[world->edges[e].u][world->edges[e].v] += world->edges[e].capacity;
        }
    }
    int shock_flow = ford_fulkerson(shock_graph, 0, sink, n);

    printf("Baseline Max Throughput (Node 0 -> Node %d): %d units/hr\n", sink, normal_flow);
    printf("Post-Shock Max Throughput (Corridor 0->2 closed): %d units/hr\n", shock_flow);

    // Dynamic drone-to-parcel feasibility graph based on payload capacity
    int bp_graph[MAX_NODES][MAX_NODES] = {{0}};
    for (int d = 0; d < world->drone_count; d++) {
        for (int p = 0; p < world->parcel_count; p++) {
            if (world->drones[d].is_active && world->parcels[p].weight <= world->drones[d].payload_capacity) {
                bp_graph[d][p] = 1;
            }
        }
    }
    int match[MAX_NODES];
    int matched = max_bipartite_matching(bp_graph, world->drone_count, world->parcel_count, match);
    printf("Feasible Drone-to-Parcel Assignments (%d drones x %d parcels): %d matched\n",
           world->drone_count, world->parcel_count, matched);
}

void simulate_stress_test(WorldState *world) {
    (void)world;
    printf("\n=== SCENARIO 3: GENERATED LARGE-DATA STRESS CASE (N=250 Parcels, V=25 Nodes) ===\n");

    // 1. Generate 250 synthetic parcels deterministically
    int n_parcels = 250;
    double *wt = (double *)malloc((size_t)n_parcels * sizeof(double));
    double *val = (double *)malloc((size_t)n_parcels * sizeof(double));
    Parcel *plist = (Parcel *)malloc((size_t)n_parcels * sizeof(Parcel));
    if (!wt || !val || !plist) { free(wt); free(val); free(plist); return; }

    for (int i = 0; i < n_parcels; i++) {
        wt[i] = 0.50 + (double)((i * 17) % 450) / 100.0; // 0.50 to 4.99 kg
        val[i] = 10.0 + (double)((i * 31) % 90);
        plist[i] = (Parcel){i + 1, 0, (i % 24) + 1, wt[i], (int)val[i], 60, false};
    }

    long long dp_ops = 0, greedy_ops = 0;
    double t0 = now_us();
    double dp_best = knapsack_01_dp(25.0, wt, val, n_parcels, &dp_ops);
    double t1 = now_us();

    double t2 = now_us();
    double gr_best = fractional_knapsack(25.0, plist, n_parcels, &greedy_ops);
    double t3 = now_us();

    printf("Stress Payload (N=%d, Cap=25.0 kg):\n", n_parcels);
    printf("  0/1 DP Knapsack:    Value=%.2f | Ops=%lld | Time=%.2f us\n", dp_best, dp_ops, t1 - t0);
    printf("  Greedy Fractional:  Value=%.2f | Ops=%lld | Time=%.2f us\n", gr_best, greedy_ops, t3 - t2);

    // 2. Huffman tree on telemetry event distribution
    char syms[] = {'N', 'W', 'B', 'R', 'E', 'C'};
    int freqs[] = {450, 210, 140, 95, 65, 40};
    huffman_coding(syms, freqs, 6);

    // 3. Strassen Matrix Multiplication on 8x8 energy transition matrix
    double A[MAX_NODES][MAX_NODES] = {{0}}, B[MAX_NODES][MAX_NODES] = {{0}}, C[MAX_NODES][MAX_NODES] = {{0}};
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++) {
            A[i][j] = (double)((i + j) % 5 + 1);
            B[i][j] = (i == j) ? 2.0 : 0.5;
        }
    double ts0 = now_us();
    strassen_matrix_mult(8, A, B, C);
    double ts1 = now_us();
    printf("[Strassen D&C Matrix Mult] 8x8 fleet energy transition computed in %.2f us (C[0][0]=%.2f)\n",
           ts1 - ts0, C[0][0]);

    free(wt);
    free(val);
    free(plist);
}

void compare_payload_strategies_and_export(WorldState *world) {
    // Non-trivial instance where Greedy Fractional and 0/1 DP DIVERGE:
    // Cap = 7.0 kg
    // Item 1: 4.0 kg, val 100 (ratio 25.0)
    // Item 2: 3.5 kg, val 65  (ratio 18.57)
    // Item 3: 3.5 kg, val 65  (ratio 18.57)
    // Greedy takes Item 1 (4.0kg, 100) + 3.0/3.5 of Item 2 (55.71) = 155.71
    // 0/1 DP cannot split Item 2, so it takes Item 2 + Item 3 (7.0kg) = 130.00
    double capacity = 7.0;
    double wt[]  = {4.0, 3.5, 3.5};
    double val[] = {100.0, 65.0, 65.0};
    int n = 3;

    Parcel greedy_copy[3];
    for (int i = 0; i < n; i++) {
        greedy_copy[i].parcel_id = i + 1;
        greedy_copy[i].weight = wt[i];
        greedy_copy[i].priority = (int)val[i];
    }

    long long greedy_ops = 0, dp_ops = 0;

    double t0 = now_us();
    double greedy_score = fractional_knapsack(capacity, greedy_copy, n, &greedy_ops);
    double t1 = now_us();

    double t2 = now_us();
    double dp_score = knapsack_01_dp(capacity, wt, val, n, &dp_ops);
    double t3 = now_us();

    double greedy_us = t1 - t0;
    double dp_us = t3 - t2;

    printf("\n=== SAME-INPUT ALGORITHMIC COMPARISON (Payload Cap: %.2f kg) ===\n", capacity);
    printf("Input Items: [(4.0kg, val=100), (3.5kg, val=65), (3.5kg, val=65)]\n");
    printf("Fractional Knapsack (Greedy): Score=%.2f | Ops=%lld | Time=%.2f us | Space=O(n) sort\n",
           greedy_score, greedy_ops, greedy_us);
    printf("0/1 Knapsack (1D DP):         Score=%.2f | Ops=%lld | Time=%.2f us | Space=O(W_scaled)\n",
           dp_score, dp_ops, dp_us);
    printf("Analytical Takeaway: Greedy overestimates (155.71) by splitting a parcel; 0/1 DP finds the true indivisible optimum (130.00).\n");

    FILE *fp = fopen("mission_report.json", "w");
    if (fp != NULL) {
        fprintf(fp, "{\n  \"project\": \"Project 38 - Drone Delivery Fleet Mission Planner\",\n");
        fprintf(fp, "  \"nodes\": %d,\n  \"drones\": %d,\n  \"parcels\": %d,\n",
                world->node_count, world->drone_count, world->parcel_count);
        fprintf(fp, "  \"comparison\": {\n");
        fprintf(fp, "    \"capacity_kg\": %.2f,\n", capacity);
        fprintf(fp, "    \"greedy_fractional_value\": %.2f,\n    \"greedy_ops\": %lld,\n    \"greedy_time_us\": %.2f,\n",
                greedy_score, greedy_ops, greedy_us);
        fprintf(fp, "    \"dp_01_value\": %.2f,\n    \"dp_01_ops\": %lld,\n    \"dp_01_time_us\": %.2f\n",
                dp_score, dp_ops, dp_us);
        fprintf(fp, "  }\n}\n");
        fclose(fp);
        printf("[Export] Saved benchmark & comparison report to mission_report.json\n");
    }
}
void run_custom_live_test(void);
void run_menu(WorldState *world) {
    int choice;
    do {
        printf("\n=========================================\n");
        printf(" DRONE DELIVERY FLEET MISSION PLANNER (PBL-38)\n");
        printf("=========================================\n");
        printf("1. Load / Reset Sample Mission World\n");
        printf("2. Run Normal Operation (Dijkstra + MergeSort + Prim/Kruskal + TSP DP)\n");
        printf("3. Run Constraint Shock (Dynamic World Flow + Bipartite Matching)\n");
        printf("4. Run Generated Large-Data Stress Test (N=250 DP + Huffman + Strassen)\n");
        printf("5. Run String Matching Telemetry Verification (Naive / Rabin-Karp / KMP)\n");
        printf("6. Compare Payload Strategies (Divergent Case + Ops + us) & Export JSON\n");
        printf("7. Run Custom Live Input Test (Payload / String Search)\n");
        printf("8. Exit\n");
        printf("Enter your choice: ");
        if (scanf("%d", &choice) != 1) break;

        switch (choice) {
            case 1: load_sample_world(world); break;
            case 2: simulate_normal_operation(world); break;
            case 3: simulate_constraint_shock(world); break;
            case 4: simulate_stress_test(world); break;
            case 5: {
                char text[] = "DRONE_LOG_SUCCESS_TAKEOFF_ALTITUDE_50M";
                char pattern[] = "SUCCESS";
                printf("Running String Matching algorithms on telemetry logs:\n");
                naive_string_match(text, pattern);
                rabin_karp_match(text, pattern, 101);
                kmp_match(text, pattern);
                break;
            }
            case 6: compare_payload_strategies_and_export(world); break;
            case 7: run_custom_live_test(); break;
            case 8: printf("Exiting Drone Fleet Mission Planner.\n"); break;
            default: printf("Invalid choice. Try again.\n");
        }
    } while (choice != 7);
}
void run_custom_live_test(void) {
    int mode;
    printf("\n=== CUSTOM LIVE INPUT TESTER ===\n");
    printf("1. Custom Payload Comparison (Fractional Greedy vs 0/1 DP)\n");
    printf("2. Custom Telemetry String Search (Naive vs Rabin-Karp vs KMP)\n");
    printf("Select test mode: ");
    if (scanf("%d", &mode) != 1) return;

    if (mode == 1) {
        double cap;
        int n;
        printf("Enter drone payload capacity (kg): ");
        if (scanf("%lf", &cap) != 1) return;
        printf("Enter number of parcels (1-%d): ", MAX_PARCELS);
        if (scanf("%d", &n) != 1 || n <= 0 || n > MAX_PARCELS) return;

        double wt[MAX_PARCELS], val[MAX_PARCELS];
        Parcel items[MAX_PARCELS];
        for (int i = 0; i < n; i++) {
            printf("Parcel %d weight (kg) and priority value: ", i + 1);
            if (scanf("%lf %lf", &wt[i], &val[i]) != 2) return;
            items[i] = (Parcel){i + 1, 0, 1, wt[i], (int)val[i], 30, false};
        }

        long long g_ops = 0, d_ops = 0;
        double g_res = fractional_knapsack(cap, items, n, &g_ops);
        double d_res = knapsack_01_dp(cap, wt, val, n, &d_ops);
        printf("Result -> Greedy Fractional: %.2f (ops=%lld) | 0/1 DP: %.2f (ops=%lld)\n",
               g_res, g_ops, d_res, d_ops);
    } else if (mode == 2) {
        char text[MAX_STRING], pat[MAX_STRING];
        printf("Enter telemetry log text (no spaces): ");
        if (scanf("%255s", text) != 1) return;
        printf("Enter pattern to search: ");
        if (scanf("%255s", pat) != 1) return;
        naive_string_match(text, pat);
        rabin_karp_match(text, pat, 101);
        kmp_match(text, pat);
    }
}