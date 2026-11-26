#include "../include/drone_system.h"
int main(int argc, char *argv[]) {
    WorldState world;
    init_world(&world);
    load_sample_world(&world);

    if (argc > 1 && strcmp(argv[1], "--benchmark") == 0) {
        printf("[Benchmark Runner] Executing automated performance suite...\n");
        simulate_normal_operation(&world);
        simulate_constraint_shock(&world);
        simulate_stress_test(&world);
        compare_payload_strategies_and_export(&world);
        return 0;
    }

    run_menu(&world);
    return 0;
}