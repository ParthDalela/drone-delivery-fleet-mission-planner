#include "../include/algorithms.h"

// --- 1. DIVIDE AND CONQUER ---
static void merge(Parcel arr[], int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;
    Parcel *L = (Parcel *)malloc((size_t)n1 * sizeof(Parcel));
    Parcel *R = (Parcel *)malloc((size_t)n2 * sizeof(Parcel));
    if (!L || !R) { free(L); free(R); return; }

    for (int i = 0; i < n1; i++) L[i] = arr[l + i];
    for (int j = 0; j < n2; j++) R[j] = arr[m + 1 + j];

    int i = 0, j = 0, k = l;
    while (i < n1 && j < n2) {
        if (L[i].priority >= R[j].priority) arr[k++] = L[i++];
        else arr[k++] = R[j++];
    }
    while (i < n1) arr[k++] = L[i++];
    while (j < n2) arr[k++] = R[j++];
    free(L);
    free(R);
}

void merge_sort_parcels(Parcel arr[], int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;
        merge_sort_parcels(arr, l, m);
        merge_sort_parcels(arr, m + 1, r);
        merge(arr, l, m, r);
    }
}

static double **alloc_mat(int sz) {
    double **m = (double **)malloc((size_t)sz * sizeof(double *));
    if (!m) return NULL;
    for (int i = 0; i < sz; i++) {
        m[i] = (double *)calloc((size_t)sz, sizeof(double));
        if (!m[i]) {
            for (int j = 0; j < i; j++) free(m[j]);
            free(m);
            return NULL;
        }
    }
    return m;
}

static void free_mat(double **m, int sz) {
    if (!m) return;
    for (int i = 0; i < sz; i++) free(m[i]);
    free(m);
}

static void mat_add(int sz, double **A, double **B, double **R) {
    for (int i = 0; i < sz; i++)
        for (int j = 0; j < sz; j++)
            R[i][j] = A[i][j] + B[i][j];
}

static void mat_sub(int sz, double **A, double **B, double **R) {
    for (int i = 0; i < sz; i++)
        for (int j = 0; j < sz; j++)
            R[i][j] = A[i][j] - B[i][j];
}

static void strassen_rec(int sz, double **A, double **B, double **C) {
    if (sz == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }
    int k = sz / 2;
    double **A11 = alloc_mat(k), **A12 = alloc_mat(k), **A21 = alloc_mat(k), **A22 = alloc_mat(k);
    double **B11 = alloc_mat(k), **B12 = alloc_mat(k), **B21 = alloc_mat(k), **B22 = alloc_mat(k);
    double **M1 = alloc_mat(k), **M2 = alloc_mat(k), **M3 = alloc_mat(k), **M4 = alloc_mat(k);
    double **M5 = alloc_mat(k), **M6 = alloc_mat(k), **M7 = alloc_mat(k);
    double **T1 = alloc_mat(k), **T2 = alloc_mat(k);

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            A11[i][j] = A[i][j];         A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];     A22[i][j] = A[i + k][j + k];
            B11[i][j] = B[i][j];         B12[i][j] = B[i][j + k];
            B21[i][j] = B[i + k][j];     B22[i][j] = B[i + k][j + k];
        }
    }

    mat_add(k, A11, A22, T1); mat_add(k, B11, B22, T2); strassen_rec(k, T1, T2, M1);
    mat_add(k, A21, A22, T1); strassen_rec(k, T1, B11, M2);
    mat_sub(k, B12, B22, T2); strassen_rec(k, A11, T2, M3);
    mat_sub(k, B21, B11, T2); strassen_rec(k, A22, T2, M4);
    mat_add(k, A11, A12, T1); strassen_rec(k, T1, B22, M5);
    mat_sub(k, A21, A11, T1); mat_add(k, B11, B12, T2); strassen_rec(k, T1, T2, M6);
    mat_sub(k, A12, A22, T1); mat_add(k, B21, B22, T2); strassen_rec(k, T1, T2, M7);

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            C[i][j]         = M1[i][j] + M4[i][j] - M5[i][j] + M7[i][j];
            C[i][j + k]     = M3[i][j] + M5[i][j];
            C[i + k][j]     = M2[i][j] + M4[i][j];
            C[i + k][j + k] = M1[i][j] - M2[i][j] + M3[i][j] + M6[i][j];
        }
    }

    free_mat(A11, k); free_mat(A12, k); free_mat(A21, k); free_mat(A22, k);
    free_mat(B11, k); free_mat(B12, k); free_mat(B21, k); free_mat(B22, k);
    free_mat(M1, k);  free_mat(M2, k);  free_mat(M3, k);  free_mat(M4, k);
    free_mat(M5, k);  free_mat(M6, k);  free_mat(M7, k);
    free_mat(T1, k);  free_mat(T2, k);
}

void strassen_matrix_mult(int n, double A[MAX_NODES][MAX_NODES], double B[MAX_NODES][MAX_NODES], double C[MAX_NODES][MAX_NODES]) {
    if (n <= 0 || n > MAX_NODES) return;
    int m = 1;
    while (m < n) m <<= 1;

    double **pA = alloc_mat(m);
    double **pB = alloc_mat(m);
    double **pC = alloc_mat(m);
    if (!pA || !pB || !pC) { free_mat(pA, m); free_mat(pB, m); free_mat(pC, m); return; }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            pA[i][j] = A[i][j];
            pB[i][j] = B[i][j];
        }
    }

    strassen_rec(m, pA, pB, pC);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = pC[i][j];

    free_mat(pA, m);
    free_mat(pB, m);
    free_mat(pC, m);
}

// --- 2. GREEDY METHODS ---
static void merge_ratio(Parcel arr[], int l, int m, int r, long long *ops) {
    int n1 = m - l + 1, n2 = r - m;
    Parcel *L = (Parcel *)malloc((size_t)n1 * sizeof(Parcel));
    Parcel *R = (Parcel *)malloc((size_t)n2 * sizeof(Parcel));
    if (!L || !R) { free(L); free(R); return; }
    for (int i = 0; i < n1; i++) L[i] = arr[l + i];
    for (int j = 0; j < n2; j++) R[j] = arr[m + 1 + j];

    int i = 0, j = 0, k = l;
    while (i < n1 && j < n2) {
        if (ops) (*ops)++;
        double r1 = (L[i].weight > 0) ? ((double)L[i].priority / L[i].weight) : 0.0;
        double r2 = (R[j].weight > 0) ? ((double)R[j].priority / R[j].weight) : 0.0;
        if (r1 >= r2) arr[k++] = L[i++];
        else arr[k++] = R[j++];
    }
    while (i < n1) arr[k++] = L[i++];
    while (j < n2) arr[k++] = R[j++];
    free(L);
    free(R);
}

static void sort_parcels_by_ratio(Parcel arr[], int l, int r, long long *ops) {
    if (l < r) {
        int m = l + (r - l) / 2;
        sort_parcels_by_ratio(arr, l, m, ops);
        sort_parcels_by_ratio(arr, m + 1, r, ops);
        merge_ratio(arr, l, m, r, ops);
    }
}

double fractional_knapsack(double capacity, Parcel parcels[], int n, long long *ops) {
    if (ops) *ops = 0;
    if (n <= 0 || capacity <= 0.0) return 0.0;

    sort_parcels_by_ratio(parcels, 0, n - 1, ops);

    double total_value = 0.0;
    for (int i = 0; i < n; i++) {
        if (ops) (*ops)++;
        if (parcels[i].weight <= 0.0) continue;
        if (capacity >= parcels[i].weight) {
            capacity -= parcels[i].weight;
            total_value += parcels[i].priority;
        } else {
            total_value += parcels[i].priority * (capacity / parcels[i].weight);
            break;
        }
    }
    return total_value;
}

typedef struct HuffNode {
    char ch;
    int freq;
    struct HuffNode *left, *right;
} HuffNode;

static HuffNode *new_huff_node(char ch, int freq, HuffNode *l, HuffNode *r) {
    HuffNode *node = (HuffNode *)malloc(sizeof(HuffNode));
    if (!node) return NULL;
    node->ch = ch;
    node->freq = freq;
    node->left = l;
    node->right = r;
    return node;
}

static void min_heapify(HuffNode **heap, int sz, int idx) {
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;
    if (left < sz && heap[left]->freq < heap[smallest]->freq) smallest = left;
    if (right < sz && heap[right]->freq < heap[smallest]->freq) smallest = right;
    if (smallest != idx) {
        HuffNode *tmp = heap[idx];
        heap[idx] = heap[smallest];
        heap[smallest] = tmp;
        min_heapify(heap, sz, smallest);
    }
}

static HuffNode *extract_min(HuffNode **heap, int *sz) {
    HuffNode *top = heap[0];
    heap[0] = heap[--(*sz)];
    min_heapify(heap, *sz, 0);
    return top;
}

static void insert_heap(HuffNode **heap, int *sz, HuffNode *node) {
    int i = (*sz)++;
    while (i > 0 && node->freq < heap[(i - 1) / 2]->freq) {
        heap[i] = heap[(i - 1) / 2];
        i = (i - 1) / 2;
    }
    heap[i] = node;
}

static void print_huff_codes(HuffNode *root, char *code, int depth, int *total_bits) {
    if (!root) return;
    if (!root->left && !root->right) {
        code[depth] = '\0';
        int bits = (depth == 0) ? 1 : depth;
        if (depth == 0) strcpy(code, "0");
        *total_bits += root->freq * bits;
        printf("  Symbol '%c' (freq=%d) -> code: %s (%d bits)\n",
               root->ch, root->freq, code, bits);
        return;
    }
    code[depth] = '0';
    print_huff_codes(root->left, code, depth + 1, total_bits);
    code[depth] = '1';
    print_huff_codes(root->right, code, depth + 1, total_bits);
}

static void free_huff_tree(HuffNode *root) {
    if (!root) return;
    free_huff_tree(root->left);
    free_huff_tree(root->right);
    free(root);
}

void huffman_coding(char data[], int freq[], int size) {
    if (size <= 0 || !data || !freq) {
        printf("[Huffman Coding] Invalid input.\n");
        return;
    }
    HuffNode **heap = (HuffNode **)malloc((size_t)size * sizeof(HuffNode *));
    if (!heap) return;

    int heap_sz = 0;
    int total_symbols = 0;
    for (int i = 0; i < size; i++) {
        heap[heap_sz++] = new_huff_node(data[i], freq[i], NULL, NULL);
        total_symbols += freq[i];
    }
    for (int i = (heap_sz - 2) / 2; i >= 0; i--) min_heapify(heap, heap_sz, i);

    while (heap_sz > 1) {
        HuffNode *left = extract_min(heap, &heap_sz);
        HuffNode *right = extract_min(heap, &heap_sz);
        HuffNode *merged = new_huff_node('$', left->freq + right->freq, left, right);
        insert_heap(heap, &heap_sz, merged);
    }

    HuffNode *root = heap[0];
    char code[MAX_STRING];
    int compressed_bits = 0;
    int fixed_bits_per_sym = (size <= 1) ? 1 : (int)ceil(log2((double)size));
    int fixed_total = total_symbols * fixed_bits_per_sym;

    printf("[Huffman Tree] Built optimal prefix tree for %d symbols:\n", size);
    print_huff_codes(root, code, 0, &compressed_bits);
    printf("  Total bits: %d (Huffman) vs %d (Fixed %d-bit) -> %.1f%% savings\n",
           compressed_bits, fixed_total, fixed_bits_per_sym,
           fixed_total > 0 ? (100.0 * (fixed_total - compressed_bits) / fixed_total) : 0.0);

    free_huff_tree(root);
    free(heap);
}

void dijkstra(WorldState *world, int src, int dist[], int prev[]) {
    int n = world->node_count;
    bool sptSet[MAX_NODES];
    for (int i = 0; i < n; i++) { dist[i] = INT_MAX; sptSet[i] = false; prev[i] = -1; }
    dist[src] = 0;

    for (int count = 0; count < n - 1; count++) {
        int min = INT_MAX, u = -1;
        for (int v = 0; v < n; v++)
            if (!sptSet[v] && dist[v] <= min) { min = dist[v]; u = v; }
        if (u == -1) break;
        sptSet[u] = true;

        for (int e = 0; e < world->edge_count; e++) {
            if (world->edges[e].blocked) continue;
            int v = -1;
            if (world->edges[e].u == u) v = world->edges[e].v;
            else if (world->edges[e].v == u) v = world->edges[e].u;
            if (v != -1 && !sptSet[v] && !world->nodes[v].is_no_fly_zone) {
                int w = (int)llround(world->edges[e].weight);
                if (dist[u] != INT_MAX && dist[u] + w < dist[v]) {
                    dist[v] = dist[u] + w;
                    prev[v] = u;
                }
            }
        }
    }
}

void prim_mst(WorldState *world) {
    int V = world->node_count;
    if (V <= 0) return;
    double key[MAX_NODES];
    bool in_mst[MAX_NODES];
    int parent[MAX_NODES];

    for (int i = 0; i < V; i++) { key[i] = 1e18; in_mst[i] = false; parent[i] = -1; }
    key[0] = 0.0;

    double total_cost = 0.0;
    for (int count = 0; count < V; count++) {
        int u = -1;
        double min_val = 1e18;
        for (int v = 0; v < V; v++) {
            if (!in_mst[v] && key[v] < min_val) { min_val = key[v]; u = v; }
        }
        if (u == -1) break;
        in_mst[u] = true;
        if (parent[u] != -1) total_cost += key[u];

        for (int e = 0; e < world->edge_count; e++) {
            if (world->edges[e].blocked) continue;
            int v = -1;
            if (world->edges[e].u == u) v = world->edges[e].v;
            else if (world->edges[e].v == u) v = world->edges[e].u;
            if (v != -1 && !in_mst[v] && world->edges[e].weight < key[v]) {
                key[v] = world->edges[e].weight;
                parent[v] = u;
            }
        }
    }
    printf("[Prim's MST] Total active relay/hub network cost: %.2f\n", total_cost);
}

static int find_set(int parent[], int x) {
    if (parent[x] == x) return x;
    return parent[x] = find_set(parent, parent[x]);
}

void kruskal_mst(WorldState *world) {
    int V = world->node_count;
    if (V <= 1) return;
    int parent[MAX_NODES];
    for (int i = 0; i < V; i++) parent[i] = i;

    int active_count = 0;
    Edge sorted_edges[MAX_NODES * MAX_NODES];
    for (int i = 0; i < world->edge_count; i++) {
        if (!world->edges[i].blocked) sorted_edges[active_count++] = world->edges[i];
    }
    for (int i = 0; i < active_count - 1; i++) {
        for (int j = i + 1; j < active_count; j++) {
            if (sorted_edges[j].weight < sorted_edges[i].weight) {
                Edge tmp = sorted_edges[i]; sorted_edges[i] = sorted_edges[j]; sorted_edges[j] = tmp;
            }
        }
    }

    double mst_cost = 0.0;
    int edges_used = 0;
    for (int i = 0; i < active_count && edges_used < V - 1; i++) {
        int ru = find_set(parent, sorted_edges[i].u);
        int rv = find_set(parent, sorted_edges[i].v);
        if (ru != rv) {
            parent[ru] = rv;
            mst_cost += sorted_edges[i].weight;
            edges_used++;
        }
    }
    printf("[Kruskal's MST] Edges in spanning tree: %d/%d | Total cost: %.2f\n",
           edges_used, V - 1, mst_cost);
}

// --- 3. DYNAMIC PROGRAMMING ---
double knapsack_01_dp(double W, double wt[], double val[], int n, long long *ops) {
    if (ops) *ops = 0;
    int scale = 100;
    int cap = (int)llround(W * scale);
    if (cap <= 0 || n <= 0) return 0.0;

    double *dp = (double *)calloc((size_t)(cap + 1), sizeof(double));
    if (!dp) {
        fprintf(stderr, "[Error] Memory allocation failed in knapsack_01_dp.\n");
        return -1.0;
    }

    for (int i = 0; i < n; i++) {
        int item_w = (int)llround(wt[i] * scale);
        if (item_w < 0 || item_w > cap) continue;
        for (int w = cap; w >= item_w; w--) {
            if (ops) (*ops)++;
            double candidate = val[i] + dp[w - item_w];
            if (candidate > dp[w]) dp[w] = candidate;
        }
    }

    double res = dp[cap];
    free(dp);
    return res;
}

double tsp_dp(int mask, int pos, int n, double cost[MAX_NODES][MAX_NODES], double **memo) {
    if (mask == (1 << n) - 1) return cost[pos][0];
    if (memo[mask][pos] >= 0.0) return memo[mask][pos];
    double ans = 1e18;
    for (int city = 0; city < n; city++) {
        if ((mask & (1 << city)) == 0) {
            double newAns = cost[pos][city] + tsp_dp(mask | (1 << city), city, n, cost, memo);
            if (newAns < ans) ans = newAns;
        }
    }
    return memo[mask][pos] = ans;
}

// --- 4. FLOW AND MATCHING ---
static bool bfs_flow(int rGraph[MAX_NODES][MAX_NODES], int s, int t, int parent[], int n) {
    bool visited[MAX_NODES] = {false};
    int queue[MAX_NODES], head = 0, tail = 0;
    queue[tail++] = s; visited[s] = true; parent[s] = -1;
    while (head < tail) {
        int u = queue[head++];
        for (int v = 0; v < n; v++) {
            if (!visited[v] && rGraph[u][v] > 0) {
                queue[tail++] = v; parent[v] = u; visited[v] = true;
            }
        }
    }
    return visited[t];
}

int ford_fulkerson(int graph[MAX_NODES][MAX_NODES], int s, int t, int n) {
    int rGraph[MAX_NODES][MAX_NODES];
    for (int u = 0; u < n; u++)
        for (int v = 0; v < n; v++)
            rGraph[u][v] = graph[u][v];
    int parent[MAX_NODES];
    int max_flow = 0;
    while (bfs_flow(rGraph, s, t, parent, n)) {
        int path_flow = INT_MAX;
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            if (rGraph[u][v] < path_flow) path_flow = rGraph[u][v];
        }
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
        max_flow += path_flow;
    }
    return max_flow;
}

static bool bpm_dfs(int bp_graph[MAX_NODES][MAX_NODES], int u, bool seen[], int match[], int n) {
    for (int v = 0; v < n; v++) {
        if (bp_graph[u][v] && !seen[v]) {
            seen[v] = true;
            if (match[v] < 0 || bpm_dfs(bp_graph, match[v], seen, match, n)) {
                match[v] = u;
                return true;
            }
        }
    }
    return false;
}

int max_bipartite_matching(int bp_graph[MAX_NODES][MAX_NODES], int m, int n, int match[]) {
    int count = 0;
    for (int i = 0; i < n; i++) match[i] = -1;
    for (int u = 0; u < m; u++) {
        bool seen[MAX_NODES] = {false};
        if (bpm_dfs(bp_graph, u, seen, match, n)) count++;
    }
    return count;
}

// --- 5. STRING MATCHING ---
void naive_string_match(char *text, char *pattern) {
    int n = (int)strlen(text);
    int m = (int)strlen(pattern);
    for (int i = 0; i <= n - m; i++) {
        int j;
        for (j = 0; j < m; j++) if (text[i + j] != pattern[j]) break;
        if (j == m) printf("[Naive Match] Found pattern at index %d\n", i);
    }
}

void rabin_karp_match(char *text, char *pattern, int q) {
    int d = 256;
    int n = (int)strlen(text);
    int m = (int)strlen(pattern);
    int p = 0, t = 0, h = 1;
    for (int i = 0; i < m - 1; i++) h = (h * d) % q;
    for (int i = 0; i < m; i++) {
        p = (d * p + pattern[i]) % q;
        t = (d * t + text[i]) % q;
    }
    for (int i = 0; i <= n - m; i++) {
        if (p == t) {
            int j;
            for (j = 0; j < m; j++) if (text[i + j] != pattern[j]) break;
            if (j == m) printf("[Rabin-Karp] Found pattern at index %d\n", i);
        }
        if (i < n - m) {
            t = (d * (t - text[i] * h) + text[i + m]) % q;
            if (t < 0) t += q;
        }
    }
}

void kmp_match(char *text, char *pattern) {
    int m = (int)strlen(pattern);
    int n = (int)strlen(text);
    int lps[MAX_STRING];
    int len = 0; lps[0] = 0;
    int i = 1;
    while (i < m) {
        if (pattern[i] == pattern[len]) { len++; lps[i] = len; i++; }
        else if (len != 0) len = lps[len - 1];
        else { lps[i] = 0; i++; }
    }
    i = 0; int j = 0;
    while (i < n) {
        if (pattern[j] == text[i]) { i++; j++; }
        if (j == m) { printf("[KMP Match] Found pattern at index %d\n", i - j); j = lps[j - 1]; }
        else if (i < n && pattern[j] != text[i]) {
            if (j != 0) j = lps[j - 1];
            else i++;
        }
    }
}