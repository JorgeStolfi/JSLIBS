#include <tabela.h>

#define KDTREE_VERSION_0 0
#define KDTREE_VERSION_1 1
#define KDTREE_VERSION_2 2

struct KdNode_t{
	int level; //level of node on the tree
	int table_entry; //
	double radius;
	struct KdNode_t* child[2]; //0 - low, 1 - high

};

typedef struct KdNode_t kdnode_t;

struct kdTree_Stats_t{
	int depth;
	int* level_population;
	int num_searches;
	double num_dist_euclids;
	double usec_time;
};

typedef struct kdTree_Stats_t kdtree_stats_t;

struct KdTree_t{
	kdnode_t* root;
	Tabela*  tab;
	kdtree_stats_t stats;
};

typedef struct KdTree_t kdtree_t;

kdnode_t* createKdNode(int level,int axis,int table_entry, double radius );
kdnode_t* buildKdSubTree(int* entries, int num_entries, Tabela* tab,int level,int num_axes,kdtree_stats_t* tree_stats);
kdtree_t* buildKdTree(Tabela* tab);
void searchKdSubTree(kdnode_t* node,Tabela* tab,r2_t box[], double* so,int* best_entry,double* best_dist,double* dist_count);
int localiza_normal_kdtree(kdtree_t* tree, Tabela* tab, double SO[], double* dist, double *albedo,double* n_euclid_evalsP);
void printfTreeStats(FILE* arq, kdtree_t* tree);
