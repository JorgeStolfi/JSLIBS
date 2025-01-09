#include <kdtree.h>
#include <float.h>
#include <r2.h>
#include <normais.h>
#include <assert.h>

kdnode_t* createKdNode(int level,int axis,int table_entry, double radius ){
	kdnode_t* novo =  (kdnode_t*)malloc(sizeof(kdnode_t));
	novo->child[0] = NULL;
	novo->child[1] = NULL;
	novo->level = level;
	novo->table_entry = table_entry;
	novo->radius = radius;
	return novo;
}

void splitEntryList(
		int* entries,
		int num_entries,
		Tabela* tab,
		int axis,
		int* best_ref,
		double* radius,
		int** lo_entries,
		int* num_lo_entries,
		int** hi_entries,
		int* num_hi_entries
		);

void splitEntryList(
		int* entries,
		int num_entries,
		Tabela* tab,
		int axis,
		int* best_ref,
		double* radius,
		int** lo_entries,
		int* num_lo_entries,
		int** hi_entries,
		int* num_hi_entries
		){
	auto int compara_entries (const void * a, const void * b);
	int compara_entries  (const void * a, const void * b){
		int* p_a = (int*) a;
		int ref_a = *p_a;
		int* p_b = (int*) b;
		int ref_b = *p_b;
		const double* sig_a = get_intdir(tab,ref_a);
		const double* sig_b = get_intdir(tab,ref_b);
		if(sig_a[axis] > sig_b[axis]) return 1;
		if(sig_a[axis] < sig_b[axis]) return -1;
		return 0;
	}

	qsort(entries, num_entries, sizeof(int), compara_entries);
	
	int ind_best = num_entries/2;
	*best_ref = entries[ind_best];

	/*compute the radius*/
	double ro = 0.0;
	int i;
	const double* sig_best = get_intdir(tab,*best_ref);
	int num_dimensions = get_num_luzes(tab);
	for(i = 0; i < num_entries; i++){
		if( i != (*best_ref) ){
			const double* sig = get_intdir(tab,entries[i]);
			double dist = dist_euclid(sig,sig_best,num_dimensions);
			if(dist > ro){
				ro = dist;
			}
		}
		
	}

	*radius = ro;
	/*split ref_list in two tables*/	
	*num_lo_entries = num_entries/2;
	*num_hi_entries = num_entries - (*num_lo_entries) - 1;
	*lo_entries = entries;
	*hi_entries = entries + (*num_lo_entries + 1);
	if((*num_lo_entries) <= 0 ){
		*lo_entries = NULL;
	}
	if((*num_hi_entries) <= 0 ){
		*hi_entries = NULL;
	}
	
	return;
}

kdnode_t* buildKdSubTree(int* entries, int num_entries, Tabela* tab,int level,int num_axes,kdtree_stats_t* tree_stats){
	if(entries == NULL){
		return NULL;
	}
	/*if we have a ref table, we must find teh limits and the mean*/	
	
	int axis = level%num_axes;
	if(level > tree_stats->depth){
		tree_stats->depth = level;
	}
	tree_stats->level_population[level]++;
	
	/*fill data of node*/

	
	/*dirty trick - alloc worst case ref list to each one*/
	int* lo_entries;
	int* hi_entries;
	int best_ref;
	double radius;
	int num_lo_entries = 0;
	int num_hi_entries = 0;
	splitEntryList(
		entries,
		num_entries,
		tab,
		axis,
		&best_ref,
		&radius,
		&lo_entries,
		&num_lo_entries,
		&hi_entries,
		&num_hi_entries
		);
	
	kdnode_t* node = createKdNode(level,axis,best_ref,radius);

	/*invoques next call to nodes*/

	node->child[0] = buildKdSubTree(lo_entries, num_lo_entries, tab,level+1,num_axes,tree_stats);
	node->child[1] = buildKdSubTree(hi_entries, num_hi_entries, tab,level+1,num_axes,tree_stats);

	/*return node*/
	return node;
}


kdtree_t* buildKdTree(Tabela* tab){
	kdtree_t* novo = (kdtree_t*)malloc(sizeof(kdtree_t));
	novo->root = NULL;
	novo->tab = tab;
	/*prepare list of pointers*/
	int num_lines = get_num_linhas(tab);
	int num_axes = get_num_luzes(tab);
	int* entries = (int*)malloc(sizeof(int)*num_lines);
	int i;
	for(i = 0; i < num_lines; i++){
		entries[i] = i;
	}
	novo->stats.depth = -1;
	novo->stats.num_dist_euclids = 0.0;
	novo->stats.num_searches = 0;
	/*alloc array  for the worst case scenario*/
	novo->stats.level_population = (int*)malloc(sizeof(int)*num_lines);
	for(i = 0; i < num_lines; i++){
		novo->stats.level_population[i] = 0;
	}
	novo->root = buildKdSubTree(entries,num_lines, tab,0,num_axes,&novo->stats);
	int* aux_level_population = (int*)malloc(sizeof(int)*(novo->stats.depth));
	for(i = 0; i < novo->stats.depth; i++){
		aux_level_population[i] = novo->stats.level_population[i];
	}
	free(novo->stats.level_population);
	novo->stats.level_population = aux_level_population;
	free(entries);
	return novo;
}

void searchKdSubTree(kdnode_t* node,Tabela* tab,r2_t box[], double* so,int* best_entry,double* best_dist,double* dist_count){
	if(node == NULL) return;

	int num_dimensions = get_num_luzes(tab);
	int axis = node->level%num_dimensions;
	//procura distancia do ponto à caixa
	double dist_box = dist_point_box(so,box,num_dimensions);
	if(dist_box >= *best_dist){
		return;
	}
	*dist_count = (*dist_count) + (1.0/((double)num_dimensions));

	const double* go = get_intdir(tab,node->table_entry);
	double dist = dist_euclid(so,go,num_dimensions);
	*dist_count = (*dist_count) + 1.0;
	if(dist < *best_dist){
		*best_dist = dist;
		*best_entry = node->table_entry;
	}
	if(*best_dist < (dist - node->radius)){
		return ;
	}
	//divide the box in two parts and recurse
	r2_t save_box;
	if(so[axis] < go[axis] ){
		save_box = box[axis];
		box[axis].c[1] = go[axis];
		searchKdSubTree(node->child[0],tab,box,so,best_entry,best_dist,dist_count);
		box[axis] = save_box;
		box[axis].c[0] = go[axis];
		searchKdSubTree(node->child[1],tab,box,so,best_entry,best_dist,dist_count);
		box[axis] = save_box;
	}else{
		save_box = box[axis];
		box[axis].c[0] = go[axis];
		searchKdSubTree(node->child[1],tab,box,so,best_entry,best_dist,dist_count);
		box[axis] = save_box;
		box[axis].c[1] = go[axis];
		searchKdSubTree(node->child[0],tab,box,so,best_entry,best_dist,dist_count);
		box[axis] = save_box;
	}
}

int localiza_normal_kdtree(kdtree_t* tree, Tabela* tab, double SO[], double* dist, double *albedo,double* n_euclid_evalsP){
	int table_entry;
	int num_luzes = get_num_luzes(tab);
	double so[num_luzes];
  	double Smag;
  	extrai_assinatura(SO,so,&Smag,num_luzes);
	r2_t box[num_luzes];
	int i;
	for(i = 0; i < num_luzes; i++){
		box[i].c[0] = 0.0;
		box[i].c[1] = 1.0;
	}
	*dist = DBL_MAX;
	table_entry = -1;
	double dist_count = 0.0;
	searchKdSubTree(tree->root,tab,box,so,&table_entry,dist,&dist_count);
	tree->stats.num_dist_euclids = tree->stats.num_dist_euclids +  dist_count;
	tree->stats.num_searches++;
	if(n_euclid_evalsP != NULL) *n_euclid_evalsP = dist_count;
	if ((*dist) > 2.0) { (*dist) = 2.0; }
	assert((table_entry != -1));
  	(*albedo) = calcula_albedo(tab, table_entry, SO);
	return table_entry;
}

void printfTreeStats(FILE* arq, kdtree_t* tree){
	int num_lines = get_num_linhas(tree->tab);
	fprintf(arq,"Depth %d\n",tree->stats.depth);
	fprintf(arq,"Level Population\n");
	int i;
	for(i = 0; i < tree->stats.depth; i++){
		int pop = tree->stats.level_population[i];
		double percent = pop/(double)num_lines;
		fprintf(arq,"[%d]:%d - %2.3f%%\n",i,pop,percent*100.0);
	}
	fprintf(arq,"Searches Performed %d\n",tree->stats.num_searches);
	if ( tree->stats.num_dist_euclids <= 0.0) puts("KKKKKKKKKKKKKKKKKKKKKKKKKKKk");
	fprintf(arq,"Processing Time Total %lf\n",tree->stats.usec_time);
	fprintf(arq,"Processing Time Avg. %4.6lf\n",tree->stats.usec_time/(double)tree->stats.num_searches);
	fprintf(arq,"Euclidean Dist. Total %lf\n",tree->stats.num_dist_euclids);
	fprintf(arq,"Euclidean Dist. Avg %4.6lf\n",tree->stats.num_dist_euclids/(double)tree->stats.num_searches);
}

