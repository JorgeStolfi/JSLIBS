#include <super_tabela.h>
#include <jsmath.h>
#include <assert.h>


/*This library contains Super Table routines that are used to speed bucked grid search relisient to shadows and highligt effects
The main ideia is use alpha distance with bucket grids with a limited number of illuminations (suficient to recover information) to ensure
a good matching. The structure contains the number of original test subject as the number of used lighs, this respective table and buckd grid and
the necessary structs to identify the light sources used
*/

/*Thanks to Bruno Azevedo to give this idea whith the checkboard problem !*/

struct SuperTabela{
  int num_luzes; //this is the total numbers used in the test subject
  int subsetSize; //this one is the number of lights used in this super table
  Tabela* tab; //the contained table
  bucketGrid* bg; //our marvelous bucket grid
  bool_t* subset; //bool vector that marks with true lights used, it have m size
  int* ind_luz; //int vector with light indices (between 0...(m-1)) with size subsetSize
  int* tableIndex; //vector which contains the correspondent index in the Main Table
  int qtd_index; //Number of rows in tableIndex
};


/*Interface Functions of SuperTabela*/


/*first output functions*/

Tabela* stGetTable(SuperTabela* st){
	return st->tab;
}

bucketGrid* stGetBucketGrid(SuperTabela* st){
	return st->bg;
}

int stGetK(SuperTabela* st){
	return st->subsetSize;
}

int stGetNumLuzes(SuperTabela* st){
	return st->num_luzes;
}

const int* stGetIndLuz(SuperTabela* st){
	return st->ind_luz;
}

/* For now we dont need any of input functions !*/

/*Creation and Destruction of ST strtucture*/

SuperTabela* stCreate(
    int num_luzes, 
    int canal, 
    int subsetSize,
    int ind_luz[],
    int gridsize,
    Tabela* mainTable
){	
	/*ensure good conditions to generate table*/
	assert((subsetSize <= num_luzes) && (subsetSize >=3));
	/*alloc structure and init simple properties*/
	SuperTabela* st = (SuperTabela*) malloc(sizeof(SuperTabela));
	st->subsetSize = subsetSize;
	st->num_luzes = num_luzes;
	
	//float_image_t * BlankImg = float_image_new(3, nx, ny);
	int i;
	
	/*init indices in ST*/
	st->ind_luz = (int*)malloc(sizeof(int)*subsetSize);
	fprintf(stderr,"Subset:");
	for(i = 0; i < subsetSize;i++){
		st->ind_luz[i] = ind_luz[i];
		fprintf(stderr,"%03d ",st->ind_luz[i]);
	}
	fprintf(stderr,"\n");
	/*init bool vector*/
	st->subset = (bool_t*)malloc(sizeof(bool_t)*num_luzes);
	for(i = 0; i < num_luzes; i++){
		st->subset[i] = 0;
	}
	for(i = 0; i < subsetSize; i++){
		st->subset[ind_luz[i]] = 1;
	}
	st->tab = criaSubTabela(mainTable,subsetSize,ind_luz, &(st->tableIndex));
	/*and our trustfull bucket grid*/
	st->bg = CriaBucketGrid(st->tab,gridsize);
	st->qtd_index = get_num_linhas(st->tab);
	
	int count_no_match = 0;
	fprintf(stderr,"\n%d Non Matched Normals of %d\n",count_no_match,st->qtd_index);
  	
	return st;
}


void stFree(SuperTabela** st){
	/* We DONT HAVE ANY FREE FUNCTIONS FOR BG AND TABLE !!!*/
	free((*st)->subset);
	free((*st)->ind_luz);
	free(*st);
}

/*TODO*/



int stSearchNormal(
    SuperTabela* st,
    const double SO[], 
    double *SO_mag,
    double *dist, 
    double *albedo,
    int  *matchedIndex,
    int* n_euclid_evalsP,
    int* n_scansP
){
	/*first we must prepare SO to use in bucket grid*/
        double SO_st[st->subsetSize];
	//double SO_st2[st->num_luzes];
	int i;
	double sum2 = 0;
	for(i = 0; i < st->subsetSize; i++){
		SO_st[i] = SO[st->ind_luz[i]];
		sum2+= SO_st[i]*SO_st[i];
		//printf("%d ",st->ind_luz[i]);
	}
	//printf("\n");
	/*
  	for(i = 0; i < st->num_luzes; i++){
		if(st->subset) SO_st2[i] = SO[i];
		else SO_st2[i] = 0;
		sum2+= SO_st2[i]*SO_st2[i];
	}*/

	*SO_mag = sqrt(sum2);
	//return localiza_normal_hash(st->bg, st->tab, SO_st, dist, albedo);
	int line  = localiza_normal_hash(st->bg, st->tab, SO_st, dist, albedo,n_euclid_evalsP,n_scansP);
	if(line != -1) *matchedIndex = st->tableIndex[line];
	return line;
}
int compareInt(const void* a, const void* b);
int compareInt(const void* a, const void* b){
	int va,vb;
	va = *((int*)a);
	vb = *((int*)b);
	return va - vb;
}

void stGenerateSubsets(int** subsets,int subset_size,int num_luzes,int num_subsets){
	int i = 0; //next candidate.
	int step = 1; //step between candidates.
	int t = 0; //number of candidates at current step.
	
	int k; //subset index.
	for(k = 0; k < num_subsets; k++){
		int r; //element index.
		for( r = 0; r< subset_size; r++){
			subsets[k][r] = i;
			t++;
			if(t >= num_luzes){
				//get the next step.
				do{
					step= (step +1)% num_luzes;
				}while(gcd(step,num_luzes) != 1);
				assert((step != 1)); //if it fails there too many sets to generate
				t = 0;
			}
			i=(i+step)%num_luzes;
		}
	}
	
}

void stGenerateSubsets_deprecated(int** subsets,int subsetSize,int m,int M){
	int i,j;
	int count_M,count_K;
	int selected_matrix[M][subsetSize];
	int sorted_set[subsetSize];
	int selected_vector[m];
	count_M = 0;
	while(count_M < M){
		count_K = 0;
		for(i=0; i < M;i++){
			selected_vector[i] = 0;
		}
		while(count_K < subsetSize){
			int sorteado = rand()%m;
			if(selected_vector[sorteado] == 0){
				selected_vector[sorteado] =1;
				sorted_set[count_K] = sorteado;
				count_K++;
				
			}
		}
		qsort(sorted_set,subsetSize,sizeof(int),compareInt);
		int test = 1;
		for(i = 0; i < count_M;i++){
			int test_vector = 1;
			for(j = 0; j < subsetSize; j++){
				if(sorted_set[j] != selected_matrix[i][j]){
					test_vector = 0;
					break;
				}
			}
			if(test_vector == 1){
				test = 0;
				break;
			}
		}
		if(test == 1){
			for(i = 0; i < subsetSize;i++){
				selected_matrix[count_M][i] = sorted_set[i];
			}
			count_M++;
		}
	}
	for(i = 0; i < M; i++){
		for(j = 0; j < subsetSize;j++){
			subsets[i][j] = selected_matrix[i][j];
		}
	}
	
}

void LiberaSuperTabela( SuperTabela* st){
 
  LiberaTabela(st->tab);
  LiberaBucketGrid(st->bg);
  free(st->subset);
  free(st->ind_luz);
  free(st->tableIndex);
  free(st);
 
};
