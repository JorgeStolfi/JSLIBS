#include <tabela.h>
#include <r3.h>
#include <math.h>
#include <assert.h>

Tabela* TrimTable(Tabela* tab, double minZ);

Tabela* TrimTable(Tabela* tab, double minZ){
	int num_linhas = get_num_linhas(tab);
	int i;
	int trim_linhas = 0;
	for(i = 0; i < num_linhas; i++){
		r3_t normal = get_normal(tab,i);
		if(normal.c[2] > minZ){
			trim_linhas++;
		}
	}

	Tabela* novaTab;
	int num_luzes = get_num_luzes(tab);
	r3_t viewDir = get_view_dir(tab);
	novaTab = aloca_tabela_vazia(num_luzes,trim_linhas,viewDir);
	int count = 0;
	for(i = 0; i < num_linhas; i++){
		r3_t normal = get_normal(tab,i);
		const double* go = get_intdir(tab,i);
		double Gmag = get_intmag(tab,i);
		if(normal.c[2] > minZ){
			set_normal(novaTab,count,normal);
			set_intmag(novaTab,count,Gmag);
			set_intdir(novaTab,count,(double*)go);
			count++;
		}
	}
	return novaTab;
}


int main(int argc, char** argv){
	if(argc != 4){
		fprintf(stderr,"Chamada do programa\nprograma <in tab> <out tab> <minZ>\n");
		return 1;
	}
	char* inFileName = argv[1];
	char* outFileName = argv[2];
	double minZ;
	sscanf(argv[3],"%lf",&minZ);
	
	
	
	Tabela* inTab ;
	LoadTable(inFileName,&inTab);
	if(inTab == 0){
	  fprintf(stderr,"Cannot load table %s\n",inFileName);
	  assert(FALSE);
	}
	Tabela* outTab = TrimTable(inTab,minZ);
	SaveTable(outFileName,outTab,TRUE);
	fprintf(stderr,"Input Table has %d, Output Table has %d\n",get_num_linhas(inTab),get_num_linhas(outTab));
	return 0;
}
