/* See pst_signature_cube.h */
/* Last edited on 2006-11-03 18:51:07 by stolfi */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <jspnm.h>
#include <pst_signature.h>
#include <pst_signature_cube.h>

/* 
ltn_cube_t *ltn_cubify_table(int size)
  {
    ltn_cube_t *cb = alloc_cube(size);
    // ..
    return cb;
  }

void verifica_cubo(Cubo* c);
struct Lista{

	//float intensidade[3][3]; //primeiro indice - arquivo, segundo indice - canal
	long int linha;
	struct Lista* prox;
};

struct Cube{
	int resolucao;
	lista****  cubo;
};

struct Vet{
	float x;
	float y;
	float z;
};

typedef struct Vet vet;

lista* insere(lista* l,long int linha){
	int ind,ind2;
	lista* novo = (lista*)malloc(sizeof(lista));
	novo->linha = linha;
	// for(ind =0;ind<3;ind++){
	// 	for(ind2 =0; ind2 < 3 ; ind2++){
	// 		novo->intensidade[ind][ind2] = i[ind][ind2];
	// 	}
	// 	// o primeiro índice refere-se ao arquivo!!!  o segundo ao canal !!!//
	// }
	novo->prox = l;
	return novo;
}

void func_U(pixel cA,pixel cB, pixel cC, vet* novo, int canal){
	// float corA = Calcula_Intesidade((float) cA.canal[canal]);
	// float corB = Calcula_Intesidade((float) cB.canal[canal]);
	// float corC = Calcula_Intesidade((float) cC.canal[canal]);
	float corA = (float)cA.canal[canal];
	float corB = (float)cB.canal[canal];
	float corC = (float)cC.canal[canal];

	novo->x = corA/(sqrt((corA*corA)+(corB*corB)+(corC*corC)));
	novo->y = corB/(sqrt((corA*corA)+(corB*corB)+(corC*corC)));
	novo->z = corC/(sqrt((corA*corA)+(corB*corB)+(corC*corC)));
	return;
}

void func_V(Tabela* tab,long int linha,vet* novo, int canal){
	float vetA = get_intensidade(tab,linha,0,canal);
	float vetB = get_intensidade(tab,linha,1,canal);
	float vetC = get_intensidade(tab,linha,2,canal);

	novo->x = vetA/(sqrt((vetA*vetA)+(vetB*vetB)+(vetC*vetC)));
	novo->y = vetB/(sqrt((vetA*vetA)+(vetB*vetB)+(vetC*vetC)));
	novo->z = vetC/(sqrt((vetA*vetA)+(vetB*vetB)+(vetC*vetC)));
//	printf("F - ");
	//imprime_vetor(*novo);
	return;
}

float distancia(vet u, vet v){
	float X = u.x - v.x;
	float Y = u.y - v.y;
	float Z = u.z - v.z;
	return sqrt((X*X)+(Y*Y)+(Z*Z));
}


Cubo* gera_cubo(int resolucao){
	lista**** nova_lista;
	int i,j,k;
	int tamanho;
	Cubo* novo = (Cubo*)malloc(sizeof(Cubo));
	novo->resolucao = resolucao;
	nova_lista = (lista****)malloc(sizeof(lista***)*resolucao);
	for(i = 0;i<resolucao;i++){
		nova_lista[i]= (lista***)malloc(sizeof(lista**)*resolucao);
		for(j = 0;j<resolucao;j++){
			nova_lista[i][j] =(lista**)malloc(sizeof(lista*)*resolucao);
			for(k=0;k<resolucao;k++){
				nova_lista[i][j][k] = NULL;
			}
		}
	}
	novo->cubo = nova_lista;
	return novo;
}


float Media_ponderada(float R, float G, float B){
	float temp;
	temp = wRed*R;
	temp = temp + (wGrn*G);
	temp = temp + (wBlu*B);
	//temp = fabs(temp);// fabs - funcao do math.h que retorna o valor absoluto de um  float
	return temp;
}

int index_cubo(float num,int r){
	return   num*r;
}

void processa_tabela(Tabela* tab, Cubo* cubo){
	int ind,n_arq,canal;
	int a,b,c;
	float** intensidade;//[arquivo] [canal];
	float dx;
	float dy;
	float dz;
	lista** lp;
	lista* l;
	//printf("chegou");
	intensidade = (float**)malloc(sizeof(float*)*3);
	for(ind = 0;ind<3;ind++){
		intensidade[ind] = (float*)malloc(sizeof(float)*3);
	}

	for(ind = 0;ind<get_tamanho_tabela(tab);ind++){
		for(n_arq =0;n_arq<3;n_arq++){
			for(canal = 0;canal<3;canal++){
				intensidade[n_arq][canal] = get_intensidade(tab,ind,n_arq, canal);
				//printf("[%d][%d]: %f\n",n_arq,canal,intensidade[n_arq][canal]);
			}
		}
		//calculando o índice da tabela
		//a = (int) floor(media_ponderada(intensidade[0][0],intensidade[0][1],intensidade[0][2])*(cubo->resolucao));
		//b = (int) floor(media_ponderada(intensidade[1][0],intensidade[1][1],intensidade[1][2])*(cubo->resolucao));
		//c = (int) floor(media_ponderada(intensidade[2][0],intensidade[2][1],intensidade[2][2])*(cubo->resolucao));
		//dx = Media_ponderada(intensidade[0][0],intensidade[0][1],intensidade[0][2]);
		//printf("Media:%f Index:%d \n",dx, index_cubo(dx,6));
		a = index_cubo(Media_ponderada(intensidade[0][0],intensidade[0][1],intensidade[0][2]),cubo->resolucao);//x
		b = index_cubo(Media_ponderada(intensidade[1][0],intensidade[1][1],intensidade[1][2]),cubo->resolucao);
		c = index_cubo(Media_ponderada(intensidade[2][0],intensidade[2][1],intensidade[2][2]),cubo->resolucao);
		//printf("Tab [%d]: A:%d B:%d C:%d X:%f Y%f Z:%f \n",ind,a,b,c,vet[0],vet[1],vet[2]);
		//printf("A:%d B:%d C:%d\n",a,b,c);
		lp = &(cubo->cubo[a][b][c]);
		*lp = insere(*lp,ind);
		// if(cubo->cubo[a][b][c] != NULL){
		// 	puts("OK");
		// }
	}
        //verifica_cubo(cubo);

}

int conta_lista(lista* l){
        int cont;
        lista* aux;
        cont =0;
        for(aux =l;aux!=NULL;aux = aux->prox){
                cont++;
        }
        return cont;
}

void verifica_cubo(Cubo* c){
        int inda,indb,indc;
	lista* l;
        long int cont =0;
        for(inda =0;inda< c->resolucao;inda++){
                for(indb =0;indb< c->resolucao;indb++){
                        for(indc =0;indc < c->resolucao;indc++){

                                cont = conta_lista(c->cubo[inda][indb][indc]) + cont;
                                if(conta_lista(c->cubo[inda][indb][indc]) != 0){
                                        printf("Cubo[%d][%d][%d]: %d\n",inda,indb,indc,conta_lista(c->cubo[inda][indb][indc]));
				        printf("-----------------------------------------------");
				        for(l = c->cubo[inda][indb][indc];l!=NULL;l = l->prox){
        					printf("%d -",l->linha);
	        			}
		        		printf("\n");
                                }

                        }
                }
        }

	printf("TOTAL:%d\n",cont);
}


int obtem_indice(int k,int ind,int tam){
	int aux;
	aux = k +ind;
	if(aux < 0){
		//printf("(MInus)");
		aux = 0;
	}
	if(aux > (tam-1)){
		//printf("(MAIS) ");
		aux = tam-1;
	}
	return aux;
}

void media_vetor(vet* dest,vet* source, int canal){
	vet d = *dest;
	vet s = *source;
	if(canal == 0){
		d.x = d.x + (wRed*s.x);
		d.y = d.y + (wRed*s.y);
		d.z = d.z + (wRed*s.z);

	}
	else if(canal == 1){
		d.x = d.x + (wGrn*s.x);
		d.y = d.y + (wGrn*s.y);
		d.z = d.z + (wGrn*s.z);

	}
	else if(canal == 2){
		d.x = d.x + (wBlu*s.x);
		d.y = d.y + (wBlu*s.y);
		d.z = d.z + (wBlu*s.z);

	}
	else{
		d.x = 0;
		d.y = 0;
		d.z = 0;

	}
}

long int busca_cubo(Tabela* tab,pixel* cores,int canal,Cubo* cubo){
	long int resposta;
	int area_varredura = 3;
	int cont,inda,indb,indc;
	float ia,ib,ic;
	int ka,kb,kc;
	int ind,inicializado;
	float menor;
	float menor_dist;
	lista* ind_lista;
	vet* u;
	vet* v;
	u = (vet*) malloc(sizeof(vet));
	v = (vet*) malloc(sizeof(vet));
	//func_U(cores[0],cores[1],cores[2],u,canal);
	ka = index_cubo(cores[0].canal[canal]/255.0,cubo->resolucao);
        kb = index_cubo(cores[1].canal[canal]/255.0,cubo->resolucao);
	kc = index_cubo(cores[2].canal[canal]/255.0,cubo->resolucao);
        //printf("KA:%d KB:%d KC:%d\n",ka,kb,kc);
	inicializado = 0;
	resposta = -1;
	ind =0;
	for(inda = -2; inda <= 2 ; inda++){
		for(indb = - 2;indb<= 2; indb++){
			for(indc = -2;indc <= 2;indc++){
				int auxa = obtem_indice(ka,inda,cubo->resolucao);
				int auxb = obtem_indice(kb,indb,cubo->resolucao);
				int auxc = obtem_indice(kc,indc,cubo->resolucao);
				//printf("A:%d B:%d C:%d\n",auxa,auxb,auxc);
				for(ind_lista = cubo->cubo[auxa][auxb][auxc];ind_lista!=NULL;ind_lista = ind_lista->prox){

					long int linha_atual = ind_lista->linha;
					func_U(cores[0],cores[1],cores[2],u,canal);
					func_V(tab,linha_atual,v,canal);
			//		printf("X:%f Y:%f Z:%f \n",v->x,v->y,v->z);
					if(inicializado == 0){
						menor = distancia(*u,*v);
						resposta = linha_atual;
						inicializado = 1;
					}
					else{
						if(distancia(*u,*v) < menor){
							menor = distancia(*u,*v);
							resposta = linha_atual;
						}
					}
					ind++;
				}
			}
		}
	}
	//printf("RESP: %d \n",resposta);
        if(resposta == -1){
          //      printf(" N� entrou");
        }
	//printf("---------------\n");
    //    printf("\n");
	return resposta;
}

*/
