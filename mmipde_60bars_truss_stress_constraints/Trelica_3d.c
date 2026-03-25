    #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mex.h>
#include <conio.h>
/*#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>*/


#define D 3
#define Dkl 3


//void trelica_3D(double * u, double * sumConst, int m, int n, int *index, int * in, int ** conectividade, int ** vector, int numConst, int numNodes, int nodesConst, int nodesLoad, double *area, double * elasticity, double * length, double ** kr, int * nodesC, double * nodesL, double * fr, double ** C, double ** coefXY, double * F);
void trelica_3D(double * u, double * sumConst, int n, int * in, int ** conectividade, int ** vector, int numConst, int numNodes, double *area, double * elasticity, double * length, double ** kr, double * nodesC, double * nodesL, double * fr, double ** C, double ** coefXY, double * F, double * comp, double * nodesCC, double ro,double * aalfa, int tamNodesCC, double *Ke, double * Kg, double ** kg, double ** mr, double *Mmat);
void reduced_matrix (int numNodes, int n, double *A, double *elasticity, double *length, int **conectividade, double **C, double **kr, int numConst, int *in, int **vector, double **c, double *fr, double * nodesC, double * Ke, double ro, double * aalfa, double ** mr, double *Mmat);
double LU (int numNodes, int numConst, double ** kr, double * fr, double * u, int * in, double * nodesC);
void tension_stress (int n, int numNodes, double * u, double ** C, double * elasticity, double * length, int ** vector, double * F, double * sumConst, double * area, int **conectividade, double * aalfa);
void calculations (int n, int ** conectividade, int ** vector, int numNodes, int numConst, double * nodesL, int * in, double * fr, double ** C, double * length, double * area, int * nosSaem, double * nodesC);
//double strstate(double ** kr, double * u, double *fr, int dim);

//void compliance(double **kr, double *u,double * nodesC, int numNodes, double * comp, double *fr);
void ccrit(int n, double * sumConst, double * length, int ** vector, int numNodes, int numConst, double * Kg, double ** kg, int * in, double * area, double * elasticity, int ** conectividade, double ** C, double * aalfa);


/**
 * Truss function
 */
void mexFunction(int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[])
{
	int i, j, n, m = D * 2, numColConect, tamNodesCC;
	int * in, ** vector, numConst, numNodes, ** conectividade;
	double * sumConst, * u, * Ke, *Kg, *Mmat;
	double * area, * elasticity, * length, ** kr, * nodesL, * fr, ** C, ** coefXY, * F, * nodesC, * nodesCC, ro, * aalfa, ** kg, **mr;


	double * CVetor;
	double * conectVetor;

	double * comp;


    //mexPrintf("Hello World1!\n");

	n = (int)mxGetScalar(prhs[0]);
	area = mxGetPr(prhs[1]);
  numNodes = (int)mxGetScalar(prhs[2]);
  numColConect = (int)mxGetScalar(prhs[3]);
  conectVetor = mxGetPr(prhs[4]);
  conectividade = (int**) mxMalloc( n * sizeof(int*));
    //mexPrintf("alocou algo");
  for(i=0; i<n; i++) {
      conectividade[i] = (int*) mxMalloc( numColConect * sizeof(int));
      for(j=0; j<numColConect; j++) {
          conectividade[i][j] = (int)conectVetor[i*numColConect+j];
          //mexPrintf("nodesCood[%d][%d]=%d\n", i, j, conectividade[i][j]);
      }
  }
   // mexPrintf("OK1");
    CVetor = mxGetPr(prhs[5]);
    C = (double**) mxMalloc( numNodes * sizeof(double*));
    for(i=0; i<numNodes; i++) {
        C[i] = (double*) mxMalloc( 3 * sizeof(double));
        for(j=0; j<3; j++) {
            C[i][j] = CVetor[i*3+j];
        }
    }
    nodesC = mxGetPr(prhs[6]);

    nodesL = mxGetPr(prhs[7]);
    numConst = (int)mxGetScalar(prhs[8]);
    elasticity = mxGetPr(prhs[9]);
    nodesCC = mxGetPr(prhs[10]);
    aalfa = mxGetPr(prhs[11]);
    ro = (double)mxGetScalar(prhs[12]);
    tamNodesCC=(int)mxGetScalar(prhs[13]);

    /*for (i=0;i<numNodes*Dkl;i++){
    	mexPrintf("%d ",nodesC[i]);
	}
	mexPrintf("\n\n");*/



    //saida
    plhs[0] = mxCreateDoubleMatrix(D * numNodes, 1, mxREAL);
	u = (double*)mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
	sumConst = (double*)mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(numNodes*Dkl-numConst,numNodes*Dkl-numConst, mxREAL);
	Ke = (double*)mxGetPr(plhs[2]);
	plhs[3] = mxCreateDoubleMatrix(numNodes*Dkl-numConst,numNodes*Dkl-numConst, mxREAL);
	Kg = (double*)mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(numNodes*Dkl-numConst,numNodes*Dkl-numConst, mxREAL);
    Mmat = (double*)mxGetPr(plhs[4]);

	kr = (double**) mxMalloc((numNodes * Dkl - numConst + 1) * sizeof(double*));
	kg = (double**) mxMalloc((numNodes * Dkl - numConst + 1) * sizeof(double*));
	mr = (double**) mxMalloc((numNodes * Dkl - numConst + 1) * sizeof(double*));

	for(i = 0; i < (numNodes*Dkl-numConst+1); i++) {
		kr[i] =  (double*)mxMalloc((numNodes*Dkl-numConst+1) * sizeof(double));
		kg[i] =  (double*)mxMalloc((numNodes*Dkl-numConst+1) * sizeof(double));
		mr[i] =  (double*)mxMalloc((numNodes*Dkl-numConst+1) * sizeof(double));
	}

	fr =(double*) mxMalloc((numNodes*Dkl-numConst+1) * sizeof(double));


  in = (int*) mxMalloc(Dkl * numNodes * sizeof(int));
	coefXY = (double**) mxMalloc(n * sizeof(double*)); // coeficientes cx, cy, ...
	for (i = 0;i < n; i++) {
		coefXY[i] = (double*) mxMalloc(D * sizeof(double));
	}
	F = (double*) mxMalloc(n * sizeof(double));
	vector = (int**) mxMalloc(n * sizeof(int*)); // Monta o vetor de orientacao global
	for(i = 0; i < n; i++) {
    	vector[i] = (int*) mxMalloc(m * sizeof(int));
	}
	length = (double*) mxMalloc(n * sizeof(double));

	//mexPrintf ("%d %d %d %d\n",numNodes, Dkl, numConst, numNodes * Dkl - numConst + 1);





	//mexPrintf("Tentando executar");

	trelica_3D(u, sumConst, n, in, conectividade, vector, numConst, numNodes, area, elasticity, length, kr, nodesC, nodesL, fr, C, coefXY, F, comp, nodesCC, ro, aalfa,tamNodesCC, Ke, Kg, kg,mr, Mmat);



    //mexPrintf("Finalizou e vai desalocar");
    //mxFree (area);

	mxFree (length);

	//mxFree (nodesL);
	mxFree (in);
	mxFree (F);
	mxFree (fr);


	for (i = 0; i < n; i++) {
		mxFree (vector[i]);
		mxFree (coefXY[i]);
		mxFree (conectividade[i]);
		//mxFree (elasticity[i]);
	}
	//mxFree (elasticity);
	mxFree (vector);
	mxFree (coefXY);
	mxFree (conectividade);

	for (i = 0; i < numNodes; i++)
		mxFree (C[i]);

	mxFree (C);

	for (i = 0; i < (numNodes*Dkl-numConst+1); i++) {
		mxFree (kr[i]);
		mxFree (mr[i]);
        mxFree (kg[i]);

	}
	mxFree (kr);
	mxFree (mr);
    mxFree (kg);
    
    //mxFree(nodesC);


	/*for (i=0; i<numNodes*Dkl;i++)
		mxFree(nodesC[i]);
	mxFree (nodesC);*/

	//mxFree (aalfa);

	//mexPrintf("Fim");
}

void trelica_3D(double * u, double * sumConst, int n, int * in, int ** conectividade, int ** vector, int numConst, int numNodes, double *area, double * elasticity, double * length, double ** kr, double * nodesC, double * nodesL, double * fr, double ** C, double ** coefXY, double * F, double * comp, double * nodesCC, double ro, double * aalfa, int tamNodesCC,double *Ke, double * Kg, double ** kg, double ** mr, double *Mmat) {
    int i,j,flagNo=0,cont=0;
    //double state = 0.0;
    //double stateTol = 2.0;
    double det=0.0;
    int nosSaem[numNodes];
    /*double nodesC[numNodes*Dkl];

    for (i=0;i<numNodes*Dkl;i++){
    	nodesC[i]=nodesC[i];
	}*/

     /*mexPrintf("NodesC antes do calculations\n");
    for (i=0;i<numNodes*Dkl;i++){       
		mexPrintf("%f ",nodesC[i]);
	}
	mexPrintf("\n\n");*/
    
    
	calculations(n, conectividade, vector, numNodes, numConst, nodesL, in, fr, C, length, area, nosSaem, nodesC);

    /*mexPrintf("NodesC depois do calculations \n");
    for (i=0;i<numNodes*Dkl;i++){        
		mexPrintf("%f ",nodesC[i]);
	}
	mexPrintf("\n\n");*/

	numConst=0;
	for (i=0;i<Dkl*numNodes;i++){
		if (nodesC[i]==1){
			numConst+=1;
		}
	}

	if (nosSaem[0]!=0) {
		//mexPrintf("Entrou aqui\n");
		for (i=0;i<numNodes;i++){
			if (nosSaem[i]!=0){
				for (j=0;j<tamNodesCC;j++)
				if (nosSaem[i] == nodesCC[j]) {
					//mexPrintf("Nï¿½ de CC que foi tirado: %d \n",nosSaem[i]);
					flagNo=1;
					break;
				}
			}
		}
	}
    
    
    //mexPrintf("tamanho do vetor nodescc = %d %d \n",sizeof((int)nodesCC),sizeof((int)nodesCC)/sizeof((int)nodesCC[0]));


	//numConst=numConst-cont;

	//mexPrintf("numConst %d\nFlagNo %d \n",numConst,flagNo);

	if (flagNo!=1){

			reduced_matrix(numNodes, n, area, elasticity, length, conectividade, C, kr, numConst, in, vector, coefXY, fr, nodesC, Ke, ro, aalfa,mr, Mmat);

            
            det = LU  (numNodes, numConst, kr, fr, u, in, nodesC);            
            //mexPrintf("Determinante = %f\n",det);
            
            if (det>0) {

				tension_stress(n, numNodes, u, C, elasticity, length, vector, F, sumConst, area, conectividade, aalfa);
				ccrit(n, sumConst, length, vector, numNodes, numConst, Kg, kg, in, area, elasticity, conectividade, C, aalfa);
				//compliance(kr, u, nodesC, numNodes, comp, fr);

				//mexPrintf("%f",*comp);

			}else{
                for (i = 0; i < (numNodes*D-numConst); i++) {
                    for (j = 0; j < (numNodes*D-numConst); j++) {
                        kr[i][j] = 0.5;
                        mr[i][j] = 10;
                        kg[i][j] = 0.6;
                    }
                }
                
                int cont=0;
                for (i=0;i<numNodes*Dkl-numConst;i++){
                    for (j=0;j<numNodes*Dkl-numConst;j++){  
                        Ke[cont]=kr[j][i];
                        Mmat[cont]=mr[j][i];
                        Kg[cont]=kg[j][i];
                        cont++;
                    }
                }

				for (i=0;i<numNodes*Dkl;i++){
					u[i]=9999.0;
				}

				for (i=0;i<n;i++){
					sumConst[i]=99999999.0;
				}

				//*comp=9999.0;
			}

	}else{
        
        for (i = 0; i < (numNodes*D-numConst); i++) {
            for (j = 0; j < (numNodes*D-numConst); j++) {
                kr[i][j] = 0.5;
                mr[i][j] = 10;
                kg[i][j] = 0.6;
            }
        }
                
        int cont=0;
        for (i=0;i<numNodes*Dkl-numConst;i++){
            for (j=0;j<numNodes*Dkl-numConst;j++){  
                Ke[cont]=kr[j][i];
                Mmat[cont]=mr[j][i];
                Kg[cont]=kg[j][i];
                cont++;
            }
        }

		for (i=0;i<numNodes*Dkl;i++){
			u[i]=99.0;
		}

		for (i=0;i<n;i++){
			sumConst[i]=999999.0;
		}

		//*comp=9999.0;
	}

	/*for (i=0;i<numNodes*Dkl;i++){
		mexPrintf("%f ",nodesC[i]);
	}
	mexPrintf("\n");
	for (i=0;i<numNodes*Dkl;i++){
		mexPrintf("%f ",nodesC[i]);
	}*/


	/*for (i=0;i<(numNodes-3);i++){
		nodesC[3*i]=0;
		nodesC[3*i+1]=1;
		nodesC[3*i+2]=0;
	}

	for (i=(numNodes-3);i<numNodes;i++){
		nodesC[3*i]=1;
		nodesC[3*i+1]=1;
		nodesC[3*i+2]=1;
	}*/


    /*free (area);
	free (elasticity);
	free (length);
	free (nodesC);
	free (nodesL);
	free (in);
	free (F);
	free (fr);
	for (i = 0; i < (numNodes*Dkl-numConst+1); i++)
		free (kr[i]);
	free (kr);
	for (i = 0; i < n; i++) {
		free (vector[i]);
		free (coefXY[i]);
		free (conectividade[i]);
	}
	free (vector);
	free (coefXY);
	free (conectividade);
	for (i = 0; i < numNodes; i++)
		free (C[i]);
	free (C);*/

}
/**
 * Reduced_matrix function
 */
void reduced_matrix (int numNodes, int n, double * area, double * elasticity, double * length, int ** conectividade, double ** C, double ** kr, int numConst, int * in, int ** vector, double ** coefXY, double * fr, double * nodesC, double * Ke, double ro, double * aalfa, double ** mr, double *Mmat) {
	int i,j;

	numConst=0;
	for (i=0;i<Dkl*numNodes;i++){
		if (nodesC[i]==1){
			numConst+=1;
		}
		//mexPrintf("nodesC[%d] = %f\n",i+1,nodesC[i]);

	}
	//mexPrintf("numConst = %f\n",numConst);
	//numConst=27;
	//mexPrintf("Entrou aqui %d\n",numConst);
	//double cx, cy, cz, ea, rl, ro;
	double cx, cy, cz, mm, kk, ea, a, l;
	//ro = 7850.0;

	//int incog = numNodes*Dkl-numConst;
	//mexPrintf("numConst2: %d\n",numConst);
	//incog = 18 * 3;
	// Inicializacao da matriz kr
	for (i = 0; i < (numNodes*D-numConst); i++) {
		for (j = 0; j < (numNodes*D-numConst); j++) {
		    kr[i][j] = 0.;
		    mr[i][j] = 0.;
		}
	}



	/*mexPrintf("Conectividades\n");
	for (i=0;i<n;i++) {
		mexPrintf("%d %d\n",conectividade[i][0],conectividade[i][1]);
	}

	mexPrintf("nodesC\n");
	for (i=0;i<numNodes;i++) {
		mexPrintf("%d\n",nodesC[i]);
	}

	mexPrintf("vector\n");
	for (i=0;i<n;i++) {
		for (j=0;j<D;j++){
			mexPrintf("%d\n",vector[i][j]);
		}
		mexPrintf("\n");
	}*/


	for(i = 0; i < n; i++) {



			//mexPrintf("\nBarra %d aalfa %f\n",i+1,aalfa[i]);

			//mexPrintf("Conectividade = %d %d\n",conectividade[i][0],conectividade[i][1]);
			//mexPrintf("%d %d %d %d %d %d\n",in[vector[i][0]],in[vector[i][1]],in[vector[i][2]],in[vector[i][3]],in[vector[i][4]],in[vector[i][5]]);
		    a = area[i];
		    l = length[i];
	        cx = (C[conectividade[i][1]-1][0] - C[conectividade[i][0]-1][0]) / length[i];
	        cy = (C[conectividade[i][1]-1][1] - C[conectividade[i][0]-1][1]) / length[i];
	        cz = (C[conectividade[i][1]-1][2] - C[conectividade[i][0]-1][2]) / length[i];
	        //mexPrintf("Cx Cy Cz = %f %f %f\n",cx,cy,cz);
	        //mexPrintf("EA = %f\n",ea);
			mm = cos(aalfa[i]);
	        kk = sin(aalfa[i]);
					//mexPrintf("sen cos = %f %f\n",kk,mm);
					//mexPrintf("ro = %f\n",ro);
	        //mm=1.0;
	        //kk=0.0;
	        ea = elasticity[i] * area[i] / length[i];




			if ((cx == 0.0) && (cz == 0.0)) {

				//mexPrintf("Entrou aqui (barra vertical)\n");

				kr[in[vector[i][0]]][in[vector[i][0]]] += 0.0;
				kr[in[vector[i][0]]][in[vector[i][1]]] += 0.0;
				kr[in[vector[i][0]]][in[vector[i][2]]] += 0.0;
				kr[in[vector[i][0]]][in[vector[i][3]]] += 0.0;
				kr[in[vector[i][0]]][in[vector[i][4]]] += 0.0;
				kr[in[vector[i][0]]][in[vector[i][5]]] += 0.0;

				kr[in[vector[i][1]]][in[vector[i][0]]] += 0.0;
				kr[in[vector[i][1]]][in[vector[i][1]]] += ea * cy*cy;
				kr[in[vector[i][1]]][in[vector[i][2]]] += 0.0;
				kr[in[vector[i][1]]][in[vector[i][3]]] += 0.0;
				kr[in[vector[i][1]]][in[vector[i][4]]] += -ea * cy*cy;
				kr[in[vector[i][1]]][in[vector[i][5]]] += 0.0;

				kr[in[vector[i][2]]][in[vector[i][0]]] += 0.0;
				kr[in[vector[i][2]]][in[vector[i][1]]] += 0.0;
				kr[in[vector[i][2]]][in[vector[i][2]]] += 0.0;
				kr[in[vector[i][2]]][in[vector[i][3]]] += 0.0;
				kr[in[vector[i][2]]][in[vector[i][4]]] += 0.0;
				kr[in[vector[i][2]]][in[vector[i][5]]] += 0.0;

				kr[in[vector[i][3]]][in[vector[i][0]]] += 0.0;
				kr[in[vector[i][3]]][in[vector[i][1]]] += 0.0;
				kr[in[vector[i][3]]][in[vector[i][2]]] += 0.0;
				kr[in[vector[i][3]]][in[vector[i][3]]] += 0.0;
				kr[in[vector[i][3]]][in[vector[i][4]]] += 0.0;
				kr[in[vector[i][3]]][in[vector[i][5]]] += 0.0;

				kr[in[vector[i][4]]][in[vector[i][0]]] += 0.0;
				kr[in[vector[i][4]]][in[vector[i][1]]] += -ea * cy*cy;
				kr[in[vector[i][4]]][in[vector[i][2]]] += 0.0;
				kr[in[vector[i][4]]][in[vector[i][3]]] += 0.0;
				kr[in[vector[i][4]]][in[vector[i][4]]] += ea * cy*cy;
				kr[in[vector[i][4]]][in[vector[i][5]]] += 0.0;

				kr[in[vector[i][5]]][in[vector[i][0]]] += 0.0;
				kr[in[vector[i][5]]][in[vector[i][1]]] += 0.0;
				kr[in[vector[i][5]]][in[vector[i][2]]] += 0.0;
				kr[in[vector[i][5]]][in[vector[i][3]]] += 0.0;
				kr[in[vector[i][5]]][in[vector[i][4]]] += 0.0;
				kr[in[vector[i][5]]][in[vector[i][5]]] += 0.0;


				// -----------------------------------------------------------------//



				mr[in[vector[i][0]]][in[vector[i][0]]] += a*(((kk*kk*l*ro)/(3))+((l*mm*mm*ro)/(3)))*cy*cy;
				mr[in[vector[i][0]]][in[vector[i][1]]] += 0.0;
				mr[in[vector[i][0]]][in[vector[i][2]]] += 0.0;
				mr[in[vector[i][0]]][in[vector[i][3]]] += a*(((kk*kk*l*ro)/(6))+((l*mm*mm*ro)/(6)))*cy*cy;
				mr[in[vector[i][0]]][in[vector[i][4]]] += 0.0;
				mr[in[vector[i][0]]][in[vector[i][5]]] += 0.0;

				mr[in[vector[i][1]]][in[vector[i][0]]] += 0.0;
				mr[in[vector[i][1]]][in[vector[i][1]]] += ((a*l*ro*cy*cy)/(3));
				mr[in[vector[i][1]]][in[vector[i][2]]] += 0.0;
				mr[in[vector[i][1]]][in[vector[i][3]]] += 0.0;
				mr[in[vector[i][1]]][in[vector[i][4]]] += ((a*l*ro*cy*cy)/(6));
				mr[in[vector[i][1]]][in[vector[i][5]]] += 0.0;

				mr[in[vector[i][2]]][in[vector[i][0]]] += 0.0;
				mr[in[vector[i][2]]][in[vector[i][1]]] += 0.0;
				mr[in[vector[i][2]]][in[vector[i][2]]] += a*(((kk*kk*l*ro)/(3))+((l*mm*mm*ro)/(3)));
				mr[in[vector[i][2]]][in[vector[i][3]]] += 0.0;
				mr[in[vector[i][2]]][in[vector[i][4]]] += 0.0;
				mr[in[vector[i][2]]][in[vector[i][5]]] += a*(((kk*kk*l*ro)/(6))+((l*mm*mm*ro)/(6)));

				mr[in[vector[i][3]]][in[vector[i][0]]] += a*(((kk*kk*l*ro)/(6))+((l*mm*mm*ro)/(6)))*cy*cy;
				mr[in[vector[i][3]]][in[vector[i][1]]] += 0.0;
				mr[in[vector[i][3]]][in[vector[i][2]]] += 0.0;
				mr[in[vector[i][3]]][in[vector[i][3]]] += a*(((kk*kk*l*ro)/(3))+((l*mm*mm*ro)/(3)))*cy*cy;
				mr[in[vector[i][3]]][in[vector[i][4]]] += 0.0;
				mr[in[vector[i][3]]][in[vector[i][5]]] += 0.0;

				mr[in[vector[i][4]]][in[vector[i][0]]] += 0.0;
				mr[in[vector[i][4]]][in[vector[i][1]]] += ((a*l*ro*cy*cy)/(6));
				mr[in[vector[i][4]]][in[vector[i][2]]] += 0.0;
				mr[in[vector[i][4]]][in[vector[i][3]]] += 0.0;
				mr[in[vector[i][4]]][in[vector[i][4]]] += ((a*l*ro*cy*cy)/(3));
				mr[in[vector[i][4]]][in[vector[i][5]]] += 0.0;

				mr[in[vector[i][5]]][in[vector[i][0]]] += 0.0;
				mr[in[vector[i][5]]][in[vector[i][1]]] += 0.0;
				mr[in[vector[i][5]]][in[vector[i][2]]] += a*(((kk*kk*l*ro)/(6))+((l*mm*mm*ro)/(6)));
				mr[in[vector[i][5]]][in[vector[i][3]]] += 0.0;
				mr[in[vector[i][5]]][in[vector[i][4]]] += 0.0;
				mr[in[vector[i][5]]][in[vector[i][5]]] += a*(((kk*kk*l*ro)/(3))+((l*mm*mm*ro)/(3)));



			} else {

				//mexPrintf("Entrou aqui (barra nao vertical))\n");

				kr[in[vector[i][0]]][in[vector[i][0]]] += cx*cx*ea;
				kr[in[vector[i][0]]][in[vector[i][1]]] += cx*cy*ea;
				kr[in[vector[i][0]]][in[vector[i][2]]] += cx*cz*ea;
				kr[in[vector[i][0]]][in[vector[i][3]]] += -cx*cx*ea;
				kr[in[vector[i][0]]][in[vector[i][4]]] += -cx*cy*ea;
				kr[in[vector[i][0]]][in[vector[i][5]]] += -cx*cz*ea;

				kr[in[vector[i][1]]][in[vector[i][0]]] += cx*cy*ea;
				kr[in[vector[i][1]]][in[vector[i][1]]] += cy*cy*ea;
				kr[in[vector[i][1]]][in[vector[i][2]]] += cy*cz*ea;
				kr[in[vector[i][1]]][in[vector[i][3]]] += -cx*cy*ea;
				kr[in[vector[i][1]]][in[vector[i][4]]] += -cy*cy*ea;
				kr[in[vector[i][1]]][in[vector[i][5]]] += -cy*cz*ea;

				kr[in[vector[i][2]]][in[vector[i][0]]] += cx*cz*ea;
				kr[in[vector[i][2]]][in[vector[i][1]]] += cy*cz*ea;
				kr[in[vector[i][2]]][in[vector[i][2]]] += cz*cz*ea;
				kr[in[vector[i][2]]][in[vector[i][3]]] += -cx*cz*ea;
				kr[in[vector[i][2]]][in[vector[i][4]]] += -cy*cz*ea;
				kr[in[vector[i][2]]][in[vector[i][5]]] += -cz*cz*ea;

				kr[in[vector[i][3]]][in[vector[i][0]]] += -cx*cx*ea;
				kr[in[vector[i][3]]][in[vector[i][1]]] += -cx*cy*ea;
				kr[in[vector[i][3]]][in[vector[i][2]]] += -cx*cz*ea;
				kr[in[vector[i][3]]][in[vector[i][3]]] += cx*cx*ea;
				kr[in[vector[i][3]]][in[vector[i][4]]] += cx*cy*ea;
				kr[in[vector[i][3]]][in[vector[i][5]]] += cx*cz*ea;

				kr[in[vector[i][4]]][in[vector[i][0]]] += -cx*cy*ea;
				kr[in[vector[i][4]]][in[vector[i][1]]] += -cy*cy*ea;
				kr[in[vector[i][4]]][in[vector[i][2]]] += -cy*cz*ea;
				kr[in[vector[i][4]]][in[vector[i][3]]] += cx*cy*ea;
				kr[in[vector[i][4]]][in[vector[i][4]]] += cy*cy*ea;
				kr[in[vector[i][4]]][in[vector[i][5]]] += cy*cz*ea;

				kr[in[vector[i][5]]][in[vector[i][0]]] += -cx*cz*ea;
				kr[in[vector[i][5]]][in[vector[i][1]]] += -cy*cz*ea;
				kr[in[vector[i][5]]][in[vector[i][2]]] += -cz*cz*ea;
				kr[in[vector[i][5]]][in[vector[i][3]]] += cx*cz*ea;
				kr[in[vector[i][5]]][in[vector[i][4]]] += cy*cz*ea;
				kr[in[vector[i][5]]][in[vector[i][5]]] += cz*cz*ea;


					// -----------------------------------------------------------------//


				mr[in[vector[i][0]]][in[vector[i][0]]] += ((-a*(kk*kk+mm*mm)*l*ro*(cy*cy-1)*cz*cz)/(3*(cx*cx+cz*cz)))+((a*l*ro*cx*cx)/(3))+((a*(kk*kk+mm*mm)*l*ro*cy*cy)/(3));
				mr[in[vector[i][0]]][in[vector[i][1]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cx*cy)/(3));
				mr[in[vector[i][0]]][in[vector[i][2]]] += ((a*(kk*kk+mm*mm)*l*ro*cx*(cy*cy-1)*cz)/(3*(cx*cx+cz*cz)))+((a*l*ro*cx*cz)/(3));
				mr[in[vector[i][0]]][in[vector[i][3]]] += ((-a*(kk*kk+mm*mm)*l*ro*(cy*cy-1)*cz*cz)/(6*(cx*cx+cz*cz)))+((a*l*ro*cx*cx)/(6))+((a*(kk*kk+mm*mm)*l*ro*cy*cy)/(6));
				mr[in[vector[i][0]]][in[vector[i][4]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cx*cy)/(6));
				mr[in[vector[i][0]]][in[vector[i][5]]] += ((a*(kk*kk+mm*mm)*l*ro*cx*(cy*cy-1)*cz)/(6*(cx*cx+cz*cz)))+((a*l*ro*cx*cz)/(6));

				mr[in[vector[i][1]]][in[vector[i][0]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cx*cy)/(3));
				mr[in[vector[i][1]]][in[vector[i][1]]] += ((a*l*ro*((kk*kk+mm*mm)*cx*cx+cy*cy+(kk*kk+mm*mm)*cz*cz))/(3));
				mr[in[vector[i][1]]][in[vector[i][2]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cy*cz)/(3));
				mr[in[vector[i][1]]][in[vector[i][3]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cx*cy)/(6));
				mr[in[vector[i][1]]][in[vector[i][4]]] += ((a*l*ro*((kk*kk+mm*mm)*cx*cx+cy*cy+(kk*kk+mm*mm)*cz*cz))/(6));
				mr[in[vector[i][1]]][in[vector[i][5]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cy*cz)/(6));

				mr[in[vector[i][2]]][in[vector[i][0]]] += ((a*(kk*kk+mm*mm)*l*ro*cx*(cy*cy-1)*cz)/(3*(cx*cx+cz*cz)))+((a*l*ro*cx*cz)/(3));
				mr[in[vector[i][2]]][in[vector[i][1]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cy*cz)/(3));
				mr[in[vector[i][2]]][in[vector[i][2]]] += ((a*(kk*kk+mm*mm)*l*ro*(cy*cy-1)*cz*cz)/(3*(cx*cx+cz*cz)))+((a*l*ro*cz*cz)/(3))+((a*(kk*kk+mm*mm)*l*ro)/(3));
				mr[in[vector[i][2]]][in[vector[i][3]]] += ((a*(kk*kk+mm*mm)*l*ro*cx*(cy*cy-1)*cz)/(6*(cx*cx+cz*cz)))+((a*l*ro*cx*cz)/(6));
				mr[in[vector[i][2]]][in[vector[i][4]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cy*cz)/(6));
				mr[in[vector[i][2]]][in[vector[i][5]]] += ((a*(kk*kk+mm*mm)*l*ro*(cy*cy-1)*cz*cz)/(6*(cx*cx+cz*cz)))+((a*l*ro*cz*cz)/(6))+((a*(kk*kk+mm*mm)*l*ro)/(6));

				mr[in[vector[i][3]]][in[vector[i][0]]] += ((-a*(kk*kk+mm*mm)*l*ro*(cy*cy-1)*cz*cz)/(6*(cx*cx+cz*cz)))+((a*l*ro*cx*cx)/(6))+((a*(kk*kk+mm*mm)*l*ro*cy*cy)/(6));
				mr[in[vector[i][3]]][in[vector[i][1]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cx*cy)/(6));
				mr[in[vector[i][3]]][in[vector[i][2]]] += ((a*(kk*kk+mm*mm)*l*ro*cx*(cy*cy-1)*cz)/(6*(cx*cx+cz*cz)))+((a*l*ro*cx*cz)/(6));
				mr[in[vector[i][3]]][in[vector[i][3]]] += ((-a*(kk*kk+mm*mm)*l*ro*(cy*cy-1)*cz*cz)/(3*(cx*cx+cz*cz)))+((a*l*ro*cx*cx)/(3))+((a*(kk*kk+mm*mm)*l*ro*cy*cy)/(3));
				mr[in[vector[i][3]]][in[vector[i][4]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cx*cy)/(3));
				mr[in[vector[i][3]]][in[vector[i][5]]] += ((a*(kk*kk+mm*mm)*l*ro*cx*(cy*cy-1)*cz)/(3*(cx*cx+cz*cz)))+((a*l*ro*cx*cz)/(3));

				mr[in[vector[i][4]]][in[vector[i][0]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cx*cy)/(6));
				mr[in[vector[i][4]]][in[vector[i][1]]] += ((a*l*ro*((kk*kk+mm*mm)*cx*cx+cy*cy+(kk*kk+mm*mm)*cz*cz))/(6));
				mr[in[vector[i][4]]][in[vector[i][2]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cy*cz)/(6));
				mr[in[vector[i][4]]][in[vector[i][3]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cx*cy)/(3));
				mr[in[vector[i][4]]][in[vector[i][4]]] += ((a*l*ro*((kk*kk+mm*mm)*cx*cx+cy*cy+(kk*kk+mm*mm)*cz*cz))/(3));
				mr[in[vector[i][4]]][in[vector[i][5]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cy*cz)/(3));

				mr[in[vector[i][5]]][in[vector[i][0]]] += ((a*(kk*kk+mm*mm)*l*ro*cx*(cy*cy-1)*cz)/(6*(cx*cx+cz*cz)))+((a*l*ro*cx*cz)/(6));
				mr[in[vector[i][5]]][in[vector[i][1]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cy*cz)/(6));
				mr[in[vector[i][5]]][in[vector[i][2]]] += ((a*(kk*kk+mm*mm)*l*ro*(cy*cy-1)*cz*cz)/(6*(cx*cx+cz*cz)))+((a*l*ro*cz*cz)/(6))+((a*(kk*kk+mm*mm)*l*ro)/(6));
				mr[in[vector[i][5]]][in[vector[i][3]]] += ((a*(kk*kk+mm*mm)*l*ro*cx*(cy*cy-1)*cz)/(3*(cx*cx+cz*cz)))+((a*l*ro*cx*cz)/(3));
				mr[in[vector[i][5]]][in[vector[i][4]]] += ((-a*(kk*kk+mm*mm-1)*l*ro*cy*cz)/(3));
				mr[in[vector[i][5]]][in[vector[i][5]]] += ((a*(kk*kk+mm*mm)*l*ro*(cy*cy-1)*cz*cz)/(3*(cx*cx+cz*cz)))+((a*l*ro*cz*cz)/(3))+((a*(kk*kk+mm*mm)*l*ro)/(3));


			}

            /*mexPrintf("Elemento %d \nCx = %f\nCy = %f\nCz = %f\nEA/L = %f\n",i+1,cx,cy,cz,ea);
			if ((cx!=0) && (cz!=0)) {
                mexPrintf("%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n",cx*cx*ea,cx*cy*ea,cx*cz*ea,-cx*cx*ea,-cx*cy*ea,-cx*cz*ea,cx*cy*ea,cy*cy*ea,cy*cz*ea,-cx*cy*ea,-cy*cy*ea,-cy*cz*ea,cx*cz*ea,cy*cz*ea,cz*cz*ea,-cx*cz*ea,-cy*cz*ea,-cz*cz*ea,-cx*cx*ea,-cx*cy*ea,-cx*cz*ea,cx*cx*ea,cx*cy*ea,cx*cz*ea,-cx*cy*ea,-cy*cy*ea,-cy*cz*ea,cx*cy*ea,cy*cy*ea,cy*cz*ea,-cx*cz*ea,-cy*cz*ea,-cz*cz*ea,cx*cz*ea,cy*cz*ea,cz*cz*ea);
            }else{
                mexPrintf("%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n",0,0,0,0,0,0,0,ea * cy*cy,0,0,-ea * cy*cy,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-ea * cy*cy,0,0,ea * cy*cy,0,0,0,0,0,0,0);
            }*/
            
            //mexPrintf("%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n",kr[in[vector[i][0]]][in[vector[i][0]]],kr[in[vector[i][0]]][in[vector[i][1]]],kr[in[vector[i][0]]][in[vector[i][2]]],kr[in[vector[i][0]]][in[vector[i][3]]],kr[in[vector[i][0]]][in[vector[i][4]]],kr[in[vector[i][0]]][in[vector[i][5]]],kr[in[vector[i][1]]][in[vector[i][0]]],kr[in[vector[i][1]]][in[vector[i][1]]],kr[in[vector[i][1]]][in[vector[i][2]]],kr[in[vector[i][1]]][in[vector[i][3]]],kr[in[vector[i][1]]][in[vector[i][4]]],kr[in[vector[i][1]]][in[vector[i][5]]],kr[in[vector[i][2]]][in[vector[i][0]]],kr[in[vector[i][2]]][in[vector[i][1]]],kr[in[vector[i][2]]][in[vector[i][2]]],kr[in[vector[i][2]]][in[vector[i][3]]],kr[in[vector[i][2]]][in[vector[i][4]]],kr[in[vector[i][2]]][in[vector[i][5]]],kr[in[vector[i][3]]][in[vector[i][0]]],kr[in[vector[i][3]]][in[vector[i][1]]],kr[in[vector[i][3]]][in[vector[i][2]]],kr[in[vector[i][3]]][in[vector[i][3]]],kr[in[vector[i][3]]][in[vector[i][4]]],kr[in[vector[i][3]]][in[vector[i][5]]],kr[in[vector[i][4]]][in[vector[i][0]]],kr[in[vector[i][4]]][in[vector[i][1]]],kr[in[vector[i][4]]][in[vector[i][2]]],kr[in[vector[i][4]]][in[vector[i][3]]],kr[in[vector[i][4]]][in[vector[i][4]]],kr[in[vector[i][4]]][in[vector[i][5]]],kr[in[vector[i][5]]][in[vector[i][0]]],kr[in[vector[i][5]]][in[vector[i][1]]],kr[in[vector[i][5]]][in[vector[i][2]]],kr[in[vector[i][5]]][in[vector[i][3]]],kr[in[vector[i][5]]][in[vector[i][4]]],kr[in[vector[i][5]]][in[vector[i][5]]]);
            // matrizes acumuladas
			//mexPrintf("%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n",mr[in[vector[i][0]]][in[vector[i][0]]],mr[in[vector[i][0]]][in[vector[i][1]]],mr[in[vector[i][0]]][in[vector[i][2]]],mr[in[vector[i][0]]][in[vector[i][3]]],mr[in[vector[i][0]]][in[vector[i][4]]],mr[in[vector[i][0]]][in[vector[i][5]]],mr[in[vector[i][1]]][in[vector[i][0]]],mr[in[vector[i][1]]][in[vector[i][1]]],mr[in[vector[i][1]]][in[vector[i][2]]],mr[in[vector[i][1]]][in[vector[i][3]]],mr[in[vector[i][1]]][in[vector[i][4]]],mr[in[vector[i][1]]][in[vector[i][5]]],mr[in[vector[i][2]]][in[vector[i][0]]],mr[in[vector[i][2]]][in[vector[i][1]]],mr[in[vector[i][2]]][in[vector[i][2]]],mr[in[vector[i][2]]][in[vector[i][3]]],mr[in[vector[i][2]]][in[vector[i][4]]],mr[in[vector[i][2]]][in[vector[i][5]]],mr[in[vector[i][3]]][in[vector[i][0]]],mr[in[vector[i][3]]][in[vector[i][1]]],mr[in[vector[i][3]]][in[vector[i][2]]],mr[in[vector[i][3]]][in[vector[i][3]]],mr[in[vector[i][3]]][in[vector[i][4]]],mr[in[vector[i][3]]][in[vector[i][5]]],mr[in[vector[i][4]]][in[vector[i][0]]],mr[in[vector[i][4]]][in[vector[i][1]]],mr[in[vector[i][4]]][in[vector[i][2]]],mr[in[vector[i][4]]][in[vector[i][3]]],mr[in[vector[i][4]]][in[vector[i][4]]],mr[in[vector[i][4]]][in[vector[i][5]]],mr[in[vector[i][5]]][in[vector[i][0]]],mr[in[vector[i][5]]][in[vector[i][1]]],mr[in[vector[i][5]]][in[vector[i][2]]],mr[in[vector[i][5]]][in[vector[i][3]]],mr[in[vector[i][5]]][in[vector[i][4]]],mr[in[vector[i][5]]][in[vector[i][5]]]);

			//mexPrintf("\n\n");




		// ----------------------------------------------------------------------------------




		// ----------------------------------------------------------------------------------


}

//for (i=0;i<numNodes*Dkl-numConst;i++) {
//    mr[i][i]+=45;
//}


int cont=0;
for (i=0;i<numNodes*Dkl-numConst;i++){
	for (j=0;j<numNodes*Dkl-numConst;j++){
		Ke[cont]=0;
		Mmat[cont]=0;
		Ke[cont]=kr[j][i];
		Mmat[cont]=mr[j][i];
		cont++;
	}
}





	//mexPrintf("%d \n",sizeof(kr));

	/*for (i=0;i<numNodes*Dkl-numConst;i++) {
		for (j=0;j<numNodes*Dkl-numConst;j++) {
			mexPrintf("%f ", kr[i][j]);
		}
		mexPrintf("\n");
	}*/

	/*for (i=0;i<numNodes*Dkl-numConst;i++) {
		for (j=0;j<numNodes*Dkl-numConst;j++) {
			mexPrintf("%f ", mr[i][j]);
		}
		mexPrintf("\n");
	}*/


	/*
	for (i=0;i<n;i++) {
		for (j=0;j<6;j++) {
			mexPrintf("%d \n", vector[i][j]);
		}
		mexPrintf("\n");
	}
  */

	//mexPrintf("%d \n\n",n);

	/*
	int elem=2;

	mexPrintf("%d %d\n",in[vector[elem][0]],in[vector[elem][0]]);
	mexPrintf("%d %d\n",in[vector[elem][0]],in[vector[elem][1]]);
	mexPrintf("%d %d\n",in[vector[elem][0]],in[vector[elem][2]]);
	mexPrintf("%d %d\n",in[vector[elem][0]],in[vector[elem][3]]);
	mexPrintf("%d %d\n",in[vector[elem][0]],in[vector[elem][4]]);
	mexPrintf("%d %d\n\n",in[vector[elem][0]],in[vector[elem][5]]);

	mexPrintf("%d %d\n",in[vector[elem][1]],in[vector[elem][0]]);
	mexPrintf("%d %d\n",in[vector[elem][1]],in[vector[elem][1]]);
	mexPrintf("%d %d\n",in[vector[elem][1]],in[vector[elem][2]]);
	mexPrintf("%d %d\n",in[vector[elem][1]],in[vector[elem][3]]);
	mexPrintf("%d %d\n",in[vector[elem][1]],in[vector[elem][4]]);
	mexPrintf("%d %d\n\n",in[vector[elem][1]],in[vector[elem][5]]);

	mexPrintf("%d %d\n",in[vector[elem][2]],in[vector[elem][0]]);
	mexPrintf("%d %d\n",in[vector[elem][2]],in[vector[elem][1]]);
	mexPrintf("%d %d\n",in[vector[elem][2]],in[vector[elem][2]]);
	mexPrintf("%d %d\n",in[vector[elem][2]],in[vector[elem][3]]);
	mexPrintf("%d %d\n",in[vector[elem][2]],in[vector[elem][4]]);
	mexPrintf("%d %d\n\n",in[vector[elem][2]],in[vector[elem][5]]);

	mexPrintf("%d %d\n",in[vector[elem][3]],in[vector[elem][0]]);
	mexPrintf("%d %d\n",in[vector[elem][3]],in[vector[elem][1]]);
	mexPrintf("%d %d\n",in[vector[elem][3]],in[vector[elem][2]]);
	mexPrintf("%d %d\n",in[vector[elem][3]],in[vector[elem][3]]);
	mexPrintf("%d %d\n",in[vector[elem][3]],in[vector[elem][4]]);
	mexPrintf("%d %d\n\n",in[vector[elem][3]],in[vector[elem][5]]);

	mexPrintf("%d %d\n",in[vector[elem][4]],in[vector[elem][0]]);
	mexPrintf("%d %d\n",in[vector[elem][4]],in[vector[elem][1]]);
	mexPrintf("%d %d\n",in[vector[elem][4]],in[vector[elem][2]]);
	mexPrintf("%d %d\n",in[vector[elem][4]],in[vector[elem][3]]);
	mexPrintf("%d %d\n",in[vector[elem][4]],in[vector[elem][4]]);
	mexPrintf("%d %d\n\n",in[vector[elem][4]],in[vector[elem][5]]);

	mexPrintf("%d %d\n",in[vector[elem][5]],in[vector[elem][0]]);
	mexPrintf("%d %d\n",in[vector[elem][5]],in[vector[elem][1]]);
	mexPrintf("%d %d\n",in[vector[elem][5]],in[vector[elem][2]]);
	mexPrintf("%d %d\n",in[vector[elem][5]],in[vector[elem][3]]);
	mexPrintf("%d %d\n",in[vector[elem][5]],in[vector[elem][4]]);
	mexPrintf("%d %d\n\n",in[vector[elem][5]],in[vector[elem][5]]);
  */


/*printf("\nKr...\n");
    for (i = 0; i < 9; i++) {
		for(j = 0; j < 9; j++) {
			printf ("%d\t\t",kr[i][j]);

		}printf("\n");
    }

FILE * out1a = fopen("out1a.txt", "w");
	if (out1a == NULL) {
        printf("\nError 'out1'..\n");
        exit(1);
	}

for (i=0; i<9; i++) {
		//fprintf(out1, "\n\t");
		for(j=0; j<9; j++){
			fprintf (out1a, "%lf\t",kr[i][j]);
		}fprintf(out1a, "\n");
	}
*/
}

/**
 * LU function
 */
double LU (int numNodes, int numConst, double ** kr, double * fr, double * u, int * in, double * nodesC) {

	int n, i, j, k, l;
	//double AAA[n][n+1], x[n], termo, m;
	double termo, m;
	double **AAA;
	double *x;

	numConst=0;
	for (i=0;i<Dkl*numNodes;i++){
		if (nodesC[i]==1){
			numConst+=1;
		}
		//mexPrintf("%f ",nodesC[i]);
	}
	//mexPrintf("\n");


	int c = 0;
	n = numNodes*D-numConst;
	AAA = (double**) mxMalloc(sizeof(double*)*n);
	x = (double*) mxMalloc (sizeof(double)*n);
	for(i = 0; i < n; i++) {
        AAA[i] = (double*) mxMalloc(sizeof(double)*(n+1));
	}

	//printf ("%lf\t xxxxxxxxxxxxxxxxxxxxxxxxxxxxx",x[n]);

	// Definir A
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			AAA[i][j] = kr[i][j];
			//mexPrintf("%f ",AAA[i][j]);
		}
		//mexPrintf("\n");
	}

	// acoplamento fr
	for (j = 0; j < n; j++) {
		AAA[j][n] = fr[j];
	}

	/*for (j = 0; j < n; j++) {
		mexPrintf("%f \n",fr[j]);
	}

	mexPrintf("\n\n\n\n");*/


	// Implementando MÃ©todo de Gauss
	for (k = 0; k < n-1; k++) {
		for (i = k+1; i < n; i++) {
			// Multiplicadores
			m = -1 * (AAA[i][k] / AAA[k][k]);
			for (j = 0; j < n+1; j++) {
				AAA[i][j] = (AAA[k][j] * m) + AAA[i][j];
			}
		}
	}
    
    double det=1.0;
    
   /* for (i=0; i<n ; i++){
        for (j=0 ; j<n ; j++){
            mexPrintf("%f ",AAA[i][j]);
        }
        mexPrintf("\n");
    }*/
    
    for (i=0;i<n;i++){
        det=det*AAA[i][i];
    }

	// Resolvendo o sistema
	for (i = 0; i < n; i++){
		termo = 0;
		l = n - i;
		for (j=l; j<n;j++){
			termo = termo + (x[j] * AAA[n-i-1][j]);
		}
		x[n-i-1] = (AAA[n-1-i][n] - termo) / AAA[n-i-1][n-i-1];
	}

	//impressao de U

	for (i = 0; i < (numNodes*D); i++) {
		u[i] = 0;
		if(in[i] < (numNodes*D-numConst)){
			u[i] = x[c];
			c++;
		}
	}

	/*for (i = 0; i < (numNodes*D); i++){
		mexPrintf("%f\n",u[i]);
	}*/


	// ---------------------------------------------------------------------------------------------------------------------------------------

	// GSL - Invertendo a matriz 'matmassar'
    /*
	int incog=8;
	double lambda[incog];
	double inversa[incog][incog];
	// Define the dimension n of the matrix and the signum s (for LU decomposition)
	int dim = 8;
	int s;

	// Define all the used matrices
	gsl_matrix * matrix = gsl_matrix_alloc (dim, dim);
	gsl_matrix * inverse = gsl_matrix_alloc (dim, dim);
	gsl_permutation * perm = gsl_permutation_alloc (dim);

	// Fill the matrix m
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			gsl_matrix_set (matrix, i, j, mr[i][j]);
		}
	}

	// Make LU decomposition of matrix m
	gsl_linalg_LU_decomp (matrix, perm, &s);

	// Invert the matrix m
	gsl_linalg_LU_invert (matrix, perm, inverse);
	//printf("Inversa\n");
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			inversa[i][j] = gsl_matrix_get(inverse, i, j);
		}
	}

	//Multiplicando a inversa de 'mr' (inversa) por 'kr'
	int cont5;
	double result[incog][incog], AUX;
	for (i = 0; i < incog; i++) {
		for (j = 0; j < incog; j++) {
			AUX = 0.;
			for (cont5 = 0; cont5 < incog; cont5++) {
				AUX += inversa[i][cont5] * kr[cont5][j];
			}
			result[i][j] = AUX;
		}
	}

	// Somando as matrizes 'mr' e 'kr'
    double matriz[incog][incog], data[incog*incog];
    for (i = 0; i < incog; i++) {
        for (j = 0; j < incog; j++) {
            matriz[i][j] = kr[i][j] + mr[i][j];
        }
    }

	// Transformando a matriz 'result[incog][incog]' em um vetor 'data[incog*incog]'
	int datacount = 0;
    for (i = 0; i < incog; i++) {
        for (k = 0; k < incog; k++) {
			data[datacount] = result[i][k];
			datacount++;
		}
	}


	// GSL - Calculo dos autovalores e autovetores
	gsl_matrix_view m2 = gsl_matrix_view_array (data, incog, incog);
	gsl_vector_complex *eval = gsl_vector_complex_alloc (incog);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc (incog, incog);
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (incog);
	gsl_eigen_nonsymmv (&m2.matrix, eval, evec, w);
	gsl_eigen_nonsymmv_free (w);
	gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
	{
		int count;

		for (count = 0; count < incog; count++){
			gsl_complex eval_i = gsl_vector_complex_get (eval, count);
			gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, count);
			//lambda[count] = eval_i;
			lambda[count] = GSL_REAL(eval_i);
			//printf ("eigenvalue = %g\n", sqrt(eval_i));
			//printf ("eigenvector = \n");
			//gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
		}
	}
	gsl_vector_complex_free (eval);
	gsl_matrix_complex_free (evec);
	gsl_permutation_free (perm);
	gsl_matrix_free (matrix);
	gsl_matrix_free (inverse);
	for (i = 0; i < incog; i++) {
		if(lambda[i] < 0)
			lambda[i] = -lambda[i];
		lambda[i] = (sqrt(lambda[i]))/ (2 * 3.1416);
	}
    */
	// -------------------------------------------------------------------------------------------------------------------------------------------------



	//printf("\nU...\n");
    //for (i = 0; i < 18; i++) {
	//		printf ("%f\t\t",u[i]);
	//		printf("\n");
    //}

    for(i = 0; i < n; i++) {
        mxFree(AAA[i]);
	}
	mxFree(AAA);
	mxFree(x);

    return det;
}

/**
 * Tension_stress function
 */
void tension_stress(int n, int numNodes, double * u, double ** C, double * elasticity, double * length, int ** vector, double * F, double * sumConst, double * area, int **conectividade, double * aalfa){

	int i,j;
	double cx,cy,cz,kk,mm,normal,ea;


	for (i = 0; i < n; i++) {

		sumConst[i] = 0.0;

		if (area[i] != 0) {


			cx = (C[conectividade[i][0]-1][0] - C[conectividade[i][1]-1][0]) / length[i];
	        cy = (C[conectividade[i][0]-1][1] - C[conectividade[i][1]-1][1]) / length[i];
	        cz = (C[conectividade[i][0]-1][2] - C[conectividade[i][1]-1][2]) / length[i];
					mm = cos(aalfa[i]);
	        kk = sin(aalfa[i]);
	        ea = elasticity[i] * area[i] / length[i];

	        if (cx==0 && cz==0) {


	        	sumConst[i] = cy*ea*u[3*conectividade[i][0]-2] - cy*ea*u[3*conectividade[i][1]-2];

	        /*
	        	R[0][0]= 0.0;
	        	R[0][1]= cy;
	        	R[0][2]= 0.0;
	        	R[0][3]= 0.0;
	        	R[0][4]= 0.0;
	        	R[0][5]= 0.0;

	        	R[1][0]= -mm*cy;
	        	R[1][1]= 0.0;
	        	R[1][2]= kk;
	        	R[1][3]= 0.0;
	        	R[1][4]= 0.0;
	        	R[1][5]= 0.0;

	        	R[2][0]= kk*cy;
	        	R[2][1]= 0.0;
	        	R[2][2]= mm;
	        	R[2][3]= 0.0;
	        	R[2][4]= 0.0;
	        	R[2][5]= 0.0;

	        	R[3][0]= 0.0;
	        	R[3][1]= 0.0;
	        	R[3][2]= 0.0;
	        	R[3][3]= 0.0;
	        	R[3][4]= cy;
	        	R[3][5]= 0.0;

	        	R[4][0]= 0.0;
	        	R[4][1]= 0.0;
	        	R[4][2]= 0.0;
	        	R[4][3]= -mm*cy;
	        	R[4][4]= 0.0;
	        	R[4][5]= kk;

	        	R[5][0]= 0.0;
	        	R[5][1]= 0.0;
	        	R[5][2]= 0.0;
	        	R[5][3]= kk*cy;
	        	R[5][4]= 0.0;
	        	R[5][5]= mm;*/

	        }else{

	        	sumConst[i] = cx*ea*u[3*conectividade[i][0]-3] + cy*ea*u[3*conectividade[i][0]-2] + cz*ea*u[3*conectividade[i][0]-1] - (cx*ea*u[3*conectividade[i][1]-3] + cy*ea*u[3*conectividade[i][1]-2] + cz*ea*u[3*conectividade[i][1]-1]);

	        	/*R[0][0]= cx;
	        	R[0][1]= cy;
	        	R[0][2]= cz;
	        	R[0][3]= 0.0;
	        	R[0][4]= 0.0;
	        	R[0][5]= 0.0;

	        	R[1][0]= (-cx*cy*mm-cz*k)/(sqrt(cx*cx+cz*cz));
	        	R[1][1]= sqrt(cx*cx+cz*cz)*mm;
	        	R[1][2]= (-cy*cz*mm+cx*kk)/(sqrt(cx*cx+cz*cz));
	        	R[1][3]= 0.0;
	        	R[1][4]= 0.0;
	        	R[1][5]= 0.0;

	        	R[2][0]= (cx*cy*kk-cz*mm)/(sqrt(cx*cx+cz*cz));
	        	R[2][1]= -sqrt(cx*cx+cz*cz)*kk;
	        	R[2][2]= (cy*cz*kk+cx*mm)/(sqrt(cx*cx+cz*cz));
	        	R[2][3]= 0.0;
	        	R[2][4]= 0.0;
	        	R[2][5]= 0.0;

	        	R[3][0]= 0.0;
	        	R[3][1]= 0.0;
	        	R[3][2]= 0.0;
	        	R[3][3]= cx;
	        	R[3][4]= cy;
	        	R[3][5]= cz;

	        	R[4][0]= 0.0;
	        	R[4][1]= 0.0;
	        	R[4][2]= 0.0;
	        	R[4][3]= (-cx*cy*mm-cz*k)/(sqrt(cx*cx+cz*cz));
	        	R[4][4]= sqrt(cx*cx+cz*cz)*mm;
	        	R[4][5]= (-cy*cz*mm+cx*kk)/(sqrt(cx*cx+cz*cz));

	        	R[5][0]= 0.0;
	        	R[5][1]= 0.0;
	        	R[5][2]= 0.0;
	        	R[5][3]= (cx*cy*kk-cz*mm)/(sqrt(cx*cx+cz*cz));
	        	R[5][4]= -sqrt(cx*cx+cz*cz)*kk;
	        	R[5][5]= (cy*cz*kk+cx*mm)/(sqrt(cx*cx+cz*cz));*/

			}



		//for(j = 0; j < D; j++)
			//sumConst[i] += coefXY[i][j] * u[vector[i][j]] - coefXY[i][j] * u[vector[i][j+D]];
		sumConst[i] = sumConst[i] / area[i];
		//F[i] = sumConst[i] * area[i];
		}
	}

    /*
	printf ("\nNode\tX-DISPLACEMENT\tY-DISPLACEMENT\tZ-DISPLACEMENT\n");
	for (i=0; i<(numNodes); i++) {
		printf("\n %d\t", i);
		for(j=0; j<D; j++)
			printf ("%lf\t",u[D*i+j]);
	}

	printf ("\n\nElement\t\tForce\t\tStress\nNumber\n");
	for (i=0; i<n; i++)
		printf (" %d\t\t%e\t%e\n",i,F[i],sumConst[i]);



 	 printf ("\nElemento_2\tcx\t\tcy\t\tcz\n");
	 for (i=0; i<n; i++) {
	 printf ("\n%d\t\t",i);
	 for (j = 0; j < D; j++)
	 printf ("%lf\t",coefXY[i][j]);
			    }

	printf("\nVetor area\n");
	for (i=0; i<n; i++)
		printf("%f\n",area[i]);
	printf("\nVetor elasticity\n");
	for (i=0; i<n; i++)
		printf("%f\n",elasticity[i]);

	printf ("\n\n");

	FILE * out1 = fopen("out1.txt", "w");
	if (out1 == NULL) {
        printf("\nError 'out1'..\n");
        exit(1);
	}
	FILE * out2 = fopen("out2.txt", "w");
	if (out2 == NULL) {
        printf("\nError 'out2'..\n");
        exit(1);
	}


	for (i=0; i<(numNodes); i++) {
		fprintf(out1, "\n\t");
		for(j=0; j<3; j++)
			fprintf (out1, "%lf\t",u[3*i+j]);
	}

	for (i=0; i<n; i++)
		fprintf (out2, " \t%e\n",sumConst[i]);

*/
    //exit(0);

}

/**
 * Calculations function
 */
void calculations (int n, int ** conectividade, int ** vector, int numNodes, int numConst, double * nodesL, int * in, double * fr, double ** C, double * length, double * area, int * nosSaem, double * nodesC) {

	int i,j,cont=0;


	for (i = 0; i < numNodes; i++) {
		nosSaem[i]=0;
	}

	// Calculo dos comprimentos
	double c1;
	int cn = 0;
	for (i = 0; i < n; i++) {
		c1 = 0;
		for(j = 0; j < D; j++) {
			c1 += pow((C[conectividade[i][0]-1][j] - C[conectividade[i][1]-1][j]), 2);
		}
		length[i] = pow(c1,0.5);
		//printf("length[%d] = %lf;\n", i, length[i]);
	}
	//exit(0);

	// Possï¿½vel retirada de nï¿½s da estrutura devido ï¿½ topologia

	for (i = 0; i < numNodes; i++) {
		int flag=0;
		for (j=0; j<n; j++){
			//mexPrintf("area %d = %f\n",j+1,area[j]);
// 			if (area[j]!=0) {
				if ( (conectividade[j][0] == i+1) || (conectividade[j][1] == i+1)){
					//mexPrintf("Conectividade = %d %d\n",conectividade[j][0],conectividade[j][1]);
					flag = 1;
					break;
				}

// 			}
		}
		if (flag==0){
			nosSaem[cont]=i+1;
			
			cont++;
		}

	}

	/*mexPrintf("numConst: %d\n",numConst);
	mexPrintf("Nos que saem:\n");
	for (i = 0; i < numNodes; i++) {
		mexPrintf("%d ",nosSaem[i]);
	}
	mexPrintf("\n");*/




	for (i = 0; i < n; i++) {
		for(j = 0; j < D; j++) {
			vector[i][j] = conectividade[i][0] * D - (D- j);
			vector[i][j+D] = conectividade[i][1] * D - (D - j);
		}
	}

	/*for (i = 0; i < n; i++) {
		for(j = 0; j < 2*D; j++) {
			mexPrintf("%d ",vector[i][j]);
		}
		mexPrintf("\n");
	}*/


	//mexPrintf("\nNodesc\n");
	for (i = 0; i < D*numNodes; i++) {
		if(nodesC[i] == 1) {
			cn++;
			in[i] = -17;
		} else {
			in[i] = i - cn;
		}
		//mexPrintf("%f ",nodesC[i]);
	}

	//mexPrintf("\n");
    //mexPrintf("in \n");
	for (i = 0; i < D*numNodes; i++) {
		if(in[i] == -17) {
			in[i] = D * numNodes - cn;
		}
    //mexPrintf("%d ",in[i]);
	}
    //mexPrintf("\n");
	//mexPrintf("fim calculations\n");

	// Vector fr
	for (i = 0; i < D*numNodes; i++) {

		if(in[i] < (numNodes*D-numConst)) {
			fr[in[i]] = nodesL[i];
		}
	}
    
     //mexPrintf("numConst == %d\n",numConst);
    /*mexPrintf("FR \n");
	for (i = 0; i < D*numNodes-numConst; i++) {
		mexPrintf("%f\n",fr[i]);
	}
    mexPrintf("\n");*/

/*printf("\nVector...\n");
    for (i = 0; i < 10; i++) {
		for(j = 0; j < 6; j++) {
			printf ("%d\t",vector[i][j]);

		}printf("\n");
	}
printf("\nFr...\n");
	  for (i = 0; i < 18; i++) {
			printf ("%lf\t\n",fr[i]);
	}
printf("\nin...\n");
	  for (i = 0; i < 18; i++) {
			printf ("%d\t\n",in[i]);
	}*/
}


/*double strstate(double ** kr, double * u, double * fr, int n)
{
	int i,j;
	double eqValue=0.0;
	double aux=0.0;


	/*for (i=0;i<n;i++){

			//aux[i]=0.0;
			//aux[i][1]=kr[i][j];
			//mexPrintf("%f ",kr[i][j]);

		mexPrintf("%f \n",u[i]);
	}

	for (i=0;i<n;i++){
		aux=0.0;
		for (j=0;j<n;j++) {
			aux+=kr[i][j]*u[j];													//Ku
			//mexPrintf("kr u kr*u = %f %f %f\n",kr[i][j],u[j],aux);
			//mexPrintf("%f ",kr[i][j]);
		}
		//mexPrintf("\n");
		aux-=fr[i];							//Ku - F
		//mexPrintf("linha de Ku-F = %f \n",aux);
		eqValue+=aux*aux;

		//mexPrintf("\n");
	}

	eqValue=sqrt(eqValue);				//norma de Ku - F

	aux=0.0;

	for (i=0;i<n;i++){
		aux+=fr[i]*fr[i];			//norma de F
	}

	aux=sqrt(aux);

	//mexPrintf("Norma de KU - F = %f \n",eqValue);
	//mexPrintf("Norma de F = %f \n",aux);

	eqValue/=aux;

	//mexPrintf("Norma de (KU - F)/F = %f \n",eqValue);



	if (eqValue>1){
		//mexPrintf("%f \n",eqValue);
		return(eqValue);
	}else{
		//mexPrintf("Not a number \n");
		return(9999);
	}




}*/

void compliance (double **kr, double *u, double * nodesC,int numNodes, double * comp, double *fr) {
	int i,j;
	int numConst=0;
	double temp3=0.0;
	*comp = 0.0;
	for (i=0;i<Dkl*numNodes;i++){
		if (nodesC[i]==0){
			numConst+=1;
		}
	}
	double temp2[numConst];


	for (i=0;i<numConst;i++){
		for (j=0;j<numConst;j++){
			temp2[i]=0.0;
			temp2[i]+=kr[i][j]*u[j];
			//mexPrintf("%f ",kr[i][j]);
		}
		//mexPrintf("\n");
	}

	//mexPrintf("\n\n\n\n");

	for (i=0;i<numConst;i++){
		*comp+=temp2[i]*u[i];
		//mexPrintf("%f ",u[i]);
	}

	*comp*=0.5;

	for (i=0;i<numConst;i++){
		temp3+=fr[i]*u[i];
	}

	*comp-=temp3;

	temp3=0.0;
	for (i=0;i<numConst;i++){
		temp3+=u[i]*u[i];
	}
	temp3*=10^-12;

	*comp+=temp3;

	//mexPrintf("Compliance = %f",temp);
	//return(temp);

}



void ccrit(int n, double * sumConst, double * length, int ** vector, int numNodes, int numConst, double * Kg, double ** kg, int * in, double * area, double * elasticity, int ** conectividade, double ** C, double * aalfa){

	int cont=0,i=0,j=0;
	double P=0.0,a,l,cx,cy,cz,ea,kk,mm;

	for (i = 0; i < (numNodes*D-numConst); i++) {
		for (j = 0; j < (numNodes*D-numConst); j++) {		    
		    kg[i][j] = 0.;
		}
	}

	for(i = 0; i < n; i++) {
		P=sumConst[i]*area[i];
		//mexPrintf("Elemento %d, Tensao = %f, lenght = %f, Normal = %f \n",i+1,sumConst[i],length[i],P);     
        
        a = area[i];
        l = length[i];
        cx = (C[conectividade[i][1]-1][0] - C[conectividade[i][0]-1][0]) / length[i];
        cy = (C[conectividade[i][1]-1][1] - C[conectividade[i][0]-1][1]) / length[i];
        cz = (C[conectividade[i][1]-1][2] - C[conectividade[i][0]-1][2]) / length[i];
        mm = cos(aalfa[i]);
        kk = sin(aalfa[i]);
        
        //mexPrintf("EA = %f\n",ea);
        
                //mexPrintf("sen cos = %f %f\n",kk,mm);
                //mexPrintf("ro = %f\n",ro);
        //mm=1.0;
        //kk=0.0;
        ea = elasticity[i] * area[i] / length[i];




        if ((cx == 0.0) && (cz == 0.0)) {

            //mexPrintf("Entrou aqui (barra vertical)\n");

                       
            //matrizes pra elementos 2-2 , 2-5, 5-2 e 5-5 iguais a 1
            /*kg[in[vector[i][0]]][in[vector[i][0]]] += ((P*cy*cy*mm*mm)/(l)) ;
            kg[in[vector[i][0]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][0]]][in[vector[i][2]]] += ((-P*cy*kk*mm)/(l)) ;
            kg[in[vector[i][0]]][in[vector[i][3]]] += ((-P*cy*cy*mm*mm)/(l));
            kg[in[vector[i][0]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][0]]][in[vector[i][5]]] += ((P*cy*kk*mm)/(l));;

            kg[in[vector[i][1]]][in[vector[i][0]]] += 0.0;
            kg[in[vector[i][1]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][1]]][in[vector[i][2]]] += 0.0;
            kg[in[vector[i][1]]][in[vector[i][3]]] += 0.0;
            kg[in[vector[i][1]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][1]]][in[vector[i][5]]] += 0.0;
            
            kg[in[vector[i][2]]][in[vector[i][0]]] += ((-P*cy*kk*mm)/(l)) ;
            kg[in[vector[i][2]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][2]]][in[vector[i][2]]] += ((P*kk*kk)/(l)) ;
            kg[in[vector[i][2]]][in[vector[i][3]]] += ((P*cy*kk*mm)/(l)) ;
            kg[in[vector[i][2]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][2]]][in[vector[i][5]]] += ((-P*kk*kk)/(l));;
            
            kg[in[vector[i][3]]][in[vector[i][0]]] += ((-P*cy*cy*mm*mm)/(l)) ;
            kg[in[vector[i][3]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][3]]][in[vector[i][2]]] += ((P*cy*kk*mm)/(l)) ;
            kg[in[vector[i][3]]][in[vector[i][3]]] += ((P*cy*cy*mm*mm)/(l)) ;
            kg[in[vector[i][3]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][3]]][in[vector[i][5]]] += ((-P*cy*kk*mm)/(l));
            
            kg[in[vector[i][4]]][in[vector[i][0]]] += 0.0;
            kg[in[vector[i][4]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][4]]][in[vector[i][2]]] += 0.0;
            kg[in[vector[i][4]]][in[vector[i][3]]] += 0.0;
            kg[in[vector[i][4]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][4]]][in[vector[i][5]]] += 0.0;
            
            kg[in[vector[i][5]]][in[vector[i][0]]] += ((P*cy*kk*mm)/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][5]]][in[vector[i][2]]] += ((-P*kk*kk)/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][3]]] += ((-P*cy*kk*mm)/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][5]]][in[vector[i][5]]] += ((P*kk*kk)/(l));*/
            
            //matrizes abaixo para diagonal full 1
            
            kg[in[vector[i][0]]][in[vector[i][0]]] += P*cy*cy*(kk*kk+mm*mm)/l;
            kg[in[vector[i][0]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][0]]][in[vector[i][2]]] += 0.0;
            kg[in[vector[i][0]]][in[vector[i][3]]] += P*cy*cy*(-kk*kk-mm*mm)/l;
            kg[in[vector[i][0]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][0]]][in[vector[i][5]]] += 0.0;

            kg[in[vector[i][1]]][in[vector[i][0]]] += 0.0;
            kg[in[vector[i][1]]][in[vector[i][1]]] += P*cy*cy/l;
            kg[in[vector[i][1]]][in[vector[i][2]]] += 0.0;
            kg[in[vector[i][1]]][in[vector[i][3]]] += 0.0;
            kg[in[vector[i][1]]][in[vector[i][4]]] += -P*cy*cy/l;
            kg[in[vector[i][1]]][in[vector[i][5]]] += 0.0;

            kg[in[vector[i][2]]][in[vector[i][0]]] += 0.0;
            kg[in[vector[i][2]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][2]]][in[vector[i][2]]] += P*(kk*kk+mm*mm)/l;
            kg[in[vector[i][2]]][in[vector[i][3]]] += 0.0;
            kg[in[vector[i][2]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][2]]][in[vector[i][5]]] += P*(-kk*kk-mm*mm)/l;

            kg[in[vector[i][3]]][in[vector[i][0]]] += P*cy*cy*(-kk*kk-mm*mm)/l;
            kg[in[vector[i][3]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][3]]][in[vector[i][2]]] += 0.0;
            kg[in[vector[i][3]]][in[vector[i][3]]] += P*cy*cy*(kk*kk+mm*mm)/l;
            kg[in[vector[i][3]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][3]]][in[vector[i][5]]] += 0.0;

            kg[in[vector[i][4]]][in[vector[i][0]]] += 0.0;
            kg[in[vector[i][4]]][in[vector[i][1]]] += -P*cy*cy/l;
            kg[in[vector[i][4]]][in[vector[i][2]]] += 0.0;
            kg[in[vector[i][4]]][in[vector[i][3]]] += 0.0;
            kg[in[vector[i][4]]][in[vector[i][4]]] += P*cy*cy/l;
            kg[in[vector[i][4]]][in[vector[i][5]]] += 0.0;

            kg[in[vector[i][5]]][in[vector[i][0]]] += 0.0;
            kg[in[vector[i][5]]][in[vector[i][1]]] += 0.0;
            kg[in[vector[i][5]]][in[vector[i][2]]] += P*(-kk*kk-mm*mm)/l;
            kg[in[vector[i][5]]][in[vector[i][3]]] += 0.0;
            kg[in[vector[i][5]]][in[vector[i][4]]] += 0.0;
            kg[in[vector[i][5]]][in[vector[i][5]]] += P*(kk*kk+mm*mm)/l;
            
            

        } else {
            
            //mexPrintf("Cx Cy Cz = %f %f %f\n mm = %f kk = %f\n",cx,cy,cz,mm,kk);
            //mexPrintf("Entrou aqui (barra nao vertical)    %f \n",mm*cx*cy+kk*cz);

           
            //matrizes pra elementos 2-2 , 2-5, 5-2 e 5-5 iguais a 1
            /*kg[in[vector[i][0]]][in[vector[i][0]]] += ((P*(mm*cx*cy+kk*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz)));
            kg[in[vector[i][0]]][in[vector[i][1]]] += ((-P*mm*(mm*cx*cy+kk*cz))/(l));
            kg[in[vector[i][0]]][in[vector[i][2]]] += ((-P*(kk*cx-mm*cy*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz)));
            kg[in[vector[i][0]]][in[vector[i][3]]] += ((-P*(mm*cx*cy+kk*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz))) ;
            kg[in[vector[i][0]]][in[vector[i][4]]] += ((P*mm*(mm*cx*cy+kk*cz))/(l));
            kg[in[vector[i][0]]][in[vector[i][5]]] += ((P*(kk*cx-mm*cy*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz)));

            kg[in[vector[i][1]]][in[vector[i][0]]] += ((-P*mm*(mm*cx*cy+kk*cz))/(l));
            kg[in[vector[i][1]]][in[vector[i][1]]] += ((P*mm*mm*(cx*cx+cz*cz))/(l));
            kg[in[vector[i][1]]][in[vector[i][2]]] += ((P*mm*(kk*cx-mm*cy*cz))/(l));
            kg[in[vector[i][1]]][in[vector[i][3]]] += ((P*mm*(mm*cx*cy+kk*cz))/(l));
            kg[in[vector[i][1]]][in[vector[i][4]]] += ((-P*mm*mm*(cx*cx+cz*cz))/(l));
            kg[in[vector[i][1]]][in[vector[i][5]]] += ((-P*mm*(kk*cx-mm*cy*cz))/(l));
            
            kg[in[vector[i][2]]][in[vector[i][0]]] += ((-P*(kk*cx-mm*cy*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz)));
            kg[in[vector[i][2]]][in[vector[i][1]]] += ((P*mm*(kk*cx-mm*cy*cz))/(l));
            kg[in[vector[i][2]]][in[vector[i][2]]] += ((P*(kk*cx-mm*cy*cz)*(kk*cx-mm*cy*cz))/(l*(cx*cx+cz*cz)));
            kg[in[vector[i][2]]][in[vector[i][3]]] += ((P*(kk*cx-mm*cy*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz)));
            kg[in[vector[i][2]]][in[vector[i][4]]] += ((-P*mm*(kk*cx-mm*cy*cz))/(l));
            kg[in[vector[i][2]]][in[vector[i][5]]] += ((-P*(kk*cx-mm*cy*cz)*(kk*cx-mm*cy*cz))/(l*(cx*cx+cz*cz)));
            
            kg[in[vector[i][3]]][in[vector[i][0]]] += ((-P*(mm*cx*cy+kk*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz)));
            kg[in[vector[i][3]]][in[vector[i][1]]] += ((P*mm*(mm*cx*cy+kk*cz))/(l));
            kg[in[vector[i][3]]][in[vector[i][2]]] += ((P*(kk*cx-mm*cy*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz)));
            kg[in[vector[i][3]]][in[vector[i][3]]] += ((P*(mm*cx*cy+kk*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz))); 
            kg[in[vector[i][3]]][in[vector[i][4]]] += ((-P*mm*(mm*cx*cy+kk*cz))/(l));
            kg[in[vector[i][3]]][in[vector[i][5]]] += ((-P*(kk*cx-mm*cy*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz)));
            
            kg[in[vector[i][4]]][in[vector[i][0]]] += ((P*mm*(mm*cx*cy+kk*cz))/(l)) ;
            kg[in[vector[i][4]]][in[vector[i][1]]] += ((-P*mm*mm*(cx*cx+cz*cz))/(l)) ;
            kg[in[vector[i][4]]][in[vector[i][2]]] += ((-P*mm*(kk*cx-mm*cy*cz))/(l)) ;
            kg[in[vector[i][4]]][in[vector[i][3]]] += ((-P*mm*(mm*cx*cy+kk*cz))/(l)) ;
            kg[in[vector[i][4]]][in[vector[i][4]]] += ((P*mm*mm*(cx*cx+cz*cz))/(l)) ;
            kg[in[vector[i][4]]][in[vector[i][5]]] += ((P*mm*(kk*cx-mm*cy*cz))/(l));    
            
            kg[in[vector[i][5]]][in[vector[i][0]]] += ((P*(kk*cx-mm*cy*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz))) ;
            kg[in[vector[i][5]]][in[vector[i][1]]] += ((-P*mm*(kk*cx-mm*cy*cz))/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][2]]] += ((-P*(kk*cx-mm*cy*cz)*(kk*cx-mm*cy*cz))/(l*(cx*cx+cz*cz))) ;
            kg[in[vector[i][5]]][in[vector[i][3]]] += ((-P*(kk*cx-mm*cy*cz)*(mm*cx*cy+kk*cz))/(l*(cx*cx+cz*cz))) ;
            kg[in[vector[i][5]]][in[vector[i][4]]] += ((P*mm*(kk*cx-mm*cy*cz))/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][5]]] += ((P*(kk*cx-mm*cy*cz)*(kk*cx-mm*cy*cz))/(l*(cx*cx+cz*cz)));*/
            
            
            
            //matrizes abaixo para diagonal full 1            
            kg[in[vector[i][0]]][in[vector[i][0]]] += ((-P*(kk*kk+mm*mm)*(cy*cy-1)*cz*cz)/(l*(cx*cx+cz*cz)))+((P*cx*cx)/(l))+((P*(kk*kk+mm*mm)*cy*cy)/(l)) ;
            kg[in[vector[i][0]]][in[vector[i][1]]] += ((P*(-kk*kk-mm*mm+1)*cx*cy)/(l));
            kg[in[vector[i][0]]][in[vector[i][2]]] += ((P*(kk*kk+mm*mm)*cx*(cy*cy-1)*cz)/(l*(cx*cx+cz*cz)))+((P*cx*cz)/(l));
            kg[in[vector[i][0]]][in[vector[i][3]]] += ((P*(kk*kk+mm*mm)*(cy*cy-1)*cz*cz)/(l*(cx*cx+cz*cz)))-((P*cx*cx)/(l))-((P*(kk*kk+mm*mm)*cy*cy)/(l)) ;
            kg[in[vector[i][0]]][in[vector[i][4]]] += ((P*(kk*kk+mm*mm-1)*cx*cy)/(l)) ;
            kg[in[vector[i][0]]][in[vector[i][5]]] += ((-P*(kk*kk+mm*mm)*cx*(cy*cy-1)*cz)/(l*(cx*cx+cz*cz)))-((P*cx*cz)/(l));

            kg[in[vector[i][1]]][in[vector[i][0]]] += ((P*(-kk*kk-mm*mm+1)*cx*cy)/(l)) ;
            kg[in[vector[i][1]]][in[vector[i][1]]] += ((P*((kk*kk+mm*mm)*cx*cx+cy*cy+(kk*kk+mm*mm)*cz*cz))/(l)) ;
            kg[in[vector[i][1]]][in[vector[i][2]]] += ((P*(-kk*kk-mm*mm+1)*cy*cz)/(l)) ;
            kg[in[vector[i][1]]][in[vector[i][3]]] += ((P*(kk*kk+mm*mm-1)*cx*cy)/(l)) ;
            kg[in[vector[i][1]]][in[vector[i][4]]] += ((-P*((kk*kk+mm*mm)*cx*cx+cy*cy+(kk*kk+mm*mm)*cz*cz))/(l)) ;
            kg[in[vector[i][1]]][in[vector[i][5]]] += ((P*(kk*kk+mm*mm-1)*cy*cz)/(l));

            kg[in[vector[i][2]]][in[vector[i][0]]] += ((P*(kk*kk+mm*mm)*cx*(cy*cy-1)*cz)/(l*(cx*cx+cz*cz)))+((P*cx*cz)/(l));
            kg[in[vector[i][2]]][in[vector[i][1]]] += ((P*(-kk*kk-mm*mm+1)*cy*cz)/(l)) ;
            kg[in[vector[i][2]]][in[vector[i][2]]] += ((P*(kk*kk+mm*mm)*(cy*cy-1)*cz*cz)/(l*(cx*cx+cz*cz)))+((P*cz*cz)/(l))+((P*(kk*kk+mm*mm))/(l)) ;
            kg[in[vector[i][2]]][in[vector[i][3]]] += ((-P*(kk*kk+mm*mm)*cx*(cy*cy-1)*cz)/(l*(cx*cx+cz*cz)))-((P*cx*cz)/(l)) ;
            kg[in[vector[i][2]]][in[vector[i][4]]] += ((P*(kk*kk+mm*mm-1)*cy*cz)/(l)) ;
            kg[in[vector[i][2]]][in[vector[i][5]]] += ((-P*(kk*kk+mm*mm)*(cy*cy-1)*cz*cz)/(l*(cx*cx+cz*cz)))-((P*cz*cz)/(l))-((P*(kk*kk+mm*mm))/(l));;

            kg[in[vector[i][3]]][in[vector[i][0]]] += ((P*(kk*kk+mm*mm)*(cy*cy-1)*cz*cz)/(l*(cx*cx+cz*cz)))-((P*cx*cx)/(l))-((P*(kk*kk+mm*mm)*cy*cy)/(l)) ;
            kg[in[vector[i][3]]][in[vector[i][1]]] += ((P*(kk*kk+mm*mm-1)*cx*cy)/(l)) ;
            kg[in[vector[i][3]]][in[vector[i][2]]] += ((-P*(kk*kk+mm*mm)*cx*(cy*cy-1)*cz)/(l*(cx*cx+cz*cz)))-((P*cx*cz)/(l)) ;
            kg[in[vector[i][3]]][in[vector[i][3]]] += ((-P*(kk*kk+mm*mm)*(cy*cy-1)*cz*cz)/(l*(cx*cx+cz*cz)))+((P*cx*cx)/(l))+((P*(kk*kk+mm*mm)*cy*cy)/(l)) ;
            kg[in[vector[i][3]]][in[vector[i][4]]] += ((P*(-kk*kk-mm*mm+1)*cx*cy)/(l));
            kg[in[vector[i][3]]][in[vector[i][5]]] += ((P*(kk*kk+mm*mm)*cx*(cy*cy-1)*cz)/(l*(cx*cx+cz*cz)))+((P*cx*cz)/(l));;

            kg[in[vector[i][4]]][in[vector[i][0]]] += ((P*(kk*kk+mm*mm-1)*cx*cy)/(l));
            kg[in[vector[i][4]]][in[vector[i][1]]] += ((-P*((kk*kk+mm*mm)*cx*cx+cy*cy+(kk*kk+mm*mm)*cz*cz))/(l)) ;
            kg[in[vector[i][4]]][in[vector[i][2]]] += ((P*(kk*kk+mm*mm-1)*cy*cz)/(l));
            kg[in[vector[i][4]]][in[vector[i][3]]] += ((P*(-kk*kk-mm*mm+1)*cx*cy)/(l))  ;
            kg[in[vector[i][4]]][in[vector[i][4]]] += ((P*((kk*kk+mm*mm)*cx*cx+cy*cy+(kk*kk+mm*mm)*cz*cz))/(l)) ;
            kg[in[vector[i][4]]][in[vector[i][5]]] += ((P*(-kk*kk-mm*mm+1)*cy*cz)/(l));

            kg[in[vector[i][5]]][in[vector[i][0]]] += ((-P*(kk*kk+mm*mm)*cx*(cy*cy-1)*cz)/(l*(cx*cx+cz*cz)))-((P*cx*cz)/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][1]]] += ((P*(kk*kk+mm*mm-1)*cy*cz)/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][2]]] += ((-P*(kk*kk+mm*mm)*(cy*cy-1)*cz*cz)/(l*(cx*cx+cz*cz)))-((P*cz*cz)/(l))-((P*(kk*kk+mm*mm))/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][3]]] += ((P*(kk*kk+mm*mm)*cx*(cy*cy-1)*cz)/(l*(cx*cx+cz*cz)))+((P*cx*cz)/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][4]]] += ((P*(-kk*kk-mm*mm+1)*cy*cz)/(l)) ;
            kg[in[vector[i][5]]][in[vector[i][5]]] += ((P*(kk*kk+mm*mm)*(cy*cy-1)*cz*cz)/(l*(cx*cx+cz*cz)))+((P*cz*cz)/(l))+((P*(kk*kk+mm*mm))/(l));
            
            
            
        }



	}



	for (i=0;i<numNodes*Dkl-numConst;i++){
		for (j=0;j<numNodes*Dkl-numConst;j++){
			Kg[cont]=0.0;
			//Mmat[cont]=0;
			Kg[cont]=kg[j][i];
			//Mmat[cont]=mr[j][i];
			cont++;
		}
	}


}
