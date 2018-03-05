#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <string.h>


/*
 * Important Notice: September 20, 2010
 *
 * In this head file, we assume that the features in the tree strucutre
 * are well ordered. That is to say, the indices of the left nodes is always less
 * than the right nodes. Ideally, this can be achieved by reordering the features.
 *
 * The advantage of this ordered features is that, we donot need to use an explicit
 * variable for recording the indices.
 *
 * To deal with the more general case when the features might not be well ordered,
 * we provide the functions in the head file "general_altra.h". Compared with the files in this head file,
 * we need an additional parameter G, which contains the indices of the nodes.
 *
 *
 */

/*
 * -------------------------------------------------------------------
 *                       Functions and parameter
 * -------------------------------------------------------------------
 *
 * altra solves the following problem
 *
 * 1/2 \|x-v\|^2 + \sum \lambda_i \|x_{G_i}\|,
 *
 * where x and v are of dimension n,
 *       \lambda_i >=0, and G_i's follow the tree structure
 *
 * It is implemented in Matlab as follows:
 *
 * x=altra(v, n, ind, nodes);
 *
 * ind is a 3 x nodes matrix.
 *       Each column corresponds to a node.
 *
 *       The first element of each column is the starting index,
 *       the second element of each column is the ending index
 *       the third element of each column corrreponds to \lambbda_i.
 *
 * -------------------------------------------------------------------
 *                       Notices:
 * -------------------------------------------------------------------
 *
 * 1. The nodes in the parameter "ind" should be given in the 
 *    either
 *           the postordering of depth-first traversal
 *    or 
 *           the reverse breadth-first traversal.
 *
 * 2. When each elements of x are penalized via the same L1 
 *    (equivalent to the L2 norm) parameter, one can simplify the input
 *    by specifying 
 *           the "first" column of ind as (-1, -1, lambda)
 *
 *    In this case, we treat it as a single "super" node. Thus in the value
 *    nodes, we only count it once.
 *
 * 3. The values in "ind" are in [1,n].
 *
 * 4. The third element of each column should be positive. The program does
 *    not check the validity of the parameter. 
 *
 *    It is still valid to use the zero regularization parameter.
 *    In this case, the program does not change the values of 
 *    correponding indices.
 *    
 *
 * -------------------------------------------------------------------
 *                       History:
 * -------------------------------------------------------------------
 *
 * Composed by Jun Liu on April 20, 2010
 *
 * For any question or suggestion, please email j.liu@asu.edu.
 *
 */


void altra(double *x, double *v, int n, double *ind, int nodes){
    
    int i, j, m;
    double lambda,twoNorm, ratio;
    
    /*
     * test whether the first node is special
     */
    if ((int) ind[0]==-1){
        
        /*
         *Recheck whether ind[1] equals to zero
         */
        if ((int) ind[1]!=-1){
            printf("\n Error! \n Check ind");
            exit(1);
        }        
        
        lambda=ind[2];
        
        for(j=0;j<n;j++){
            if (v[j]>lambda)
                x[j]=v[j]-lambda;
            else
                if (v[j]<-lambda)
                    x[j]=v[j]+lambda;
                else
                    x[j]=0;
        }
        
        i=1;
    }
    else{
        memcpy(x, v, sizeof(double) * n);
        i=0;
    }
            
    /*
     * sequentially process each node
     *
     */
	for(;i < nodes; i++){
        /*
         * compute the L2 norm of this group         
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];        
        twoNorm=sqrt(twoNorm);
        
        lambda=ind[3*i+2];
        if (twoNorm>lambda){
            ratio=(twoNorm-lambda)/twoNorm;
            
            /*
             * shrinkage this group by ratio
             */
            for(j=(int) ind[3*i]-1;j<(int) ind[3*i+1];j++)
                x[j]*=ratio;            
        }
        else{
            /*
             * threshold this group to zero
             */
            for(j=(int) ind[3*i]-1;j<(int) ind[3*i+1];j++)
                x[j]=0;
        }
	}
}



/*
 * altra_mt is a generalization of altra to the 
 * 
 * multi-task learning scenario (or equivalently the multi-class case)
 *
 * altra_mt(X, V, n, k, ind, nodes);
 *
 * It applies altra for each row (1xk) of X and V
 *
 */


void altra_mt(double *X, double *V, int n, int k, double *ind, int nodes){
    int i, j;
    
    double *x=(double *)malloc(sizeof(double)*k);
    double *v=(double *)malloc(sizeof(double)*k);
    
    for (i=0;i<n;i++){
        /*
         * copy a row of V to v
         *         
         */
        for(j=0;j<k;j++)
            v[j]=V[j*n + i];
        
        altra(x, v, k, ind, nodes);
        
        /*
         * copy the solution to X         
         */        
        for(j=0;j<k;j++)
            X[j*n+i]=x[j];
    }
    
    free(x);
    free(v);
}




/*
 * compute
 *  lambda2_max=computeLambda2Max(x,n,ind,nodes);
 *
 * compute the 2 norm of each group, which is divided by the ind(3,:),
 * then the maximum value is returned
 */
    /*
     *This function does not consider the case ind={[-1, -1, 100]',...}
     *
     *This functions is not used currently.
     */

void computeLambda2Max(double *lambda2_max, double *x, int n, double *ind, int nodes){
    int i, j, m;
    double lambda,twoNorm;
    
    *lambda2_max=0;
    
    for(i=0;i < nodes; i++){
        /*
         * compute the L2 norm of this group         
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];        
        twoNorm=sqrt(twoNorm);
        
        twoNorm=twoNorm/ind[3*i+2];
        
        if (twoNorm >*lambda2_max )
            *lambda2_max=twoNorm;        
	}
}

/*
 * -------------------------------------------------------------------
 *                       Function and parameter
 * -------------------------------------------------------------------
 *
 * treeNorm compute
 *
 *        \sum \lambda_i \|x_{G_i}\|,
 *
 * where x is of dimension n,
 *       \lambda_i >=0, and G_i's follow the tree structure
 *
 * The file is implemented in the following in Matlab:
 *
 * tree_norm=treeNorm(x, n, ind,nodes);
 */


void treeNorm(double *tree_norm, double *x, int n, double *ind, int nodes){
    
    int i, j, m;
    double twoNorm, lambda;
    
    *tree_norm=0;
    
    /*
     * test whether the first node is special
     */
    if ((int) ind[0]==-1){
        
        /*
         *Recheck whether ind[1] equals to zero
         */
        if ((int) ind[1]!=-1){
            printf("\n Error! \n Check ind");
            exit(1);
        }        
        
        lambda=ind[2];
        
        for(j=0;j<n;j++){
            *tree_norm+=fabs(x[j]);
        }
        
        *tree_norm=*tree_norm * lambda;
        
        i=1;
    }
    else{
        i=0;
    }
            
    /*
     * sequentially process each node
     *
     */
	for(;i < nodes; i++){
        /*
         * compute the L2 norm of this group         
         */
		twoNorm=0;
		for(j=(int) ind[3*i]-1;j< (int) ind[3*i+1];j++)
			twoNorm += x[j] * x[j];        
        twoNorm=sqrt(twoNorm);
        
        lambda=ind[3*i+2];
        
        *tree_norm=*tree_norm + lambda*twoNorm;
	}
}


/*
 * -------------------------------------------------------------------
 *                       Function and parameter
 * -------------------------------------------------------------------
 *
 * findLambdaMax compute
 * 
 * the lambda_{max} that achieves a zero solution for
 *
 *     min  1/2 \|x-v\|^2 +  \lambda_{\max} * \sum  w_i \|x_{G_i}\|,
 *
 * where x is of dimension n,
 *       w_i >=0, and G_i's follow the tree structure
 *
 * The file is implemented in the following in Matlab:
 *
 * lambdaMax=findLambdaMax(v, n, ind,nodes);
 */

void findLambdaMax(double *lambdaMax, double *v, int n, double *ind, int nodes){
 
    int i, j;
    double lambda=0,squaredWeight=0, lambda1,lambda2;
    double *x=(double *)malloc(sizeof(double)*n);
    double *ind2=(double *)malloc(sizeof(double)*nodes*3);
    int num=0;
       
    for(i=0;i<n;i++){
        lambda+=v[i]*v[i];
    }
    
    if ( (int)ind[0]==-1 )
        squaredWeight=n*ind[2]*ind[2];
    else
        squaredWeight=ind[2]*ind[2];
    
    for (i=1;i<nodes;i++){
        squaredWeight+=ind[3*i+2]*ind[3*i+2];
    }
    
    /* set lambda to an initial guess
     */
    lambda=sqrt(lambda/squaredWeight);
    
    /*
    printf("\n\n   lambda=%2.5f",lambda);
    */
    
    /*
     *copy ind to ind2,
     *and scale the weight 3*i+2
     */
    for(i=0;i<nodes;i++){
        ind2[3*i]=ind[3*i];
        ind2[3*i+1]=ind[3*i+1];
        ind2[3*i+2]=ind[3*i+2]*lambda;
    }
    
    /* test whether the solution is zero or not
     */
    altra(x, v, n, ind2, nodes);    
    for(i=0;i<n;i++){
        if (x[i]!=0)
            break;
    }
    
    if (i>=n) {
        /*x is a zero vector*/
        lambda2=lambda;
        lambda1=lambda;
        
        num=0;
        
        while(1){
            num++;
            
            lambda2=lambda;
            lambda1=lambda1/2;
            /* update ind2
             */
            for(i=0;i<nodes;i++){
                ind2[3*i+2]=ind[3*i+2]*lambda1;
            }
            
            /* compute and test whether x is zero
             */
            altra(x, v, n, ind2, nodes);
            for(i=0;i<n;i++){
                if (x[i]!=0)
                    break;
            }
            
            if (i<n){
                break;
                /*x is not zero
                 *we have found lambda1
                 */
            }
        }

    }
    else{
        /*x is a non-zero vector*/
        lambda2=lambda;
        lambda1=lambda;
        
        num=0;
        while(1){
            num++;            
            
            lambda1=lambda2;
            lambda2=lambda2*2;
            /* update ind2
             */
            for(i=0;i<nodes;i++){
                ind2[3*i+2]=ind[3*i+2]*lambda2;
            }
            
            /* compute and test whether x is zero
             */
            altra(x, v, n, ind2, nodes);
            for(i=0;i<n;i++){
                if (x[i]!=0)
                    break;
            }
            
            if (i>=n){
                break;
                /*x is a zero vector
                 *we have found lambda2
                 */
            }
        }
    }    
    
    /*
    printf("\n num=%d, lambda1=%2.5f, lambda2=%2.5f",num, lambda1,lambda2);
    */
    
    while ( fabs(lambda2-lambda1) > lambda2 * 1e-10 ){
        
        num++;
        
        lambda=(lambda1+lambda2)/2;
        
        /* update ind2
         */
        for(i=0;i<nodes;i++){
            ind2[3*i+2]=ind[3*i+2]*lambda;
        }
        
        /* compute and test whether x is zero
         */
        altra(x, v, n, ind2, nodes);
        for(i=0;i<n;i++){
            if (x[i]!=0)
                break;
        }
        
        if (i>=n){
            lambda2=lambda;
        }
        else{
            lambda1=lambda;
        }
        
       /*
        printf("\n lambda1=%2.5f, lambda2=%2.5f",lambda1,lambda2);
        */
    }
    
   /*
    printf("\n num=%d",num);
    
    printf("   lambda1=%2.5f, lambda2=%2.5f",lambda1,lambda2);
     
    */
    
    *lambdaMax=lambda2;
    
    free(x);
    free(ind2);
}


/*
 * findLambdaMax_mt is a generalization of findLambdaMax to the 
 * 
 * multi-task learning scenario (or equivalently the multi-class case)
 *
 * lambdaMax=findLambdaMax_mt(X, V, n, k, ind, nodes);
 *
 * It applies findLambdaMax for each row (1xk) of X and V
 *
 */


void findLambdaMax_mt(double *lambdaMax, double *V, int n, int k, double *ind, int nodes){
    int i, j;
    
    double *v=(double *)malloc(sizeof(double)*k);
    double lambda;
    
    *lambdaMax=0;
    
    for (i=0;i<n;i++){
        /*
         * copy a row of V to v
         *         
         */
        for(j=0;j<k;j++)
            v[j]=V[j*n + i];
        
        findLambdaMax(&lambda, v, k, ind, nodes);
        
        /*
        printf("\n   lambda=%5.2f",lambda);        
         */
        
        if (lambda>*lambdaMax)
            *lambdaMax=lambda;
    }
    
    /*
    printf("\n *lambdaMax=%5.2f",*lambdaMax);
     */
    
    free(v);
}

