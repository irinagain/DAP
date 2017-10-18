#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//Irina Gaynanova, July 2017

//function for maximum
double max(double a, double b){
  if (a>=b){
    return a;
  }else{
    return b;
  }
}

//function for inner product of two vectors
double mul(double *a, double *b, int* length){
  int index;
  double result=0;
  
  for (index = 0; index < *length; index++){
    result+=a[index] * b[index];
  }
  return result;
}

//function for one round of row update, here V is a p by 2 matrix
void rowUpdateProj_withS(double *X1, double *X2, double *V, double *lambda, int *p, int *n1, int *n2, double *errV, double *R1, double *R2){
  int k,j;
  double normt,v_old[2],t[2];

  *errV=0;
  for (k=0; k<*p; k++){
    //calculate t
    t[0]=mul(&X1[k*(*n1)],&R1[0],n1)/(*n1)+V[k];
    t[1]=mul(&X2[k*(*n2)],&R2[0],n2)/(*n2)+V[k+(*p)];
    
    //calculate L_2 norm of t
    normt=sqrt(pow(t[0],2)+ pow(t[1],2));
    
    // Save current values
    v_old[0]=V[k];
    v_old[1]=V[k+(*p)];
    
    //update kth row of V
    //if normt is less than lambda, V_{1i} = 0
    if (normt <= *lambda){
      // Perform the update of V
      V[k]=0;
      V[k+(*p)]=0;
      
      // Update residuals 
      for (j=0; j<*n1;j++){
        R1[j]+=v_old[0]*X1[j+k*(*n1)];
      }
      for (j=0; j<*n2;j++){
        R2[j]+=v_old[1]*X2[j+k*(*n2)];
      }
      
      //error
      *errV=max(max(*errV,fabs(v_old[0])), fabs(v_old[1])); 
    }else{
      // Perform the update of V
      V[k]=(1-*lambda/normt)*t[0];
      V[k+(*p)]=(1-*lambda/normt)*t[1];
    
      // Update residuals 
      for (j=0; j<*n1;j++){
        R1[j]+=(v_old[0]-V[k])*X1[j+k*(*n1)];
      }
      for (j=0; j<*n2;j++){
        R2[j]+=(v_old[1]-V[k+(*p)])*X2[j+k*(*n2)];
      }
    
      //error
      *errV=max(max(*errV,fabs(v_old[0] - V[k])), fabs(v_old[1] - V[k+(*p)])); 
    }
  }
}

//function for initializing R1 and R2, two residual vectors, input V is also a p by 2 matrix
void intlzR12(double *X1, double *X2, double *V, int *p, int *n1, int *n2, double *R1, double *R2){
  
  int j,k;
  
  for (j=0; j<*n1;j++){
    R1[j]=1;
    for (k=0; k<*p; ++k){
      R1[j]-=X1[j+k*(*n1)]*V[k];
    }
  }
  
  for (j=0; j<*n2;j++){
    R2[j]=-1;
    for (k=0; k<*p; ++k){
      R2[j]-=X2[j+k*(*n2)]*V[k+(*p)];
    }
  }
}


//Complete solver for standardized projections with lasso
void solveProj_withS(double *X1, double *X2, double *V, double *lambda, int *p, int *n1, int *n2, double *eps, int *maxiter, int *niter){

  double errV;
  double R1[(*n1)], R2[(*n2)]; 

  //calculate R1 and R2, residuals
  intlzR12(X1,X2,V,p,n1,n2,&R1[0],&R2[0]);

  *niter=0;
  do{
    ++*niter;
    rowUpdateProj_withS(X1,X2,V,lambda,p,n1,n2,&errV,&R1[0],&R2[0]);
  }while ((errV>*eps)&&(*niter<*maxiter));
}
