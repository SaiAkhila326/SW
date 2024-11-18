#include <stdio.h>
#include <math.h>

#define MAX_ITER 10000  
#define TOL 1e-6

void copy(int n,double A[n][n], double B[n][n]) {
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            B[i][j]=A[i][j];
        }
    }
}

void multiply(int n, double A[n][n], double B[n][n], double C[n][n]) {
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
           C[i][j] = 0; 
            for (int k=0;k<n;k++) {
                C[i][j]+= A[i][k]*B[k][j]; 
            }
        }
    }
}
void transpose(int n, double A[n][n], double B[n][n]) {
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            B[j][i]= A[i][j];
        }
    }
}
void gs(int n,double a[n][n],double q[n][n],double r[n][n]){
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            r[i][j]=0;
        }
    }
    for(int j=0;j<n;j++)
    {
        double norm=0;
        for(int i=0;i<n;i++)
        {
            q[i][j]=a[i][j];
        }
        for(int k=0;k<j;k++)
        {
            for(int i=0;i<n;i++)
            {
                r[k][j]+=q[i][j]*q[i][k];
            }
            for(int i=0;i<n;i++)
            {
                q[i][j]-=r[k][j]*q[i][k];
            }
        }
        for(int i=0;i<n;i++)
        {
         norm+= q[i][j]*q[i][j];
        }
        norm=sqrt(norm);
        r[j][j]=norm;
        if(norm>1e-8)
        for(int i=0;i<n;i++)
        {
            q[i][j]/=norm;
        }
    }
}
/*void qr(int n, double A[n][n], double eigenvalues[n]) {
    double Q[n][n], R[n][n], A_next[n][n];
    copy_matrix(n, A, A_next);

    for (int iter = 0; iter < MAX_ITER; iter++) {
        // Perform QR decomposition: A_k = Q * R
        gs(n, A_next, Q, R);

        // Compute A_k+1 = R * Q
        matrix_multiply(n, R, Q, A_next);

        // Check for convergence (off-diagonal elements close to zero)
        int converged = 1;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && fabs(A_next[i][j]) > TOL) {
                    converged = 0;
                    break;
                }
            }
            if (!converged) break;
        }
        if (converged) break;
    }}*/
void hessenberg_reduction(int n, double A[n][n]) {
    for (int k = 0; k < n - 2; k++) {
        double norm = 0.0;
        for (int i = k + 1; i < n; i++) {
            norm += A[i][k] * A[i][k];
        }
        norm = sqrt(norm);
        if (fabs(norm) < TOL) 
        continue;

        double u[n];
        for (int i = 0;i < n;i++) {
            u[i]=0.0;
        }
        u[k+1] = A[k+1][k] + (A[k+1][k] > 0?norm:-norm);
        for (int i = k + 2; i < n; i++) {
            u[i]= A[i][k];
        }
        double u_norm = 0.0;
        for (int i = 0; i < n; i++) {
            u_norm+= u[i] * u[i];
        }
        u_norm =sqrt(u_norm);

        for (int i =0; i <n;i++) {
            u[i]/= u_norm;
        }

        double H[n][n];
        for (int i= 0; i < n; i++) {
            for (int j = 0;j <n;j++) {
                H[i][j] = -2.0 * u[i] * u[j];
            }
            H[i][i] += 1.0;
        }

        double temp[n][n];
        multiply(n, H, A, temp);
        multiply(n, temp, H, A);
    }
}
int main(void)
{
    int n;
    scanf("%d",&n);
    double a[n][n],q[n][n],r[n][n];double an[n][n];
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            scanf("%lf",&a[i][j]);
            //r[i][j]=0;
        }
    }
    hessenberg_reduction(n, a);
    //gs(n,a,q,r);
    double eigenvalues[n];copy(n,a,an);
    for(int i=0;i<MAX_ITER;i++)
    {
        //copy(n,a,an);
        gs(n,an,q,r);
        multiply(n,r,q,an);int count=0;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                if(i!=j&&fabs(an[i][j])<=TOL)
                count++;
            }
        }
        if(count==n*n-n)
        break;
        
    }
   /* for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            printf("%lf ",an[i][j]);
        }
        printf("\n");
    }*/
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i==j)
            eigenvalues[i]=an[i][i];
        }
    }
    

    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%.6f\n", eigenvalues[i]);
    }
   /*for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            printf("%lf ",q[i][j]);
        }
        printf("\n");
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            printf("%lf ",r[i][j]);
        }
        printf("\n");
    }*/
}
