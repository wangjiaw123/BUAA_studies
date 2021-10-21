/**************************************************
 * Author : 王佳文
 * Date : 2020/12/22
 * Function Name : sfls_type2_C
 * Description : 计算单值二型模糊逻辑系统输出的C语言程序
 *         
**************************************************/

#include "stdio.h"
#include "mex.h"
#include "matrix.h"
#include "math.h"

void gausstype2(double *u,double *l,double Xt_gs,double sigma_gs,double M1_gs,double M2_gs);
double bubble(double * a, int n,int * idx);
void rightpoint(double *r_out_r,int *I2l_r,int *I2u_r,double *wr_r,double *c2_r,double LL_r[],double UU_r[],int N);
void leftpoint(double *l_out_l,int *I1l_l,int *I1u_l,double *wl_l,double *c1_l,double LL_l[],double UU_l[],int N);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])  
{
    if(nrhs != 6)
    {
        mexErrMsgIdAndTxt("MATLAB:Input Error!","There must be six input,which are X,M1,M2,c1,c2,sigma.");
    }
    if(nlhs != 3)
    {
        mexErrMsgIdAndTxt("MATLAB:Output Error!","There must be three output,which are R1_out,R2_out,R_out.");
    }   
    int lr_sum=0;
    for(int lr=0;lr<6;lr++)         
    {
        lr_sum+=mxIsDouble(prhs[lr]);
    }    
    if(lr_sum != 6)
    {
        mexErrMsgIdAndTxt("MATLAB:Input Error!","The six input must be stores its data as double-precision,floating-point number.");
    }

    double * X, * M1, * M2 , * c1, * c2, * sigma;
    double * R1_out, * R2_out, * R_out;
    int L,N,n;
    
    L=mxGetM(prhs[0]);        
    N=mxGetM(prhs[2]);        
    n=mxGetN(prhs[2]);       

    X=mxGetPr(prhs[0]);     
    M1=mxGetPr(prhs[1]);
    M2=mxGetPr(prhs[2]);
    c1=mxGetPr(prhs[3]);    
    c2=mxGetPr(prhs[4]);    
    sigma=mxGetPr(prhs[5]);    

    plhs[0]=mxCreateDoubleMatrix(L,1,mxREAL);  
    plhs[1]=mxCreateDoubleMatrix(L,1,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(L,1,mxREAL);

    R1_out=mxGetPr(plhs[0]);     
    R2_out=mxGetPr(plhs[1]);
    R_out=mxGetPr(plhs[2]);

    double UU[N],LL[N];
    double Uu=1.0,Ll=1.0;
    double uu,ll;
    uu=ll=1;
    double *Uuu,*Lll;
    Uuu=&uu;
    Lll=&ll;
    
    for(int i=0;i<L;i++)
    {
        //printf("i=%d \n",i);
        for(int j=0;j<N;j++)
        {
            Uu=1.0;
            Ll=1.0;
            for(int m=0;m<n;m++)
            {
                gausstype2(Uuu,Lll,*(X+m*L+i),*(sigma+m*N+j),*(M1+m*N+j),*(M2+m*N+j));
                Uu=Uu* *Uuu;
                Ll=Ll* *Lll;
                //printf("*Uuu=%f,*Lll=%f \n",*Uuu,*Lll);
            }
            UU[j]=Uu;
            LL[j]=Ll;
        }
        
       
        double * r_out_p,* wr_p;  
        double r_out,wr[N];   
        int I2l_r[N],I2u_r[N],* I2l_r_p,* I2u_r_p;
        r_out_p=&r_out;
        I2l_r_p=I2l_r;     
        I2u_r_p=I2u_r;
        wr_p=wr;
      
        for(int ij=0;ij<N;ij++)     
        {
            I2l_r[ij]=0;
            I2u_r[ij]=0;
            //printf("LL[i]=%f,UU[i]=%f,*(c2+i)=%f \n",*(c2+i),LL[i],UU[i]);
        }                
        rightpoint(r_out_p,I2l_r_p,I2u_r_p,wr_p,c2,LL,UU,N);  
     
        double * l_out_p,* wl_p;  
        double l_out,wl[N];   
        int I1l_l[N],I1u_l[N],* I1l_l_p,* I1u_l_p;
        l_out_p=&l_out;
        I1l_l_p=I1l_l;
        I1u_l_p=I1u_l;
        wl_p=wl;     
        
        for(int ij=0;ij<N;ij++)
        {
            I1l_l[ij]=0;
            I1u_l[ij]=0;
        }          
                
        leftpoint(l_out_p,I1l_l_p,I1u_l_p,wl_p,c1,LL,UU,N);
        
       *(R1_out+i)=*l_out_p;
       *(R2_out+i)=*r_out_p;
       *(R_out+i)=(*l_out_p+*r_out_p)/2.0;
          
    } // end for(int i=0;i<L;i++)
     
}   //end mexFunction



//////////////////////////////////

void gausstype2(double *u,double *l,double Xt_gs,double sigma_gs,double M1_gs,double M2_gs)      
{ 
    double m1,m2,mu1,mu2;
    sigma_gs=fabs(sigma_gs)+0.000001;
    m1=(M1_gs<M2_gs?M1_gs:M2_gs);
    m2=(M1_gs>M2_gs?M1_gs:M2_gs);
    
    if((Xt_gs>=m1)&&(Xt_gs<=(m1+m2)/2.0))
    {
        mu1=1.0;
        mu2=exp(-pow(Xt_gs-m2,2.0)/(2.0*pow(sigma_gs,2.0)));
    }else if((Xt_gs>((m1+m2)/2.0))&&(Xt_gs<=m2))
    {
        mu1=1.0;
        mu2=exp(-pow(Xt_gs-m1,2.0)/(2.0*pow(sigma_gs,2.0)));        
    }else
    {
        mu1=exp(-pow(Xt_gs-m1,2.0)/(2.0*pow(sigma_gs,2.0)));
        mu2=exp(-pow(Xt_gs-m2,2.0)/(2.0*pow(sigma_gs,2.0)));  
    }
    *u=(mu1>mu2?mu1:mu2);
    *l=(mu1<mu2?mu1:mu2); 
}



////////////////////////////////////////

double bubble(double * a, int n,int * idx) 
{
    int IDx[n];
    for(int k=0;k<n;k++)
    {
        IDx[k]=k+1;
    }
    
	for(int i=0; i<n-1; i++) 
    {
		for(int j=0; j<n-1-i; j++)
        {
			if(*(a+j) > *(a+j+1)) 
            {
				double t = *(a+j);
				*(a+j) = *(a+j+1);
			    *(a+j+1)=t;
                int tt=IDx[j];
                IDx[j]=IDx[j+1];
                IDx[j+1]=tt;
            }
        }
    }
    
    for(int k=0;k<n;k++)
    {
       *(idx+k)= IDx[k];
    }
}

///////////////////////////////

void rightpoint(double *r_out_r,int *I2l_r,int *I2u_r,double *wr_r,double *c2_r,double LL_r[],double UU_r[],int N)
{
    double b2[N],L[N],U[N],* b2_p;
    int I2[N],* I2_p;
    I2_p=I2;
    b2_p=b2;
    
    for(int j=0;j<N;j++)
    {
        *(b2_p+j)=*(c2_r+j);
    }
    bubble(b2_p,N,I2_p);
    
  
    for(int i=0;i<N;i++)
    {
        L[i]=LL_r[(int)(*(I2+i)-1)];
        U[i]=UU_r[(int)(*(I2+i)-1)];
        //printf("(*(I2+i)-1)=%f,(*(I2+i)-1)=%f",(int)(*(I2+i)-1),(int)(*(I2+i)-1));
        //printf("L[i]=%f,U[i]=%f \n",L[i],U[i]);
    }
    
    double r_out,s,s1;
    r_out=s=s1=0;
    for(int i1=0;i1<N;i1++)
    {
        s1+=U[i1];
        s+=b2[i1]*U[i1];
    }
    r_out=s/s1;
    
    int I2l_help1[N],I2u_help1[N];
    for(int d1=0;d1<N;d1++)
    {
        I2l_help1[d1]=0;
        I2u_help1[d1]=0;
    }
    
    for(int i2=0;i2<N;i2++)
    {
        s=s-b2[i2]*U[i2]+b2[i2]*L[i2];
        s1=s1-U[i2]+L[i2];
        //printf("s=%f,s1=%f,s/s1=%f \n",s,s1,s/s1);
        if(s/s1>r_out)
        {
            r_out=s/s1;
            I2l_help1[i2]=I2[i2];
        }else
        {
            I2u_help1[i2]=I2[i2];
        } 
    }
    
    *r_out_r=r_out;
    
    int count1=0;
    for(int i3=0;i3<N;i3++)
    {
        if((int)I2u_help1[i3]!=0){count1++;}
    }
    int I2u_help2[count1];

    int j1=0;
    for(int i4=0;i4<N;i4++)
    {
        if(I2u_help1[i4]!=0)
        {
            I2u_help2[j1]=I2u_help1[i4];
            j1++;
        }
    }
    
    int count2=0;
    for(int i5=0;i5<N;i5++)
    {
        if((int)I2l_help1[i5]!=0){count2++;}
    }
    int I2l_help2[count2];

    int j2=0;
    for(int i6=0;i6<N;i6++)
    {
        if(I2l_help1[i6]!=0)
        {
            I2l_help2[j2]=I2l_help1[i6];
            j2++;
        }
    }   
    
    
    for(int i7=0;i7<count2;i7++)
    {
        *(wr_r+(I2l_help2[i7]-1))=LL_r[(I2l_help2[i7]-1)]; 
        *(I2l_r+i7)=I2l_help2[i7];
        //printf("right i7=%d,*(I2l_r+i7)=%d \n",i7,*(I2l_r+i7));        
    }
    for(int i8=0;i8<count1;i8++)
    {
        *(wr_r+(I2u_help2[i8]-1))=UU_r[(I2u_help2[i8]-1)];
        *(I2u_r+i8)=I2u_help2[i8];
    } 

}

////////////////////////////////////

void leftpoint(double *l_out_l,int *I1l_l,int *I1u_l,double *wl_l,double *c1_l,double LL_l[],double UU_l[],int N)
{
    double b2[N],L[N],U[N],* b2_p;
    int I2[N],* I2_p;
    I2_p=I2;
    b2_p=b2;
    
    for(int j=0;j<N;j++)
    {
        *(b2_p+j)=*(c1_l+j);
    }
    bubble(b2_p,N,I2_p);
    
    for(int i=0;i<N;i++)
    {
        L[i]=LL_l[(int)(*(I2+i)-1)];
        U[i]=UU_l[(int)(*(I2+i)-1)];
    }
    
    double l_out,s,s1;
    l_out=s=s1=0;
    for(int i1=0;i1<N;i1++)
    {
        s1+=L[i1];
        s+=b2[i1]*L[i1];
    }
    l_out=s/s1;
    
    int I1l_help1[N],I1u_help1[N];
    for(int d1=0;d1<N;d1++)
    {
        I1l_help1[d1]=0;
        I1u_help1[d1]=0;
    }
    
    for(int i2=0;i2<N;i2++)
    {
        s=s-b2[i2]*L[i2]+b2[i2]*U[i2];
        s1=s1+U[i2]-L[i2];
        if(s/s1<l_out)
        {
            l_out=s/s1;
            I1u_help1[i2]=I2[i2];
        }else
        {
            I1l_help1[i2]=I2[i2];
        } 
    }
    
    *l_out_l=l_out;
    
    int count1=0;
    for(int i3=0;i3<N;i3++)
    {
        if((int)I1u_help1[i3]!=0){count1++;}
    }
    int I1u_help2[count1];

    int j1=0;
    for(int i4=0;i4<N;i4++)
    {
        if(I1u_help1[i4]!=0)
        {
            I1u_help2[j1]=I1u_help1[i4];
            j1++;
        }
    }
    
    int count2=0;
    for(int i5=0;i5<N;i5++)
    {
        if((int)I1l_help1[i5]!=0){count2++;}
    }
    int I1l_help2[count2];

    int j2=0;
    for(int i6=0;i6<N;i6++)
    {
        if(I1l_help1[i6]!=0)
        {
            I1l_help2[j2]=I1l_help1[i6];
            j2++;
        }
    }    
    
    
    for(int i7=0;i7<count1;i7++)
    {
        *(wl_l+(I1u_help2[i7]-1))=UU_l[(I1u_help2[i7]-1)];  
        *(I1u_l+i7)=I1u_help2[i7];      
    }
    for(int i8=0;i8<count2;i8++)
    {
        *(wl_l+(I1l_help2[i8]-1))=LL_l[(I1l_help2[i8]-1)];
        *(I1l_l+i8)=I1l_help2[i8];
    } 
    
}



