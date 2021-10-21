/**************************************************
 * Author : 王佳文
 * Date : 2020/12/21
 * Function Name : Train_sfls_type2_C
 * Description : 单值二型模糊逻辑系统参数训练的C语言程序
 *         
**************************************************/

#include "stdio.h"
#include "mex.h"
#include "matrix.h"
#include "math.h"

double SUM(double A[],int len_A);
void gausstype2(double *u,double *l,double Xt_gs,double sigma_gs,double M1_gs,double M2_gs);
double bubble(double * a, int n,int * idx);
void rightpoint(double *r_out_r,int *I2l_r,int *I2u_r,double *wr_r,double *c2_r,double LL_r[],double UU_r[],int N);
void leftpoint(double *l_out_l,int *I1l_l,int *I1u_l,double *wl_l,double *c1_l,double LL_l[],double UU_l[],int N);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])   
{
    if(nrhs != 8)
    {
        mexErrMsgIdAndTxt("MATLAB:Input Error!","There must be six input,which are X_train,Y_train,M1,M2,sigma,c1,c2,alpha.");
    }
    if(nlhs != 5)
    {
        mexErrMsgIdAndTxt("MATLAB:Output Error!","There must be five output,which are M1_out,M2_out,c1_out,c2_out,sigma_out.");
    }   
    int lr_sum=0;
    for(int lr=0;lr<8;lr++)          
    {
        lr_sum+=mxIsDouble(prhs[lr]);
    }    
    if(lr_sum != 8)
    {
        mexErrMsgIdAndTxt("MATLAB:Input Error!","The six input must be stores its data as double-precision,floating-point number.");
    }
    
    double *Xt,*Dt,*M1,*M2,*sigma,*c1,*c2,alpha;         
    double *M1_out,*M2_out,*c1_out,*c2_out,*sigma_out; 
    int L,N,n;
    
    L=mxGetM(prhs[0]);        
    N=mxGetM(prhs[2]);        
    n=mxGetN(prhs[2]);        
    
    Xt=mxGetPr(prhs[0]);    
    Dt=mxGetPr(prhs[1]);
    M1=mxGetPr(prhs[2]);
    M2=mxGetPr(prhs[3]);
    c1=mxGetPr(prhs[4]);    
    c2=mxGetPr(prhs[5]);
    sigma=mxGetPr(prhs[6]);
    alpha=mxGetScalar(prhs[7]); 
    
    //printf("*********************************************");
    //printf("M1=%f,M2=%f,c1=%f,c2=%f,sigma=%f \n",*(M1+N+1),*(M2+N+1),*(c1+1),*(c2+1),*(sigma+N+1));
    
    plhs[0]=mxCreateDoubleMatrix(N,n,mxREAL);  
    plhs[1]=mxCreateDoubleMatrix(N,n,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(N,1,mxREAL);
    plhs[3]=mxCreateDoubleMatrix(N,1,mxREAL);
    plhs[4]=mxCreateDoubleMatrix(N,n,mxREAL);
    
    M1_out=mxGetPr(plhs[0]);     
    M2_out=mxGetPr(plhs[1]);
    c1_out=mxGetPr(plhs[2]);
    c2_out=mxGetPr(plhs[3]);
    sigma_out=mxGetPr(plhs[4]);
   
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
                //printf("*(Xt+m*L+i)=%f,*(sigma+m*N+j)=%f,*(M1+m*N+j)=%f,*(M2+m*N+j)=%f \n",*(Xt+m*L+i),*(sigma+m*N+j),*(M1+m*N+j),*(M2+m*N+j));
                gausstype2(Uuu,Lll,*(Xt+m*L+i),*(sigma+m*N+j),*(M1+m*N+j),*(M2+m*N+j));
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
        
        
        //printf("l_out=%f,r_out=%f,l_out+r_out=%f \n",l_out,l_out+r_out); 
        
       
        
  /*
        for(int i=0;i<N;i++)
        {
            printf("%d*",I1l_l[i]);
        }   
        printf("***********  \n");  
       for(int i=0;i<N;i++)
        {
            printf("%d*",I1u_l[i]);
        }           
    */      
        
        int len_I2u=0;
        for(int s1=0;s1<N;s1++)
        {
            if(*(I2u_r+s1)!=0){len_I2u++;}
        }
        int I2u[len_I2u];
        int j1=0;
        for(int s2=0;s2<N;s2++)
        {
            if(I2u_r[s2]!=0)
            {
                I2u[j1]=I2u_r[s2]-1;
                j1++;
            }
        }   
       
        int len_I2l=0;
        for(int s3=0;s3<N;s3++)
        {
            if(*(I2l_r+s3)!=0){len_I2l++;}
        }
        int I2l[len_I2l];
        int j2=0;
        for(int s4=0;s4<N;s4++)
        {
            if(I2l_r[s4]!=0)
            {
                I2l[j2]=I2l_r[s4]-1;
                j2++;
            }
        } 
        //printf("j1+j2=%d \n",j1+j2);
      
        int len_I1u=0;
        for(int s5=0;s5<N;s5++)
        {
            if(*(I1u_l+s5)!=0){len_I1u++;}
        }
        int I1u[len_I1u];
        int j3=0;
        for(int s6=0;s6<N;s6++)
        {
            if(I1u_l[s6]!=0)
            {
                I1u[j3]=I1u_l[s6]-1;
                j3++;
            }
        }           
        
        int len_I1l=0;
        for(int s7=0;s7<N;s7++)
        {
            if(*(I1l_l+s7)!=0){len_I1l++;}
        }
        int I1l[len_I1l];
        int j4=0;
        for(int s8=0;s8<N;s8++)
        {
            if(I1l_l[s8]!=0)
            {
                I1l[j4]=I1l_l[s8]-1;
                j4++;
            }
        } 
        //printf("j3+j4=%d \n",j3+j4);
        
       

        //printf("len_I2u+len_I2l=%d,len_I1u+len_I1l=%d \n",len_I2u+len_I2l,len_I1u+len_I1l);
   
        double f,e,fa1[N],fa2[N],sum_wl,sum_wr;
        f=(l_out+r_out)/2.0;
        e=*(Dt+i)-f;
        sum_wl=SUM(wl,N);
        sum_wr=SUM(wr,N);
        for(int fa12=0;fa12<N;fa12++)
        {
            fa1[fa12]=wl[fa12]/sum_wl;
            fa2[fa12]=wr[fa12]/sum_wr;
        }
        //printf("############################## \n");
 
        for(int t=0;t<len_I2u;t++)
        {
            for(int k=0;k<n;k++)
            {
                int len=I2u[t];
                if(*(Xt+k*L+i)< *(M1+k*N+len))
                {
                    *(M1+k*N+len)=*(M1+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M1+k*N+len))     
                                  /pow(*(sigma+k*N+len),2.0)*(*(c2+len)-r_out)*(wr[len])/sum_wr;
                    *(sigma+k*N+len)=*(sigma+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M1+k*N+len))
                                  /pow(*(sigma+k*N+len),3.0)*(*(c2+len)-r_out)*(wr[len])/sum_wr;
                    
                }else if(*(Xt+k*L+i)> *(M2+k*N+len))
                {
                    *(M2+k*N+len)=*(M2+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len)) 
                                  /pow(*(sigma+k*N+len),2.0)*(*(c2+len)-r_out)*(wr[len])/sum_wr;
                    *(sigma+k*N+len)=*(sigma+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len))
                                  /pow(*(sigma+k*N+len),3.0)*(*(c2+len)-r_out)*(wr[len])/sum_wr;                    
                }
            }
        }
        
        for(int t=0;t<len_I2l;t++)
        {
            for(int k=0;k<n;k++)
            {
                int len=I2l[t];
                if(*(Xt+k*L+i)< (*(M2+k*N+len)+*(M1+k*N+len))/2.0)
                {
                    *(M2+k*N+len)=*(M2+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len))     
                                  /pow(*(sigma+k*N+len),2.0)*(*(c2+len)-r_out)*(wr[len])/sum_wr;
                    *(sigma+k*N+len)=*(sigma+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len))
                                  /pow(*(sigma+k*N+len),3.0)*(*(c2+len)-r_out)*(wr[len])/sum_wr;
                    
                }else 
                {
                    *(M2+k*N+len)=*(M2+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len)) 
                                  /pow(*(sigma+k*N+len),2.0)*(*(c2+len)-r_out)*(wr[len])/sum_wr;
                    *(sigma+k*N+len)=*(sigma+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len))
                                  /pow(*(sigma+k*N+len),3.0)*(*(c2+len)-r_out)*(wr[len])/sum_wr;                    
                }
            }
        }      
        
         for(int t=0;t<len_I1l;t++)
        {
            for(int k=0;k<n;k++)
            {
                int len=I1l[t];
                if(*(Xt+k*L+i)< (*(M2+k*N+len)+*(M1+k*N+len))/2.0)
                {
                    *(M2+k*N+len)=*(M2+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len))     
                                  /pow(*(sigma+k*N+len),2.0)*(*(c1+len)-l_out)*(wl[len])/sum_wl;
                    *(sigma+k*N+len)=*(sigma+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len))
                                  /pow(*(sigma+k*N+len),3.0)*(*(c1+len)-l_out)*(wl[len])/sum_wl;
                    
                }else 
                {
                    *(M2+k*N+len)=*(M2+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M1+k*N+len)) 
                                  /pow(*(sigma+k*N+len),2.0)*(*(c1+len)-l_out)*(wl[len])/sum_wl;
                    *(sigma+k*N+len)=*(sigma+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M1+k*N+len))
                                  /pow(*(sigma+k*N+len),3.0)*(*(c1+len)-l_out)*(wl[len])/sum_wl;                    
                }
            }
         }
             
        for(int t=0;t<len_I1u;t++)
        {
            for(int k=0;k<n;k++)
            {
                int len=I1u[t];
                if(*(Xt+k*L+i)< *(M1+k*N+len))
                {
                    *(M1+k*N+len)=*(M1+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M1+k*N+len))     
                                  /pow(*(sigma+k*N+len),2.0)*(*(c1+len)-l_out)*(wl[len])/sum_wl;
                    *(sigma+k*N+len)=*(sigma+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M1+k*N+len))
                                  /pow(*(sigma+k*N+len),3.0)*(*(c1+len)-l_out)*(wl[len])/sum_wl;
                    
                }else if(*(Xt+k*L+i)> *(M2+k*N+len))
                {
                    *(M2+k*N+len)=*(M2+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len)) 
                                  /pow(*(sigma+k*N+len),2.0)*(*(c1+len)-l_out)*(wl[len])/sum_wl;
                    *(sigma+k*N+len)=*(sigma+k*N+len)+alpha*0.5*e*(*(Xt+k*L+i)-*(M2+k*N+len))
                                  /pow(*(sigma+k*N+len),3.0)*(*(c1+len)-l_out)*(wl[len])/sum_wl;                    
                }
            }
        } 
        

      for(int k1=0;k1<N;k1++)
      {
          for(int k2=0;k2<n;k2++)
          {
              if(*(sigma+k2*N+k1) < 0)
              {
                  *(sigma+k2*N+k1)=fabs(*(sigma+k2*N+k1));
              }
          }
      }
    
    for(int k3=0;k3<N;k3++)
    {
        *(c1+k3)=*(c1+k3)+alpha*e*fa1[k3]/2.0;
        *(c2+k3)=*(c2+k3)+alpha*e*fa2[k3]/2.0;
    }
  
   
    
    } 
 
 
    // 使用mxDestroyArray()函数
    
    /*
    plhs[0]=mxGetPr(M1);
    plhs[1]=mxGetPr(M2);
    plhs[2]=mxGetPr(c1);
    plhs[3]=mxGetPr(c2);
    plhs[4]=mxGetPr(sigma); 
   
    plhs[0]=prhs[2];
    plhs[1]=mxGetPr(prhs[3]);
    plhs[2]=mxGetPr(prhs[4]);
    plhs[3]=mxGetPr(prhs[5]);
    plhs[4]=mxGetPr(prhs[6]); 
     */
    
    // 返回值
    
    for(int k4=0;k4<N;k4++)
    {
        for(int k5=0;k5<n;k5++)
        {
            *(M1_out+k5*N+k4)=*(M1+k5*N+k4);
            *(M2_out+k5*N+k4)=*(M2+k5*N+k4);
            *(sigma_out+k5*N+k4)=*(sigma+k5*N+k4);
        }
        *(c1_out+k4)=*(c1+k4);
        *(c2_out+k4)=*(c2+k4);        
    }
    /*
    for(int k4=0;k4<N;k4++)
    {
        for(int k5=0;k5<n;k5++)
        {
            *(M1_out+k5*N+k4)=M1_h[k5*N+k4];
            *(M2_out+k5*N+k4)=M2_h[k5*N+k4];
            *(sigma_out+k5*N+k4)=sigma_h[k5*N+k4];
        }
        *(c1_out+k4)=c1_h[k4];
        *(c2_out+k4)=c2_h[k4];        
    }
    */

}  //end mexFunction



/////////////////////////////////////

double SUM(double A[],int len_A)
{
    double a=0;
    for(int i=0;i<len_A;i++)
    {
        a+=A[i];
    }
    return a;
}

//////////////////////////////////

void gausstype2(double *u,double *l,double Xt_gs,double sigma_gs,double M1_gs,double M2_gs)      
{
    //printf("gausstype2 done! \n");
    
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
    //printf("bubble done! \n");
    
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
    //printf("compute rightpoint done! \n");
    
    double b2[N],L[N],U[N],* b2_p;
    int I2[N],* I2_p;
    I2_p=I2;
    b2_p=b2;
    
    for(int j=0;j<N;j++)
    {
        *(b2_p+j)=*(c2_r+j);
        //printf("*(b2_p+j)=%f \n",*(b2_p+j));
    }
    bubble(b2_p,N,I2_p);
    
   /* 
        for(int i=0;i<N;i++)
        {
            printf("%d*",I2[i]);
        }   
    */
  
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
        //printf("right i8=%d,*(I2u_r+i8)=%d \n",*(I2u_r+i8));
    } 
    
    //printf("leftpoint count1=%d,count2=%d count1+count2=%d \n",count1,count2,count1+count2);
}

////////////////////////////////////

void leftpoint(double *l_out_l,int *I1l_l,int *I1u_l,double *wl_l,double *c1_l,double LL_l[],double UU_l[],int N)
{
    //printf("compute leftpoint done! \n");
   
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
        //printf("left i7=%d,*(I1l_l+i7)=%d \n",i7,*(I1l_l+i7));        
    }
    for(int i8=0;i8<count2;i8++)
    {
        *(wl_l+(I1l_help2[i8]-1))=LL_l[(I1l_help2[i8]-1)];
        *(I1l_l+i8)=I1l_help2[i8];
        //printf("left i8=%d,*(I1u_l+i8)=%d \n",i8,*(I1u_l+i8));        
    } 
    //printf("leftpoint count1=%d,count2=%d count1+count2=%d \n",count1,count2,count1+count2);
}

////////////////////////////////////////








