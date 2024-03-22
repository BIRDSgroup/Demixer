#include <omp.h>
#include<iostream>
#include <vector>
#include <iostream>
#include<ctime>
#include <stdio.h>
#include<gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include<math.h>
#include <map>
#include<cmath>
#include<chrono>
#include<random>
#include<string>
#include <openssl/rand.h>
#include <cstdint>
#include <array>
#include <cstdint>
#include <vector>
#include <unistd.h>
#include<cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include "globals.h"

int m=0;
int K=0;
long long int V=0;
char buff[256];
FILE *latfile;

int ii;
std::random_device rd;
char *newWorkingDirectory;
float *alpha;
float *beta;
long int epochtime[1000]={0};

short ***n_m_t;
short *Docs_1;

short *idf;
int known;
int *snptype;
int *snppos;
unsigned int **n_m_z, **n_z_t;
map<int, int**> variant_dict;
map<int, int**> *variant_map;
typedef pair<int, int**> variant_pair;

//initialization of variables
void set()
{
    int x1,x2;
    
    if (chdir(newWorkingDirectory) == 0) {
        // Directory change successful
        std::cout << "Working directory changed to: " << newWorkingDirectory << std::endl;
    }
    else
        std::cout << "Working directory not changed";   
    
    
    std::ifstream file("config.txt");   
    file>>m>>K>>V>>x1>>x2;
    file.close();    
    Docs_1=new short[m*V];
    snptype=new int[V];
    snppos=new int[V];
    idf=new short[V];
    alpha=new float[K];
    beta=new float[K];
    
    known=K-2;
    snprintf(buff,sizeof(buff),"%s","Docs_int_onlySnpgap.dat");
    latfile=fopen(buff,"rb");
    fread(Docs_1,sizeof(short),m*V,latfile);
    fclose(latfile);

	 for(int j=0;j<V;j++)
    	{
     if(j<x1)
        {
            snptype[j]=2;
            snppos[j]=0;}
        else if(j<x2)
        {
            snptype[j]=3;
            snppos[j]=x1;}
        else
        { 
            snptype[j]=4;
            snppos[j]=x2;}
    }
    
    sprintf(buff,"%s","idf.dat");

    latfile=fopen(buff,"rb");
    fread(idf,sizeof(short),V,latfile);
    fclose(latfile);
   
   
    snprintf(buff,sizeof(buff),"%s","alpha.dat");
    latfile=fopen(buff,"rb");
    fread(alpha,sizeof(float),K,latfile);
    fclose(latfile);
    
    snprintf(buff,sizeof(buff),"%s","beta.dat");
    latfile=fopen(buff,"rb");
    fread(beta,sizeof(float),K,latfile);
    fclose(latfile);
    
    /*
    for(int k=0;k<K;k++)
    {
        cout<<k<<"\talpha\t"<<alpha[k]<<"\t"<<"beta\t"<<beta[k]<<"\n";
        
    }
    */
    variant_map= new  map<int, int**> [m];
   
}

//initialization of variables for testing a new sample
void test_set()
{
    int x1,x2;

    cout<<newWorkingDirectory<<endl;
    
    if (chdir(newWorkingDirectory) == 0) {
        // Directory change successful
        std::cout << "Working directory changed to: " << newWorkingDirectory << std::endl;
    }
    else
        std::cout << "Working directory not changed";   
    std::ifstream file1("config.txt");   
    file1>>m>>K>>V>>x1>>x2;
    file1.close();
    cout<<m<<"\t"<<K;
    
    std::ifstream file2("../test/testconfig.txt"); 
    file2>>m;
    file2.close();
    
    Docs_1=new short[m*V];
    snptype=new int[V];
    snppos=new int[V];
    idf=new short[V];
    alpha=new float[K];
    beta=new float[K];
   
    known=K-2;
    snprintf(buff,sizeof(buff),"%s","../test/Docs_test.dat");
    latfile=fopen(buff,"rb");
    fread(Docs_1,sizeof(short),m*V,latfile);
    fclose(latfile);
	 for(int j=0;j<V;j++)
    	{
     if(j<x1)
        {
            snptype[j]=2;
            snppos[j]=0;}
        else if(j<x2)
        {
            snptype[j]=3;
            snppos[j]=x1;}
        else
        { 
            snptype[j]=4;
            snppos[j]=x2;}
    }
    
    sprintf(buff,"%s","idf.dat");
    latfile=fopen(buff,"rb");
    fread(idf,sizeof(short),V,latfile);
    fclose(latfile);
    
    for(int i=0;i<K;i++)
    {
        alpha[i]=1;
        beta[i]=0.01;
    }
    /*
    sprintf(buff,"%s","alpha.dat");
    latfile=fopen(buff,"rb");
    fread(alpha,sizeof(float),K,latfile);
    fclose(latfile);
    
   
    sprintf(buff,"%s","beta.dat");
    latfile=fopen(buff,"rb");
    fread(beta,sizeof(float),K,latfile);
    fclose(latfile);
    */
    
    unsigned int *n_z_t_temp= new unsigned int [K*V]();    
    
    n_m_t=new short **[m];
    n_m_z=new unsigned int *[m];
    n_z_t=new unsigned int *[K];
    
    for(int i=0;i<m;i++){
	n_m_t[i]=new short *[V];
	for(int j=0;j<V;j++)
		n_m_t[i][j]=new short [K]();
      }
    for(int i=0;i<m;i++)
        n_m_z[i]=new unsigned int [K]();

    for(int i=0;i<K;i++)
        n_z_t[i]=new unsigned int [V]();
    
    sprintf(buff,"%s","n_z_t0.dat");

    latfile=fopen(buff,"r");
    fread(n_z_t_temp,K * V * sizeof(unsigned int),K*V,latfile);
    fclose(latfile);
    
       for(int i=0;i<K;i++)
      {
          for(int j=0;j<V;j++)
          {
              n_z_t[i][j]=*(n_z_t_temp+i*V+j);
          }
      }  
}


//initialization of model parameters
void set_variables()
{
    n_m_t=new short **[m];
    n_m_z=new unsigned int *[m];
    n_z_t=new unsigned int *[K];
    
    for(int i=0;i<m;i++){
	n_m_t[i]=new short *[V];
	for(int j=0;j<V;j++)
		n_m_t[i][j]=new short [K]();
      }
   
    for(int i=0;i<m;i++)
        n_m_z[i]=new unsigned int [K]();

    for(int i=0;i<K;i++)
        n_z_t[i]=new unsigned int [V]();
    
}
//gsl_ran_multinomial sampling function modified to return short integers
void gsl_ran_multinomial_new (const gsl_rng * r, const size_t K, const unsigned int N, const double p[], short n[])
{
  size_t k;
  double norm = 0.0;
  double sum_p = 0.0;

  unsigned int sum_n = 0;

  for (k = 0; k < K; k++)
    {
      norm += p[k];
    }

  for (k = 0; k < K; k++)
    {
      if (p[k] > 0.0)
        {
	  unsigned int temp=gsl_ran_binomial (r, p[k] / (norm - sum_p), N - sum_n);
          n[k] = short(temp);
	
        }
      else
        {
          n[k] = 0;
        }

      sum_p += p[k];
      sum_n += n[k];
    }

}

//Initialization of model parameters
void initialize(std::uint64_t seed)
{
   //construction of sample-specific dictionary
   const char *fileName="temp2.txt";
   ifstream paramFile;
   paramFile.open(fileName);
   string line;
   int variant_id,vlength,length,sample_id;
   paramFile.seekg(0, std::ios::end);
   if (paramFile.tellg() == 0) {
        std::cout << "The file is empty.\n";
    } 
    else {
      paramFile.seekg(0, std::ios::beg);
       while ( paramFile.good() ){
            getline(paramFile, line);
            istringstream ss(line);
        ss>>sample_id>>length;
        for(int j=0;j<length;j++){
            getline(paramFile, line,' ');
            istringstream ss1(line);
            getline(paramFile, line,' ');
                    istringstream ss2(line);


            ss1>>variant_id;
            ss2>>vlength;

            int *variant_array=new int[vlength]();

                int *var1=new int();
            *var1=vlength;

            getline(paramFile, line,'\n');
                    istringstream ss3(line);

                for(int i=0;i<vlength;i++)

            {	
                ss3>>variant_array[i];

            }

            int **temp=new int*[2];
            temp[0]=var1;
            temp[1]=variant_array;	
            variant_map[sample_id].insert(variant_pair(variant_id,temp));
        }}
   }     
    
   std::cout<<"reading done"; 
   paramFile.close();
   double pmf[K];
   std::fill_n (pmf, K, double(1.0)/K);
    
   vector<int> z(K,0);
   time_t start = time(&start); 
   
   std::seed_seq seq{seed};
   std::vector<std::uint64_t> seeds(m);
 
   seq.generate(seeds.begin(), seeds.end());
   
   #pragma omp parallel for
   for(int i=0;i<m;i++)	
   { 
    int **val;
    double p,p1;
    unsigned int no_read;
    double *weightslist=new double[K]();
       
    gsl_rng * r;
    gsl_rng_env_setup();   

    mt19937_64 rng(seeds[i]);
   
    const gsl_rng_type * T;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, rng()); 
    std::fill_n(n_m_z[i], K, 0);
    for(int j=0;j<V;j++)
    {
        no_read=*(Docs_1+i*V+j);
        if(no_read>0){
	
        if(variant_map[i].count(j) > 0)
	{
		val=variant_map[i].at(j);

        p=double(0.00001);
        p1=0;
	    std::fill_n (weightslist, K, 0);
            
        if(*val[0]==1)
            weightslist[val[1][0]]=double(0.9999);
        else{
            weightslist[val[1][0]]=double(0.9500);
            //p1=double(0.9999)/((*val[0]));
            p1=double(0.0500)/((*val[0])-1);

        }  
		for(int k=1;k<(*val[0]);k++)
        {
            
			weightslist[val[1][k]]=p1;
        }
     
        if(known!=K)
        {
            weightslist[known]=p/2;
            weightslist[known+1]=p/2;
        }
        gsl_ran_multinomial_new(r, K, no_read, weightslist,n_m_t[i][j]);   
        
		for(int k=0;k<K;k++)
                {
				n_m_z[i][k]=n_m_z[i][k]+(static_cast<unsigned int>(idf[j]+1)*static_cast<unsigned int>(n_m_t[i][j][k]));
		}
          
		
	}

	else
	{
	gsl_ran_multinomial_new(r, K, no_read, pmf,n_m_t[i][j]);

	for(int k=0;k<K;k++)
                {              
                                n_m_z[i][k]=n_m_z[i][k]+(static_cast<unsigned int>(idf[j]+1)*static_cast<unsigned int>(n_m_t[i][j][k]));
                                       
                }
	}
	
	}
     }

    }
    

    #pragma omp parallel for
    for(int i=0;i<V;i++)
    {
        for(int j=0;j<m;j++)
        {
             for(int k=0;k<K;k++)
             {
                 n_z_t[k][i]=n_z_t[k][i]+static_cast<unsigned int>(n_m_t[j][i][k]);
             }
        }
        
    }
    
    	
    time_t end = time(&end);
    cout<<end-start<<" seconds"<<endl;
    //delete latfile;
    //delete variant_map;
       
}

void test_initialize(std::uint64_t e_seed)
{   
     
    time_t start = time(&start);

    std::seed_seq seq{e_seed};
    std::vector<std::uint64_t> seeds(m);
 
    seq.generate(seeds.begin(), seeds.end());
  
    #pragma omp parallel for
    for(int i=0;i<m;i++)
    {
  
    unsigned int no_read;
     
    double denom_b[K];
    double b[K][4];
    
    double *weightslist=new double[K]();
    double sum=0;
    int **dict;

    int length;
    double *weightslist1=new double[K]();
    double sum1=0;
       
    gsl_rng * r;
    gsl_rng_env_setup();   

    mt19937_64 rng(seeds[i]);
   
    const gsl_rng_type * T;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, rng());
    
  
    int st=0;
    for(int j=0;j<V;j++)
    {  
  
	 no_read=*(Docs_1+i*V+j);
   
     int type=snptype[j];
    
    
	if((j-snppos[j])%type==0)
	{
        
		for(int k=0;k<K;k++)
	        {
                
                denom_b[k]=0;
           
                for(int x=0;x<type;x++){
                    b[k][x]=n_z_t[k][j+x]*(static_cast<unsigned int>(idf[j+x]+1));

                    denom_b[k]=(denom_b[k]+b[k][x])+beta[k]; //beta
              
                }

           }
    st=0;
	}
	//cin>>temp;
	if(no_read>0)
	{
	        sum=0; 
      
                for(int k=0;k<K;k++)
                {       weightslist[k]=0;                    
                        weightslist[k]=((b[k][st]+beta[k])/(denom_b[k]));
                        sum=sum+weightslist[k];
                }

        
		if(variant_dict.count(j) <= 0 ){
		        
			if(sum>0)
	                {
        	        for(int k=0;k<K;k++)
         	       	{   
                        	weightslist[k]=weightslist[k]/sum;

                	}	
                
                	}
			gsl_ran_multinomial_new(r, K, no_read, weightslist,n_m_t[i][j]);

			for(int k=0;k<K;k++)
	                {
                     
		              	n_m_z[i][k]=n_m_z[i][k]+(static_cast<unsigned int>(idf[j]+1)* static_cast<unsigned int>(n_m_t[i][j][k]));
                        
        	        }
			
		}
		
	
		else 
		
		{
            dict=variant_dict.at(j); 
             
			length=*dict[0];
			sum1=0;
			fill_n (weightslist1, K, 0);
			for(int x=0;x<length;x++)
			{
			int k=dict[1][x];
			weightslist1[k]=weightslist[k];
			
			sum1=sum1+weightslist1[k];
			}	
			if(sum1>0)
			{
			for(int x=0;x<length;x++)
                        {
				int k=dict[1][x];
	            weightslist1[k]=weightslist1[k]/sum1;}
			}
			gsl_ran_multinomial_new(r, K, no_read, weightslist1,n_m_t[i][j]);
			for(int k=0;k<K;k++)
            {
                                n_m_z[i][k]=n_m_z[i][k]+(static_cast<unsigned int>(idf[j]+1)*static_cast<unsigned int>(n_m_t[i][j][k]));
                        }

					
		}
	}
	 st=st+1;
    }
   
   
    }
    
  
    for (int i = 0; i < K; i++) {
    for (int j = 0; j < V; j++) {
        n_z_t[i][j] = 0;
    }
}

    #pragma omp parallel for
    for(int i=0;i<V;i++)
    {
        
        for(int j=0;j<m;j++)
        {
             for(int k=0;k<K;k++)
             {
                 n_z_t[k][i]=n_z_t[k][i]+static_cast<unsigned int>(n_m_t[j][i][k]); 
             }

        }

    }

    time_t end = time(&end);
   
    epochtime[ii-1]=static_cast<long int>(end-start);

}


//Updation of n_m_z and n_z_t by reassigning strains to each read of samples

void epochs(std::uint64_t e_seed,int iter)
{   
    
  time_t start = time(&start);
  //cout<<"inside function\n";
  //cin>>temp;
    
    std::seed_seq seq{e_seed};
    std::vector<std::uint64_t> seeds(m);
 
    seq.generate(seeds.begin(), seeds.end());

   //cout<<iter;
   #pragma omp parallel for
   for(int i=0;i<m;i++)
   {
    
    unsigned int no_read;
    double a[K];  
    double denom_b[K];
    double b[K][4];
    
    double *weightslist=new double[K]();
    double sum=0;
    int **dict;

    int length;
    double *weightslist1=new double[K]();
    double sum1=0;

    gsl_rng * r;
   gsl_rng_env_setup();   

   mt19937_64 rng(seeds[i]);
   
   const gsl_rng_type * T;
   T = gsl_rng_default;
   r = gsl_rng_alloc (T);
   gsl_rng_set(r, rng());
    
    
    int st=0;
    for(int j=0;j<V;j++)
    {  
  
	 no_read=*(Docs_1+i*V+j);

     int type=snptype[j];
   
	if((j-snppos[j])%type==0)
	{
		for(int k=0;k<K;k++)
	        {

                denom_b[k]=0.0;
                a[k]=n_m_z[i][k];
                for(int x=0;x<type;x++){
                    b[k][x]=n_z_t[k][j+x]*(static_cast<unsigned int>(idf[j+x]+1));
                    n_m_z[i][k]=n_m_z[i][k]-(static_cast<unsigned int>(idf[j+x]+1)*static_cast<unsigned int>(n_m_t[i][j+x][k]));
                    //cout<<i<<"\t"<<beta[k]<<"\n";
                    denom_b[k]=(denom_b[k]+b[k][x])+beta[k];
                    n_m_t[i][j+x][k]=0;
                }

           }
    st=0;
	}

	//cin>>temp;
	if(no_read>0)
	{
	        sum=0; 

                for(int k=0;k<K;k++)
                {       weightslist[k]=0;   
                        weightslist[k]=(a[k]+alpha[k])*(((b[k][st])+beta[k])/(denom_b[k]));
                        sum=sum+weightslist[k];
                }
                              
        
		if(variant_dict.count(j) <= 0 ){
		        
			if(sum>0)
	                {
        	        for(int k=0;k<K;k++)
         	       	{   
                        	weightslist[k]=weightslist[k]/sum;

                	}	
                
                	}
			gsl_ran_multinomial_new(r, K, no_read, weightslist,n_m_t[i][j]);

			for(int k=0;k<K;k++)
	                {
                     
		              	n_m_z[i][k]=n_m_z[i][k]+(static_cast<unsigned int>(idf[j]+1)* static_cast<unsigned int>(n_m_t[i][j][k]));
                      
        	        }

			
		}
		
	
		else 
		
		{
            /*    //uncomment for simulation 2 quanttb - simDemixer       
            if(variant_map[i].count(j)>0)
            	dict=variant_map[i].at(j); 
			else
            */

            dict=variant_dict.at(j); 
             
			length=*dict[0];
			sum1=0;
			fill_n (weightslist1, K, 0);
			for(int x=0;x<length;x++)
			{
			int k=dict[1][x];
			weightslist1[k]=weightslist[k];
			
			sum1=sum1+weightslist1[k];
			}	
			if(sum1>0)
			{
			for(int x=0;x<length;x++)
                        {
				int k=dict[1][x];
	            weightslist1[k]=weightslist1[k]/sum1;}
			}
			gsl_ran_multinomial_new(r, K, no_read, weightslist1,n_m_t[i][j]);
			for(int k=0;k<K;k++)
            {
                                n_m_z[i][k]=n_m_z[i][k]+(static_cast<unsigned int>(idf[j]+1)*static_cast<unsigned int>(n_m_t[i][j][k]));
                        }

					
		}
	}
	 st=st+1;
    }
   
   
    }
    
  
    for (int i = 0; i < K; i++) {
    for (int j = 0; j < V; j++) {
        n_z_t[i][j] = 0;
    }
}
    //cin>>temp;
    
    #pragma omp parallel for
    for(int i=0;i<V;i++)
    {
        
        for(int j=0;j<m;j++)
        {
             for(int k=0;k<K;k++)
             {
                 n_z_t[k][i]=n_z_t[k][i]+static_cast<unsigned int>(n_m_t[j][i][k]); 
             }

        }

    }

    
    time_t end = time(&end);
   
    epochtime[ii-1]=static_cast<long int>(end-start);
    //cout<<ii<<"\t"<<loglik[ii-1]<<"\t"<<epochtime[ii-1]<<endl;
}
// construction of global dictionary
void construct_dict()
{

   typedef pair<int, int**> variant_pair;
   const char *fileName="temp1.txt";
   ifstream paramFile;
   paramFile.open(fileName);
   string line;
   int variant_id,vlength;
   paramFile.seekg(0, std::ios::end);
   if (paramFile.tellg() == 0) {
        std::cout << "The file is empty.\n";
    } 
    else {
       paramFile.seekg(0, std::ios::beg);
    while ( paramFile.good() ){
	getline(paramFile, line,' ');
	istringstream ss1(line);
	getline(paramFile, line,' ');
        istringstream ss2(line);	
	ss1>>variant_id;
	ss2>>vlength;
	
	int *variant_array=new int[vlength]();
		
        int *var1=new int();
	*var1=vlength;

	getline(paramFile, line,'\n');
        istringstream ss3(line);
        for(int i=0;i<vlength;i++)	
		ss3>>variant_array[i];
	int **temp=new int*[2];
	temp[0]=var1;
	temp[1]=variant_array;	
	variant_dict.insert(variant_pair(variant_id,temp));
	}

    }
    paramFile.close();

}
