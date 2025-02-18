#include "globals.h"
#include <omp.h>

int main(int argc,char* argv[])
{

    newWorkingDirectory=argv[1];
    cout<<newWorkingDirectory; 
    
    omp_set_num_threads(std::atoi(argv[2]));
    set();
    construct_dict();

    fstream myfile;
    int no_iter=1000;
    unsigned int a,b,c;
    std::ofstream seedFile("output.txt"); //writing seeds to output file
    for(int r_iter=0;r_iter<1;r_iter++)
    {
    cout<<"iter_num:   "<<r_iter<<endl;
    a=rd();
    b=rd();
    c=rd();

    //a=2178893639;
    //b=1228918992;
    //c=2548831504;
    a=1613498695;
    b=432825264;
    c=460035106;
    seedFile<<a<<"\t"<<b<<"\t"<<c<<endl;
    seed_seq seq{a,b,c};
    std::vector<std::uint64_t> seeds(1001);
    seq.generate(seeds.begin(), seeds.end());
    set_variables();     

    initialize(seeds[0]);

    unsigned int *n_m_z_temp= new unsigned int [m*K]();  
    unsigned int *n_z_t_temp= new unsigned int [K*V]();    
    
    
    time_t start = time(&start);
    for(ii=1;ii<no_iter;ii++){

        cout<<ii<<"\n";
    	epochs(seeds[ii],ii);
        
       if(ii%100==0)
        {
		myfile.open("n_m_z"+std::to_string(r_iter)+".txt",fstream::out);

    		for(int x=0;x<m;x++)
            {   double sum=0;
        		for(int y=0;y<K;y++)
        		{
                	sum=sum+n_m_z[x][y];
        		}
        		for(int y=0;y<K;y++)
        		{
                	myfile<<round(n_m_z[x][y]/sum * 100) / 100<<"\t";
        		}
        	myfile<<endl;


            }
		myfile.close();

		unsigned int *n_m_z_temp= new unsigned int [m*K]();
		for (int i = 0; i < m; i++) {
		for (int j = 0; j < K; j++) {
                    *(n_m_z_temp+i*K+j)=n_m_z[i][j];
                        }}
                        std::ofstream outFile1("n_m_z"+std::to_string(r_iter)+".dat", std::ios::binary);
                        outFile1.write(reinterpret_cast<const char*>(n_m_z_temp), m*K*sizeof(unsigned int));
                          outFile1.close();

        } 
    }
    time_t end = time(&end);
    cout<<"Total time "<<static_cast<long int>(end-start);


    //writing output to file   
    myfile.open("n_m_z"+std::to_string(r_iter)+".txt",fstream::out);

        
        
    //#pragma omp parallel for
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<K;j++)
        {
             n_m_z[i][j]=0;
             for(int k=0;k<V;k++)
             {
                 n_m_z[i][j]=n_m_z[i][j]+static_cast<unsigned int>(n_m_t[i][k][j]); 
             }

        }

    }    
        
    for(int x=0;x<m;x++)
    {   double sum=0;
        for(int y=0;y<K;y++)
        {
		sum=sum+n_m_z[x][y];
        }
        for(int y=0;y<K;y++)
        {
                myfile<<round(n_m_z[x][y]/sum * 100) / 100<<"\t";
        }
        myfile<<endl;


    }
    myfile.close();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < K; j++) {
        *(n_m_z_temp+i*K+j)=n_m_z[i][j];
       
            
    }}
    std::ofstream outFile1("n_m_z"+std::to_string(r_iter)+".dat", std::ios::binary);
    outFile1.write(reinterpret_cast<const char*>(n_m_z_temp), m*K*sizeof(unsigned int));
    outFile1.close(); 
               
    for (int i = 0; i < K; i++) {
        for (long long int j = 0; j < V; j++) {
            //cout<<i<<"\t"<<j<<endl;
        *(n_z_t_temp+i*V+j)=n_z_t[i][j];      
    }}  

    std::ofstream outFile2("n_z_t"+std::to_string(r_iter)+".dat", std::ios::binary);
    outFile2.write(reinterpret_cast<char *>(n_z_t_temp), K * V * sizeof(unsigned int));
    outFile2.close();
        
    }
        

    std::ofstream outFile5("n_m_t.dat", std::ios::out | std::ios::binary);
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < V; j++) {
            for(int k=0;k<K;k++){
       outFile5.write(reinterpret_cast<const char*>(&n_m_t[i][j][k]), sizeof(short));
            }
    }}
    outFile5.close();    
    
    seedFile.close();
    return 0;
}
