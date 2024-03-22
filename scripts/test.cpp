#include "globals.h"


int main(int argc,char* argv[])
{

newWorkingDirectory="finaloutput/Cryptic_rerun/";
cout<<newWorkingDirectory; 
int no_iter;
test_set();

construct_dict();
cout<<"dict done\n";
fstream myfile;
if(m==1)
    no_iter=100;
else
    no_iter=500;
    
unsigned int a,b,c;
    std::ofstream seedFile("../test/output.txt");
    for(int r_iter=0;r_iter<1;r_iter++)
    {
    cout<<"iter_num:   "<<r_iter<<endl;
    //a=2989879115;
    //b=1980914546;
    //c=1229027252;
    a=1985875050;
    b=3586869036;
    c=2086457268;
    //a=rd();
    //b=rd();
    //c=rd();
    seedFile<<a<<"\t"<<b<<"\t"<<c<<endl;
    seed_seq seq{a,b,c};
    std::vector<std::uint64_t> seeds(1001);
    seq.generate(seeds.begin(), seeds.end());
  

    test_initialize(seeds[0]);    
    unsigned int *n_m_z_temp= new unsigned int [m*K]();  
    unsigned int *n_z_t_temp= new unsigned int [K*V]();    
   
    
    time_t start = time(&start);
    for(ii=1;ii<no_iter;ii++){
    	epochs(seeds[ii],ii);
        
       if(ii%10==0)
	{
		myfile.open("../test/n_m_z"+std::to_string(r_iter)+".txt",fstream::out);

    		for(int x=0;x<m;x++)
		{        double sum=0;
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
                        std::ofstream outFile1("../test/n_m_z"+std::to_string(r_iter)+".dat", std::ios::binary);
                        outFile1.write(reinterpret_cast<const char*>(n_m_z_temp), m*K*sizeof(unsigned int));
                          outFile1.close();

	} 
	}
    time_t end = time(&end);
    cout<<"Total time "<<static_cast<long int>(end-start);


        
myfile.open("../test/n_m_z"+std::to_string(r_iter)+".txt",fstream::out);

    for(int x=0;x<m;x++)
{        double sum=0;
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
std::ofstream outFile1("../test/n_m_z"+std::to_string(r_iter)+".dat", std::ios::binary);
outFile1.write(reinterpret_cast<const char*>(n_m_z_temp), m*K*sizeof(unsigned int));
  outFile1.close(); 
        
           
    for (int i = 0; i < K; i++) {
        for (long long int j = 0; j < V; j++) {
            //cout<<i<<"\t"<<j<<endl;
        *(n_z_t_temp+i*V+j)=n_z_t[i][j];      
    }}  

std::ofstream outFile2("../test/n_z_t"+std::to_string(r_iter)+".dat", std::ios::binary);
outFile2.write(reinterpret_cast<char *>(n_z_t_temp), K * V * sizeof(unsigned int));
outFile2.close();
    
    }
        
std::ofstream outFile5("../test/n_m_t.dat", std::ios::out | std::ios::binary);
    
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
