#include "../../scripts/globals.h"

int main(int argc, char* argv[])
{
    newWorkingDirectory="Demixer_cpp/input";
    
    set();
    set_variables();    
    cout<<m<<"\t"<<K<<endl;
    
    for(int j=0;j<V;j++)
    {
       snptype[j]=4;
       snppos[j]=0;
    }
        
    unsigned int a=rd();
    unsigned int b=rd();
    unsigned int c=rd();
    cout<<a<<"\t"<<b<<"\t"<<c<<endl;
    seed_seq seq{a,b,c};
    std::vector<std::uint64_t> seeds(1001);
    seq.generate(seeds.begin(), seeds.end());

    known=0;
    /*
    if(atoi(argv[2])==1)
        known=1;
    */
    if(atoi(argv[2])>0)
        known=atoi(argv[2]);
    
    if(atoi(argv[1])==1)
    {
      sprintf(buff,"%s","n_z_t.dat");

      latfile=fopen(buff,"r");


       for(int i=0;i<K;i++)
      {
          unsigned int *temp_K=new unsigned int [V]();
          fread(temp_K,sizeof(unsigned int),V,latfile);
          n_z_t[i]=temp_K;  

      }  

      fclose(latfile);


      sprintf(buff,"%s","n_m_t.dat");

      latfile=fopen(buff,"r");

      for(int i=0;i<m;i++){
          n_m_t[i]=new short *[V];

        for(int j=0;j<V;j++){
            short *temp_K=new short [K]();
            fread(temp_K,sizeof(short),K,latfile);
            n_m_t[i][j]=temp_K;}
      }
      fclose(latfile);

      sprintf(buff,"%s","n_m_z.dat");

      latfile=fopen(buff,"r");

      for(int i=0;i<m;i++)
      {
          unsigned int *temp_K=new unsigned int [K]();
          fread(temp_K,sizeof(unsigned int),K,latfile);
          n_m_z[i]=temp_K;  

      }
      fclose(latfile);  

        /*
        for(int i=0;i<5;i++){
        for(int j=0;j<8;j++){
            cout<<n_m_z[i][j]<<"\t";
        }
        cout<<"\n";
        }  

        */
        //cout<<"Reading Done for NMF \n";
        
    }
    else
            initialize(seeds[0]);

    //cout<<"dict done\n";
    construct_dict();
    //cout<<"dict done\n";
    fstream myfile;
    int no_iter=1001;

        for(ii=1;ii<no_iter;ii++){
            epochs(seeds[ii],ii);

            }

    unsigned int n_z_t_temp[K][V];
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < V; j++) {
            n_z_t_temp[i][j]=n_z_t[i][j];


        }}    

    std::ofstream outFile;
    outFile.open("../output/n_z_t.dat", std::ios::binary);
    outFile.write(reinterpret_cast<const char*>(n_z_t_temp), sizeof(n_z_t_temp));
    outFile.close();

    unsigned int n_m_z_temp[m][K];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < K; j++) {
            n_m_z_temp[i][j]=n_m_z[i][j];


        }}
    outFile.open("../output/n_m_z.dat", std::ios::binary);
    outFile.write(reinterpret_cast<const char*>(n_m_z_temp), sizeof(n_m_z_temp));
    outFile.close();  

        return 0;
    }
