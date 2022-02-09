#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include "ParallelTimer.hpp"

int main(int argc, char *argv[])
{

    MPI_Init(&argc,&argv);
    ParallelTimer lala;
    lala.createTimer("main_code");
    lala.createTimer("only_loop");
    lala.createTimer("prueba_loop");
    lala.start("main_code");

    //Init shit 
    double pre=0.0,suma=0.0,post=0.0;
    for(int i=0;i< 10000;i++)
        pre= pre + cos(2.0*i);

    for(int k=0;k<5;k++){


        lala.start("only_loop");
        for(int i=0;i< 10000;i++)
            suma= suma + sin(2.0*i);

        lala.stop("only_loop");


        //Init shit 
        for(int i=0;i< 10000;i++)
            post= post + cos(2.0*i);


        //lala.stop("prueba_loop");


        lala.printTimerFull("only_loop");

    }
        lala.stop("main_code");


        lala.printTimerFull("main_code");

    cout<<"pre "<<pre<<" suma "<<suma<<" post "<<post<<endl;




//    cout<<"lala timer start"<<lala.getStart("main_code")<<endl;
//    cout<<"lala timer stop "<<lala.getStop("main_code")<<endl;
//    cout<<"lala timer diff "<<lala.getTime("main_code")<<endl;
//
//    cout<<"lala timer start"<<lala.getStart("only_loop")<<endl;
//    cout<<"lala timer stop "<<lala.getStop("only_loop")<<endl;
//    cout<<"lala timer diff "<<lala.getTime("only_loop")<<endl;


    MPI_Finalize();

}
