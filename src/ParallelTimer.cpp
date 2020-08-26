#include "ParallelTimer.hpp"

double ParallelTimer::getTime_cpu()
{
     struct timeval now;

     if (gettimeofday(&now, (struct timezone *) 0) != 0) {
         exit(1);
     }

     return(((double) now.tv_sec + (double) now.tv_usec / 1000000) +1.0);

}

void ParallelTimer::createTimer(string str)
{
   time_array[str][0] = 0.0;
   time_array[str][1] = 0.0;
   time_array[str][2] = 0.0;
   time_array[str][3] = 0.0;
   time_array[str][4] = 0.0;
   time_array[str][5] = 0.0;
   time_array[str][6] = 0.0;
   time_array[str][7] = 0.0;
   time_array[str][8] = 0.0;
   time_array[str][9] = 0.0;

   checker[str] = false;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
}

void ParallelTimer::start(string str)
{
    time_array[str][0] = getTime_cpu();    
    checker[str] = true;
}

void ParallelTimer::stop(string str)
{

    MPI_Barrier(MPI_COMM_WORLD);
    time_array[str][1] = getTime_cpu();
    time_array[str][2] = time_array[str][1] - time_array[str][0];      

    if(checker[str] == false)
    {
        cout<<"WARNING: Trying to stop the timer '"<<str<<"' that does not have a start() call"<<endl; 
    }

    int nsize;
    MPI_Comm_size(MPI_COMM_WORLD,&nsize);
    double maxtime,mintime,avgtime,ltime;

    ltime = time_array[str][2];

    MPI_Allreduce(&ltime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&ltime,&mintime,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&ltime,&avgtime,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    time_array[str][3] = maxtime;
    time_array[str][4] = mintime;
    time_array[str][5] = avgtime/nsize;


    time_array[str][6] = time_array[str][6] + ltime; 

    double acctime = time_array[str][6];  

    MPI_Allreduce(&acctime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&acctime,&mintime,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&acctime,&avgtime,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    time_array[str][7] = maxtime;
    time_array[str][8] = mintime;
    time_array[str][9] = avgtime/nsize;


    checker[str] = false;
}

void ParallelTimer::printTimerFull(string str)
{
    if( rank == 0) {
        cout<<endl;
        cout<<"Timer : "<<str<<endl;
        cout<<"------------------------------"<<endl;
        cout<<"Last Time        :"<<time_array[str][2]<<endl;
        cout<<"Last Time Max    :"<<time_array[str][3]<<endl;
        cout<<"Last Time Min    :"<<time_array[str][4]<<endl;
        cout<<"Last Time Avg    :"<<time_array[str][5]<<endl;
       
        cout<<"Accumulated Time :"<<time_array[str][6]<<endl;
        cout<<"Accum. Time Max  :"<<time_array[str][7]<<endl;
        cout<<"Accum. Time Min  :"<<time_array[str][8]<<endl;
        cout<<"Accum. Time Avg  :"<<time_array[str][9]<<endl;
        cout<<"------------------------------"<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void ParallelTimer::printTimer(string str)
{
    if( rank == 0) {
        cout<<endl;
        cout<<"Timer : "<<str<<endl;
        cout<<"------------------------------"<<endl;
        cout<<"Last Time Max     :"<<time_array[str][3]<<endl;
        cout<<"Accum. Time Max   :"<<time_array[str][7]<<endl;
        cout<<"------------------------------"<<endl;
   }
    MPI_Barrier(MPI_COMM_WORLD);
}

void ParallelTimer::printTimers()
{

    if( rank == 0) {
        map<string, double[10]>::iterator it_d;
        for (it_d = time_array.begin(); it_d != time_array.end(); it_d++)
        {
            cout<<endl;
            cout<<"Timer : "<<it_d->first<<endl;
            cout<<"------------------------------"<<endl;
            cout<<"Last Time Max     :"<<it_d->second[3]<<endl;
            cout<<"Accum. Time Max   :"<<it_d->second[7]<<endl;
            cout<<"------------------------------"<<endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

void ParallelTimer::printTimers(string str)
{




    if( rank == 0) {
        char filename[100];
        sprintf(filename,"%s",str.c_str());

        ofstream outfile;

        outfile.open(filename,std::ofstream::app);


        map<string, double[10]>::iterator it_d;
        for (it_d = time_array.begin(); it_d != time_array.end(); it_d++)
        {
            outfile<<endl;
            outfile<<"Timer : "<<it_d->first<<endl;
            outfile<<"------------------------------"<<endl;
            outfile<<"Last Time Max     :"<<it_d->second[3]<<endl;
            outfile<<"Accum. Time Max   :"<<it_d->second[7]<<endl;
            outfile<<"------------------------------"<<endl;
        }

        outfile.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

