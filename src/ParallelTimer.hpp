#ifndef _ParallelTimer_
#define _ParallelTimer_
#include<mpi.h>
#include<iostream>
#include<map>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>

using namespace std;

class ParallelTimer{

    public:
        ParallelTimer(){};
        void   createTimer(string);
        void   start(string);
        void   stop(string);
        double getStart(string str){ return time_array[str][0];};
        double getStop(string str){ return time_array[str][1];};
        double getTime(string str){ return time_array[str][2];};
        double getMaxTime(string str){ return time_array[str][3]; };
        double getMinTime(string str){ return time_array[str][4]; };
        double getAvgTime(string str){ return time_array[str][5]; };
        double getAccumulatedTime(string str){ return time_array[str][6];};
        double getAccumulatedMaxTime(string str){ return time_array[str][7];};
        double getAccumulatedMinTime(string str){ return time_array[str][8];};
        double getAccumulatedAvgTime(string str){ return time_array[str][9];};

        void   printTimer(string str); 
        void   printTimerFull(string str); 



    private:
        //             0      1      2      3    4    5       6         7         8         9
        //double [10] start, stop, ltime,  max, min, avg, accum_time, acc_max, acc_min, accum_avg
        map<string, double[10]> time_array;
        map<string, bool> checker;
        double getTime_cpu();
        int rank;        
};

#endif
