#ifndef _ParallelTimer_
#define _ParallelTimer_
#include<mpi.h>
#include<iostream>
#include<fstream>
#include<map>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include<vector>

class ParallelTimer{

    public:
        ParallelTimer(){};
        void   createTimer(std::string);
        void   start(std::string);
        void   stop(std::string);
        double getStart(std::string str){ return time_array[str][0];};
        double getStop(std::string str){ return time_array[str][1];};
        double getTime(std::string str){ return time_array[str][2];};
        double getMaxTime(std::string str){ return time_array[str][3]; };
        double getMinTime(std::string str){ return time_array[str][4]; };
        double getAvgTime(std::string str){ return time_array[str][5]; };
        double getAccumulatedTime(std::string str){ return time_array[str][6];};
        double getAccumulatedMaxTime(std::string str){ return time_array[str][7];};
        double getAccumulatedMinTime(std::string str){ return time_array[str][8];};
        double getAccumulatedAvgTime(std::string str){ return time_array[str][9];};

        void   printTimers(); 
        void   printTimers(std::string); 
        void   printTimer(std::string str); 
        void   printTimerFull(std::string str); 



    private:
        //             0      1      2      3    4    5       6         7         8         9
        //double [10] start, stop, ltime,  max, min, avg, accum_time, acc_max, acc_min, accum_avg
	std::map<std::string, std::vector<double>> time_array;
	std::map<std::string, bool> checker;
        double getTime_cpu();
        int rank;        
};

#endif
