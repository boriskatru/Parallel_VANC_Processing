#include <iostream>
#include <omp.h>
#include <math.h>
#include <vector>
#include <filesystem>
#include <chrono>
#include <fstream>
#include <string>
#include <numeric>
#include <iterator>
#include <algorithm>
#include "wait_bh.h"

#define THREADS 2


using namespace std;

inline double avrg(vector<double> vec, int start, int len) {
    //cout << reduce(vec.begin() + start, vec.begin() + start + len) / len << endl;
    return (reduce(vec.begin() + start, vec.begin() + start + len) / len);
}

inline int local_max(vector<double> vec, int start = 0,  int len = INT_MAX) {
    double max_ = DBL_MIN;    
    int max_pos = 0;
    min(len, (int)vec.size());
    int length = min(len, (int)vec.size());;
    for (int pos = start; pos < start + length; pos++)  
    {       
        if (vec[pos] > max_) {
            max_ = vec[pos];
            max_pos = pos;
        }       
    }
    return max_pos;
}
inline int local_min(vector<double> vec, int start = 0, int len = INT_MAX) {
    double min_ = DBL_MAX;
    int min_pos = 0;
    int length = min(len, (int)vec.size());
    for (int pos = start; pos < start + length; pos++)   
    {
        if (vec[pos] < min_) {
            min_ = vec[pos];
            min_pos = pos;
        }
    }
    return min_pos;
}

inline pair<int, double> find_nearest(vector<double> vec, double val) {
    int pos;
    double dev = DBL_MAX;
    int size = vec.size();
    
    for (int i = 0; i < size; i++) {
        if (abs(vec[i] - val) < abs(dev)) {
            dev = vec[i] - val;
            pos = i;
        }        
    }
    return make_pair(pos, dev);
}
inline  vector<vector<double>>  decimate_vanc(vector<vector<double>> data, int p_num) {
    vector<vector<double>>  data_new(3, vector<double>(p_num, 0));
    double range = abs(data[0][0] - data[0][data[0].size()-1]);   
    int step = data[0].size() / p_num;
  /*  if (!is_fw) reverse(data.begin(), data.end());*/

    for (int i = 0; i < p_num; i++) {
        //cout << i << endl;
        data_new[0][i] = avrg(data[0], i * step, step);
        data_new[1][i] = avrg(data[1], i * step, step);
        data_new[2][i] = avrg(data[2], i * step, step);
    }

    return data_new;
}
inline vector<vector<double>> avrg_curve(vector<vector<vector<double>>> data) {
    int VAC_cnt = data.size();
    int VAC_size = data[0][0].size();
    vector<vector<double>> avrg_vac(3, vector<double>(VAC_size, 0));
    for (int k = 0; k < VAC_size; k++) {
        for (int i = 0; i < VAC_cnt; i++) {
            avrg_vac[0][k] += data[i][0][k];
            avrg_vac[1][k] += data[i][1][k];
            avrg_vac[2][k] += data[i][2][k];
        }
        avrg_vac[0][k] = avrg_vac[0][k] / VAC_cnt;
        avrg_vac[1][k] = avrg_vac[1][k] / VAC_cnt;
        avrg_vac[2][k] = avrg_vac[2][k] / VAC_cnt;
    }
    return avrg_vac;
}

inline double interpolate(vector<vector<double>> vac, double target) {
    auto nearest = find_nearest(vac[0], target);
    int n_max, n_min;
    if (nearest.second == 0) return vac[1][nearest.first];
    if (nearest.second > 0) {
        n_max = nearest.first;
        n_min = max(n_max - 1,0);
    }
    else {
        n_min = nearest.first;
        n_max = min(n_min + 1, (int)vac[0].size()-1);
    }
    double dev = target - vac[0][n_min];
    double step= vac[0][n_max]- vac[0][n_min];
    //cout << step << "  " << dev << endl;
    if (step!=0)
        return (dev * vac[1][n_max] + (step - dev) * vac[1][n_min]) / 2 / step;
    else return (vac[1][n_max] +  vac[1][n_min]) / 2 ;
}
inline double mean_sq_deviation(vector<vector<double>> vac, vector<vector<double>> ref) {

    int VAC_size = vac[0].size();
    double deviation = 0;
    double tmp = 0;
    for (int k = 0; k < VAC_size; k++) {
        tmp = interpolate(vac, ref[0][k]) - ref[1][k];
        deviation += tmp * tmp;
    }
    return deviation/ VAC_size;
}
inline pair<vector<vector<vector<double>>>, vector<vector<vector<double>>>> separate_vancs(vector<vector<double>> data, double crit_dev,int phase_shift=20, int period=1000, int p_num = 100) {
    int min_1pos = local_min(data[0], 0, 9 * period / 10);
    if (min_1pos == 0) min_1pos = period;
    int max_1pos = local_max(data[0], 0, 9 * period / 10);
    if (max_1pos == 0) max_1pos = period;
    //cout << min_1pos << endl << max_1pos << endl;
    int start = min(min_1pos, max_1pos);
    bool fw_first = min_1pos < max_1pos;
    printf(" %i thread.  Start: %i      Start val: %f \n",  omp_get_thread_num(), start, data[0][start]);

    int VAC_cnt = floor((data[0].size() - start) / period);
    vector<vector<double>> bw_tmp(3, vector<double>(period / 2, 0));
    vector<vector<double>> fw_tmp(3, vector<double>(period / 2, 0));
    vector<vector<vector<double>>> bw_vancs(VAC_cnt, vector<vector<double>>(3, vector<double>(p_num, 0)));
    vector<vector<vector<double>>> fw_vancs(VAC_cnt, vector<vector<double>>(3, vector<double>(p_num, 0)));
    int offset_fw = (fw_first ? 0 : period / 2) + start;
    int offset_bw = (fw_first ? period / 2 : 0) + start;
    for (int i = 0; i < VAC_cnt; i++) {
        for (int k = 0; k < period / 2; k++) {
           // cout << "K= " << k << endl;
            fw_tmp[0][k] = data[0][period * i + k + offset_fw];
            fw_tmp[1][k] = data[1][period * i + k + offset_fw + phase_shift];
            fw_tmp[2][k] = data[2][period * i + k + offset_fw + phase_shift];
            bw_tmp[0][k] = data[0][period * i + k + offset_bw];
            bw_tmp[1][k] = data[1][period * i + k + offset_bw + phase_shift];
            bw_tmp[2][k] = data[2][period * i + k + offset_bw + phase_shift];
           // cout <<"K= "<< k << endl;
        }
        //cout << i << endl;
        reverse(bw_tmp[0].begin(), bw_tmp[0].end());
        reverse(bw_tmp[1].begin(), bw_tmp[1].end());
        reverse(bw_tmp[2].begin(), bw_tmp[2].end());
        //cout << i << endl;
        bw_vancs[i] = decimate_vanc(bw_tmp, p_num);
        fw_vancs[i] = decimate_vanc(fw_tmp, p_num);  
    }

    auto avrg_fw = avrg_curve(fw_vancs);
    auto avrg_bw = avrg_curve(bw_vancs);
    auto iter_bw = bw_vancs.begin();
    auto iter_fw = fw_vancs.begin();
    int i = 0;
    std::cout <<"fw_bw deviation: " << mean_sq_deviation(avrg_bw, avrg_fw) << endl;
    while (iter_bw != bw_vancs.end())   
    {
        
        if (mean_sq_deviation(*iter_bw, avrg_bw) > crit_dev)   bw_vancs.erase(iter_bw); 
        else  iter_bw++;        
    }
   
    while (iter_fw != fw_vancs.end())    
    {       
        if (mean_sq_deviation(*iter_fw, avrg_fw) > crit_dev) fw_vancs.erase(iter_fw);
        else iter_fw++;
    }
    
    return make_pair(fw_vancs, bw_vancs);
}
inline void print_avrg_VANCS(pair<vector<vector<vector<double>>>, vector<vector<vector<double>>>> VANCS) {
    auto avrg_fw = avrg_curve(VANCS.first);
    auto avrg_bw = avrg_curve(VANCS.second);
    ofstream  fw_vanc, bw_vanc;
    fw_vanc.open("fw_vanc.dat");
    bw_vanc.open("bw_vanc.dat");
    for (int i = 0; i < avrg_fw[0].size(); i++) {
        for (int k = 0; k < 3; k++) {
            fw_vanc << avrg_fw[k][i] << "  ";
            bw_vanc << avrg_bw[k][i] << "  ";
        }
        fw_vanc << endl;
        bw_vanc << endl;
    }
}




int main()
{
    omp_set_num_threads(THREADS);
    
    int start = 2;
    int stop = 99;
    int phase_shift = 20;
    int period = 1000;
    double crit_dev = 0.03;
    int p_num = 50;
    string name = "C:/Users/boris/OneDrive/Рабочий стол/17_33/VANC_";
    vector<vector<vector<double>>> fw_vancs, bw_vancs;
    Timer total_tmr;
    total_tmr.set_to_zero();
#pragma omp parallel for 
    for (int i = start; i < stop; i++)
    {
        Timer tmr;
        tmr.set_to_zero();
        int size = -1;
        string line;
        std::ifstream file;
        file.open(name + to_string(i)+".dat");       
        while (!file.eof())
        {
            std:getline(file, line);
            size++;
        }
        file.close();
        //printf("%i VAC doing by %i thread. Size of file %i \n", i, omp_get_thread_num(), size);

        file.open(name + to_string(i) + ".dat");
                
        vector<vector<double> > data(3, vector<double>(size, 0));

        for (int i = 0; i < size; i++) {
            file >> data[0][i];
            file >> data[1][i];
            file >> data[2][i];
        }
        auto VANCS = separate_vancs(data, crit_dev, phase_shift, period, p_num);
        //system("pause");
        file.close();
        fw_vancs.insert(fw_vancs.end(), VANCS.first.begin(), VANCS.first.end());
        bw_vancs.insert(bw_vancs.end(), VANCS.second.begin(), VANCS.second.end());
        printf("%i VAC done by %i thread. Time per file %f s  FW count: %i    BW count: %i\n", i, omp_get_thread_num(), tmr.get_loop_interval()/1000000
                                                                                                                , VANCS.first.size(), VANCS.second.size());
    }
    auto VANCS = make_pair(fw_vancs, bw_vancs);
    std::cout << "total time, s: " << total_tmr.get_full_interval()/1000000;
    print_avrg_VANCS(VANCS);
    
}
