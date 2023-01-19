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
#include <windows.h>
#include "wait_bh.h"

#define THREADS 2


using namespace std;

inline float avrg(vector<float>& vec, int start, int len) {
    //cout << reduce(vec.begin() + start, vec.begin() + start + len) / len << endl;
    return (reduce(vec.begin() + start, vec.begin() + start + len) / len);
}

inline int local_max(vector<float>& vec, int start = 0,  int len = INT_MAX) {
    float max_ = DBL_MIN;    
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
inline int local_min(vector<float>& vec, int start = 0, int len = INT_MAX) {
    float min_ = DBL_MAX;
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

inline pair<int, float> find_nearest(vector<float>& vec, float val) {
    int pos;
    float dev = DBL_MAX;
    int size = vec.size();
    
    for (int i = 0; i < size; i++) {
        if (abs(vec[i] - val) < abs(dev)) {
            dev = vec[i] - val;
            pos = i;
        }        
    }
    return make_pair(pos, dev);
}
inline  vector<vector<float>>  decimate_vanc(vector<vector<float>>& data, int p_num) {
    vector<vector<float>>  data_new(3, vector<float>(p_num, 0));
    float range = abs(data[0][0] - data[0][data[0].size()-1]);   
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
inline vector<vector<float>> avrg_curve(vector<vector<vector<float>>>& data) {
    int VAC_cnt = data.size();
    int VAC_size = data[0][0].size();
    vector<vector<float>> avrg_vac(3, vector<float>(VAC_size, 0));
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

inline float interpolate(vector<vector<float>>& vac, float target) {
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
    float dev = target - vac[0][n_min];
    float step= vac[0][n_max]- vac[0][n_min];
    //cout << step << "  " << dev << endl;
    if (step!=0)
        return (dev * vac[1][n_max] + (step - dev) * vac[1][n_min]) / 2 / step;
    else return (vac[1][n_max] +  vac[1][n_min]) / 2 ;
}
inline float mean_sq_deviation(vector<vector<float>>& vac, vector<vector<float>>& ref) {

    int VAC_size = vac[0].size();
    float deviation = 0;
    float tmp = 0;
    for (int k = 0; k < VAC_size; k++) {
        tmp = interpolate(vac, ref[0][k]) - ref[1][k];
        deviation += tmp * tmp;
    }
    return deviation/ VAC_size;
}
inline int reject_bad_vancs(vector<vector<vector<float>>>& vancs, vector<vector<float>>& ref, float crit_dev) {
    auto iter= vancs.begin();
    int i = 0;
    while (iter != vancs.end()) {

        if (mean_sq_deviation(*iter, ref) > crit_dev) {
            vancs.erase(iter); 
            i++;
        }
        else {
            iter++;
            //cout <<"std dev:" << mean_sq_deviation(*iter, ref) << endl;
        }
    }

    return i;
}
inline pair<vector<vector<vector<float>>>, vector<vector<vector<float>>>> process_vancs(vector<vector<float>>& data, float crit_dev, int phase_shift_cur = 20, int phase_shift_n = 15, int period=1000, int p_num = 100) {
    int min_1pos = local_min(data[0], 0, 9 * period / 10);
    if (min_1pos == 0) min_1pos = period;
    int max_1pos = local_max(data[0], 0, 9 * period / 10);
    if (max_1pos == 0) max_1pos = period;
    //cout << min_1pos << endl << max_1pos << endl;
    int start = min(min_1pos, max_1pos);
    bool fw_first = min_1pos < max_1pos;
    printf(" %i thread.  Start: %i      Start val: %f \n",  omp_get_thread_num(), start, data[0][start]);

    int VAC_cnt = floor((data[0].size() - start) / period);
    vector<vector<float>> bw_tmp(3, vector<float>(period / 2, 0));
    vector<vector<float>> fw_tmp(3, vector<float>(period / 2, 0));
    vector<vector<vector<float>>> bw_vancs(VAC_cnt, vector<vector<float>>(3, vector<float>(p_num, 0)));
    vector<vector<vector<float>>> fw_vancs(VAC_cnt, vector<vector<float>>(3, vector<float>(p_num, 0)));
    int offset_fw = (fw_first ? 0 : period / 2) + start;
    int offset_bw = (fw_first ? period / 2 : 0) + start;
    for (int i = 0; i < VAC_cnt; i++) {
        for (int k = 0; k < period / 2; k++) {
           // cout << "K= " << k << endl;
            fw_tmp[0][k] = data[0][period * i + k + offset_fw];
            fw_tmp[1][k] = data[1][period * i + k + offset_fw + phase_shift_cur];
            fw_tmp[2][k] = data[2][period * i + k + offset_fw + phase_shift_n];
            bw_tmp[0][k] = data[0][period * i + k + offset_bw];
            bw_tmp[1][k] = data[1][period * i + k + offset_bw + phase_shift_cur];
            bw_tmp[2][k] = data[2][period * i + k + offset_bw + phase_shift_n];
           // cout <<"K= "<< k << endl;
        }
        //cout << i << endl;
        std::reverse(bw_tmp[0].begin(), bw_tmp[0].end());
        std::reverse(bw_tmp[1].begin(), bw_tmp[1].end());
        std::reverse(bw_tmp[2].begin(), bw_tmp[2].end());
        //cout << i << endl;
        bw_vancs[i] = decimate_vanc(bw_tmp, p_num);
        fw_vancs[i] = decimate_vanc(fw_tmp, p_num);  
    }

    auto avrg_fw = avrg_curve(fw_vancs);
    auto avrg_bw = avrg_curve(bw_vancs);

    std::cout <<"fw_bw deviation: " << mean_sq_deviation(avrg_bw, avrg_fw) << endl;

    reject_bad_vancs(fw_vancs, avrg_fw, crit_dev);
    reject_bad_vancs(bw_vancs, avrg_bw, crit_dev);

    avrg_fw = avrg_curve(fw_vancs);
    avrg_bw = avrg_curve(bw_vancs);
    std::cout << "fw_bw deviation: " << mean_sq_deviation(avrg_bw, avrg_fw) << endl;

    return make_pair(fw_vancs, bw_vancs);
}
inline void print_avrg_VANCS(pair<vector<vector<vector<float>>>, vector<vector<vector<float>>>>& VANCS) {
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
    setlocale(LC_ALL, "Russian");
    SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS);
    SetThreadPriority(GetCurrentProcess(), THREAD_PRIORITY_TIME_CRITICAL);
    cout << GetPriorityClass(GetCurrentProcess()) << endl;

    cout << "Введите номер первой кривой:" << endl;
    int start = 1;
    cin >> start;

    cout << "Введите номер последней кривой:" << endl;
    int stop =10;
    cin >> stop;

    //cout << "Введите максимальное количество точек в файле:" << endl;
    int cnt = 2000000;
    //cin >> cnt;

    //cout << "Введите сдвиг фазы напряжения, в точках:" << endl;
    int phase_shift_cur = 25;
    //cin >> phase_shift_cur;

    //cout << "Введите сдвиг фазы шума, в точках:" << endl;
    int phase_shift_n = 10;
    //cin >> phase_shift_n;
  
    cout << "Введите желаемое количество точек на  обработанной кривой:" << endl;
    int p_num = 50;
    cin >> p_num;
   
    cout << "Введите критическое отклонение:" << endl;
    float crit_dev = 0.03;
    cin >> crit_dev;

    int period = 1000;
    string filename = "VANC_";
    string filetype = ".bin";
    vector<vector<vector<float>>> fw_vancs, bw_vancs;
    Timer total_tmr;
    total_tmr.set_to_zero();
#pragma omp parallel for 
    for (int i = start; i < stop; i++)
    {

        Timer tmr;
        tmr.set_to_zero();
        


        

        FILE* fileV;
        FILE* fileA;
        FILE* fileN;

        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << filename + to_string(i)  + filetype << endl;
        cout << "open voltage file: " << fopen_s(&fileV, (filename + to_string(i) + "V" + filetype).data(), "rb") << endl;
        cout << "open current file: " << fopen_s(&fileA, (filename + to_string(i) + "A" + filetype).data(), "rb") << endl;
        cout << "open noise file: " << fopen_s(&fileN, (filename + to_string(i) + "N" + filetype).data(), "rb") << endl;
    
        std::vector<float> Vbuf(cnt); // underlying storage of std::vector is also an array
        std::vector<float> Abuf(cnt); // underlying storage of std::vector is also an array
        std::vector<float> Nbuf(cnt); // underlying storage of std::vector is also an array
        std::size_t size = std::fread(Vbuf.data(), sizeof Vbuf[0], Vbuf.size(), fileV);
        std::fread(Abuf.data(), sizeof Abuf[0], Abuf.size(), fileA);
        std::fread(Nbuf.data(), sizeof Nbuf[0], Nbuf.size(), fileN);
        vector<vector<float> > data(3, vector<float>(size, 0));

        for (int i = 0; i < size; i++) {
            data[0][i] = Vbuf[i];
            data[1][i] = Abuf[i];
            data[2][i] = Nbuf[i];
        }
        fclose(fileV);
        fclose(fileA);
        fclose(fileN);

        auto VANCS = process_vancs(data, crit_dev, phase_shift_cur, phase_shift_n, period, p_num);
        //system("pause");
       
        fw_vancs.insert(fw_vancs.end(), VANCS.first.begin(), VANCS.first.end());
        bw_vancs.insert(bw_vancs.end(), VANCS.second.begin(), VANCS.second.end());
        printf("%i VAC done by %i thread. Time per file %f s  FW count: %i    BW count: %i\n", i, omp_get_thread_num(), tmr.get_loop_interval()/1000000
                                                                                                                , VANCS.first.size(), VANCS.second.size());
    }
    auto avrg_fw = avrg_curve(fw_vancs);
    reject_bad_vancs(fw_vancs, avrg_fw, crit_dev);
    auto avrg_bw = avrg_curve(bw_vancs);
    reject_bad_vancs(bw_vancs, avrg_bw, crit_dev);
    auto VANCS = make_pair(fw_vancs, bw_vancs);
    
    print_avrg_VANCS(VANCS);
    std::cout << "Done! Total time, s: " << total_tmr.get_full_interval() / 1000000;
    getchar(); getchar();
    
}
