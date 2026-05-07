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

// ============================================================================
// ПАРАЛЛЕЛЬНАЯ ОБРАБОТКА ВОЛЬТ-АМПЕРНЫХ ХАРАКТЕРИСТИК (VANC)
// Туннельный микроскоп STM (Scanning Tunneling Microscope) регистрирует
// туннельный ток при различных напряжениях смещения. Это позволяет получить
// информацию о локальной плотности состояний образца.
//
// Программа обрабатывает измерения в двух направлениях:
//   - Forward (FW): прямое сканирование напряжений
//   - Backward (BW): обратное сканирование напряжений
//
// Для каждого измерения регистрируются три канала:
//   - V (Voltage): напряжение смещения
//   - A (Amplitude/Current): туннельный ток
//   - N (Noise): шум туннельного тока
// ============================================================================

#define THREADS 4
#define PI 3.14159f
#define EMPTY_VANC -1000  // Маркер для отклоненных/невалидных VANC
using namespace std;

const std::string MAIN_FOLDER = "C:/Users/Tunnel Noise/Desktop/STM/";
const std::string SETTINGS_FOLDER = "Settings/";
const std::string LAST_VANC_FILE_NAME = SETTINGS_FOLDER + "LAST_VANC.txt";
const std::string SESSION_FILE_NAME = SETTINGS_FOLDER + "SESSION_DATA.txt";
const std::string VANC_LIST_NAME = "\\vac_list.txt";
const std::string PARAM_FILE_NAME= "VAC_PROCESSING.txt";

/// Читает последнюю использованную директорию с VANC из конфига
inline string ReadLastVANCDiretory()
{
    ifstream file;
    string tmp1, tmp2, ans;
    file.open(MAIN_FOLDER + LAST_VANC_FILE_NAME, std::ios::in);
    file >> tmp1;
    file >> tmp2;
    ans = tmp1 + " " + tmp2;
    file.close();
    return ans;
}

/// Читает директорию текущей сессии из конфига
inline string ReadSessionDirectory()
{
    ifstream file;
    string tmp1, tmp2,ans;
    file.open(MAIN_FOLDER + SESSION_FILE_NAME, std::ios::in);
    file >> tmp1;
    file >> tmp2;
    ans = tmp1 + " " + tmp2;
    file.close();
    return ans;
}

/// Добавляет имя обработанного файла в список сессии
inline void AddFileToSession(string filename, string path = VANC_LIST_NAME )
{
    string session_folder = ReadSessionDirectory();
    cout << "Filename " << filename << " added to session path " << session_folder + path << endl;
    ofstream list;
    list.open(session_folder + path, std::ios_base::app);
    list << filename << endl;
    list.close();
}

/// Вычисляет среднее арифметическое элементов вектора в диапазоне [start, start+len)
/// Используется для децимации (понижения частоты дискретизации) VANC данных
inline float avrg(vector<float>& vec, int start, int len) {
    //cout << reduce(vec.begin() + start, vec.begin() + start + len) / len << endl;
    return (reduce(vec.begin() + start, vec.begin() + start + len) / len);
}

/// Находит индекс максимального элемента в диапазоне [start, start+len)
/// Используется для поиска характерных точек в VANC (например, пиков)
inline int local_max(vector<float>& vec, int start = 0,  int len = INT_MAX) {
    float max_ = FLT_MIN;
    int max_pos = 0;
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

/// Находит индекс минимального элемента в диапазоне [start, start+len)
/// Используется для поиска характерных точек в VANC (например, впадин)
inline int local_min(vector<float>& vec, int start = 0, int len = INT_MAX) {
    float min_ = FLT_MAX;
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

/// Находит ближайший элемент в векторе к заданному значению val
/// Возвращает пару: (индекс, отклонение от target)
/// Используется при интерполяции VANC на равномерную сетку напряжений
inline pair<int, float> find_nearest(vector<float>& vec, float val) {
    int pos = 0;
    float dev = FLT_MAX;
    int size = vec.size();
    
    for (int i = 0; i < size; i++) {
        if (abs(vec[i] - val) < abs(dev)) {
            dev = vec[i] - val;
            pos = i;
        }        
    }
    return make_pair(pos, dev);
}

/// Находит две ближайшие точки к значению val (используется редко, в основном закомментирована)
inline pair<int, int> find_two_nearest(vector<float>& vec, float val) {
    int pos1 = 0, pos2 = 0;
    float dev = FLT_MAX;
    int size = vec.size();
    for (int i = 0; i < size; i++) {
        if (abs(vec[i] - val) < abs(dev)) {
            dev = vec[i] - val; 
            pos1 = i;
        }       
    }
    for (int i = 0; i < size; i++) {
        if ((abs(vec[i] - val) < abs(dev)) && (i != pos1)) {            
            dev = vec[i] - val;
            pos2 = i;
        }

    }
    return make_pair(pos1, pos2);
}

/// ДЕЦИМАЦИЯ VANC - понижение частоты дискретизации с усреднением
/// Входные данные: data[3] - три канала (V, A, N), каждый содержит множество точек
/// Выход: data_new[3] - три канала по p_num точкам каждый
/// 
/// Алгоритм: исходные данные разбиваются на p_num блоков одинакового размера,
/// каждый блок усредняется в одну точку. Это позволяет:
/// - Снизить шум в данных
/// - Уменьшить объем данных для дальнейшей обработки
/// - Упростить анализ структур в VANC
inline  vector<vector<float>>  decimate_vanc(vector<vector<float>>& data, int p_num) {
    vector<vector<float>>  data_new(3, vector<float>(p_num, 0));
    int step = data[0].size() / p_num;

    for (int i = 0; i < p_num; i++) {
        //cout << i << endl;
        data_new[0][i] = avrg(data[0], i * step, step);
        data_new[1][i] = avrg(data[1], i * step, step);
        data_new[2][i] = avrg(data[2], i * step, step);
    }

    return data_new;
}

/// ВЫРАВНИВАНИЕ VANC ПО НАПРЯЖЕНИЮ (не используется в основном цикле)
/// Преобразует VANC так, чтобы ось напряжения была равномерно расположена от минимума до максимума
/// !!!!!!!! В РАЗРАБОТКЕ !!!!!!!!
inline  vector<vector<float>>  equalize_vanc_step(vector<vector<float>>& data, int p_num) {

    vector<vector<float>>  data_new(3, vector<float>(p_num, 0));
    vector<float>  step_axis(p_num, 0);

    float range = data[0][data[0].size() - 1] - data[0][0];
    int step = static_cast<int>(range / p_num);

    for (int i = 0; i < p_num; i++) {
        data_new[0][i] = data[0][0] + step * i;
        data_new[1][i] = avrg(data[1], i * step, step);
        data_new[2][i] = avrg(data[2], i * step, step);
    }

    return data_new;
}

/// УСРЕДНЕНИЕ КРИВЫХ - вычисляет среднюю VANC из набора множественных измерений
/// Входные данные: data[n][3][m] - n кривых, каждая с 3 каналами и m точками
/// Выход: avrg_vac[3][m] - усредненная VANC с 3 каналами и m точками
///
/// Это является основной операцией для:
/// - Снижения влияния случайного шума на отдельные измерения
/// - Выявления воспроизводимых структур в VANC
/// - Улучшения соотношения сигнал/шум
inline vector<vector<float>> avrg_curve(vector<vector<vector<float>>>& data) {
    int VAC_cnt = data.size();
    if (VAC_cnt == 0) { 
        vector<vector<float>> null;
        return null;
    }
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

/// ЛИНЕЙНАЯ ИНТЕРПОЛЯЦИЯ VANC
/// Находит значение в канале channel при напряжении target, используя линейную интерполяцию
/// между ближайшими точками измерения.

/// Входные параметры:
///   vac: VANC с регулярной сеткой напряжений (channel 0)
///   target: требуемое напряжение для интерполяции
///   channel: какой канал интерполировать (0=напряжение, 1=ток, 2=шум)
///
/// Используется при:
/// - Перепроецировании VANC на стандартную сетку напряжений
/// - Сравнении VANC с разными шагами дискретизации
/// - Вычислении среднеквадратичного отклонения между VANC
inline float interpolate(vector<vector<float>>& vac, float target, int channel = 1) {
    auto nearest = find_nearest(vac[0], target);
    
    int n_max, n_min;
    if (nearest.second == 0) return vac[channel][nearest.first];
    if (nearest.second > 0) {
        n_max = min(nearest.first, (int)vac[0].size() - 1);
        n_min = max(n_max - 1,0);
    }
    else {
        n_min = nearest.first;
        n_max = min(n_min + 1, (int)vac[0].size()-1);
    }
    float dev = target - vac[0][n_min];
    float step= vac[0][n_max]- vac[0][n_min];
    
    if (step!=0)
        return (dev * vac[channel][n_max] + (step - dev) * vac[channel][n_min])  / step;
    else return (vac[channel][n_max] +  vac[channel][n_min]) / 2 ;
   /* auto nearests = find_two_nearest(vac[0], target);
    float step = vac[0][nearests.first] - vac[0][nearests.second];
    float dev = target - vac[0][nearests.first];
    if (step != 0)
        return (dev * vac[channel][nearests.first] + (step - dev) * vac[channel][nearests.second])  / step;
    else return (vac[channel][nearests.first] +  vac[channel][nearests.second]) / 2 ;*/

}

/// УСРЕДНЕНИЕ С ИНТЕРПОЛЯЦИЕЙ НА РАВНОМЕРНУЮ СЕТКУ
/// Выполняет более сложное усреднение, учитывающее что различные VANC могут иметь
/// разные "точки приложения" на оси напряжения.

/// Алгоритм:
/// 1. Вычисляет среднюю VANC из всех входных кривых
/// 2. Создает равномерную сетку напряжений от мин до макс со step = range/p_num
/// 3. Для каждой точки каждой входной VANC находит ближайший бин в сетке
/// 4. Интерполирует значение на центр бина и добавляет к этому бину
/// 5. Нормирует каждый бин по количеству попадших в него точек

/// Это более надёжный метод, чем простое усреднение, так как он корректирует
/// небольшие смещения между VANC из-за дрейфа напряжения смещения.
inline vector<vector<float>> avrg_interpl_curve(vector<vector<vector<float>>>& data, int p_num) {

    int VAC_cnt = data.size();    
    auto avrg_ = avrg_curve(data);
        
    float range = avrg_[0][avrg_[0].size() - 1] - avrg_[0][0];
    float step = range / p_num;
    cout << "Range:" << range << endl;
    cout << "Step:" << step << endl;
    vector<float>  step_axis(p_num, 0);
    vector<int> local_count(p_num, 0);
    for (int i = 0; i < p_num; i++) {
        step_axis[i] = avrg_[0][0] + step * i;        
    }
    int pos;
    cout << "VAC count: " << VAC_cnt << endl;
    vector<vector<float>> avrg_vac(3, vector<float>(p_num, 0));
    for (int k = 0; k < data[0][0].size(); k++) {
        for (int i = 0; i < VAC_cnt; i++) {   
            pos = find_nearest(step_axis, data[i][0][k]).first;
            local_count[pos]++;
            
            avrg_vac[0][pos] += data[i][0][k];
            avrg_vac[1][pos] += data[i][1][k];
            avrg_vac[2][pos] += data[i][2][k];
        }        
    }
    for (int k = 0; k < p_num; k++) {
        avrg_vac[0][k] = avrg_vac[0][k] / local_count[k];
        avrg_vac[1][k] = avrg_vac[1][k] / local_count[k];
        avrg_vac[2][k] = avrg_vac[2][k] / local_count[k];
    }
    return avrg_vac;
}

/// СРЕДНЕКВАДРАТИЧНОЕ ОТКЛОНЕНИЕ МЕЖДУ ДВУМЯ VANC
/// Вычисляет качество соответствия тестовой ВАНC (vac) к эталонной (ref)

/// Формула: MSD = (1/N) * Σ(interpolate(vac, ref[0][k]) - ref[1][k])²

/// Используется для:
/// - Определения "плохих" VANC (с большим шумом или артефактами)
/// - Фильтрации VANC по критерию quality
/// - Оптимизации параметров измерения

/// Отклонение интерпретируется как:
/// - Низкое (< 0.01): высокое качество, VANC надежна
/// - Среднее (0.01-0.1): приемлемое качество
/// - Высокое (> 0.1): VANC заметно отличается от типичной, может быть отклонена
inline float mean_sq_deviation(vector<vector<float>>& vac, vector<vector<float>>& ref) {

    int VAC_size = ref[0].size();
    
    float deviation = 0;
    float tmp = 0;
    for (int k = 0; k < VAC_size; k++) {
        tmp = interpolate(vac, ref[0][k]) - ref[1][k];
        deviation += tmp * tmp;       
    }
    return (deviation / static_cast<float>(VAC_size));
}

/// ОТКЛОНЕНИЕ ПЛОХИХ VANC - удаляет VANC, которые слишком отличаются от эталонной
/// 
/// Входные параметры:
///   vancs: вектор VANC для фильтрации
///   ref: эталонная VANC для сравнения
///   crit_dev: критическое среднеквадратичное отклонение (порог отклонения)

/// Алгоритм: перебирает все VANC, вычисляет их отклонение от ref и удаляет те,
/// у которых отклонение превышает crit_dev. Это помогает избавиться от:
/// - Измерений с аномально высоким шумом
/// - VANC с артефактами (скачки контакта, помехи)
/// - Выбросов из-за наконечника микроскопа

/// Возвращает: количество удаленных VANC
inline int reject_bad_vancs(vector<vector<vector<float>>>& vancs, vector<vector<float>>& ref, float crit_dev) {
    auto iter= vancs.begin();
    int i = 0;
    
    while (iter != vancs.end()) {
        //cout << i << "   :  " << mean_sq_deviation(*iter, ref) << endl;
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

/// ОСНОВНАЯ ФУНКЦИЯ ОБРАБОТКИ VANC
/// Выполняет комплексную обработку измеренных данных из трех каналов (V, A, N)

/// Входные параметры:
///   data: измеренные данные (3 канала x множество точек)
///   crit_dev: критерий для отклонения плохих VANC
///   phase_shift_cur: сдвиг фазы для канала тока (в точках)
///   phase_shift_n: сдвиг фазы для канала шума (в точках)
///   period: количество точек в одном периоде биполярного сканирования (V+A+ и V-A-)
///   p_num: требуемое количество точек в выходной VANC после децимации
///   VAC_num: номер VANC (для вывода)
///   tmr: таймер для отсчета времени

/// Алгоритм обработки:
/// 1. РАЗБИЕНИЕ НА ПОЛУПЕРИОДЫ
///    - Данные содержат "биполярные" циклы: полупериод 1 (напряжение растет),
///      полупериод 2 (напряжение падает). Эта структура соответствует треугольному
///      сигналу развертки напряжений в STM.
///    - Функция находит начало этой периодичности, анализируя первый период
///      и определяя, какой полупериод идет первым (min или max напряжения).

/// 2. РАЗДЕЛЕНИЕ НА ПРЯМОЕ (FW) И ОБРАТНОЕ (BW) СКАНИРОВАНИЕ
///    - FW: полупериод с растущим (или падающим) напряжением - содержит
///      информацию о VANC при нарастании напряжения
///    - BW: другой полупериод - информация при спадании напряжения
///    - Обратные кривые переворачиваются (reverse) для прямого сравнения с прямыми

/// 3. КОРРЕКЦИЯ ФАЗОВЫХ СДВИГОВ
///    - Туннельный ток и особенно шум могут иметь задержку относительно напряжения
///    - Применяются phase_shift_cur и phase_shift_n для синхронизации каналов

/// 4. ДЕЦИМАЦИЯ
///    - Каждый полупериод децимируется до p_num точек с усреднением
///    - Снижается шум и объем данных

/// 5. ФИЛЬТРАЦИЯ ПО КАЧЕСТВУ
///    - Вычисляется средняя (эталонная) VANC для FW и BW отдельно
///    - VANC, слишком отличающиеся от средней, удаляются
///    - Эталонная VANC пересчитывается на отфильтрованных данных

/// Выходные данные: пара усредненных VANC (FW и BW)
inline pair<vector<vector<float>>, vector<vector<float>>> process_vancs(
    vector<vector<float>>& data,    // Входные данные (3 канала)
    float crit_dev,                 // Критерий качества VANC
    int phase_shift_cur,            // Фазовый сдвиг канала тока
    int phase_shift_n,              // Фазовый сдвиг канала шума
    int period,                     // Полный период (2 полупериода)
    int p_num,                      // Выходное разрешение VANC
    int VAC_num,                    // Номер VANC для отчета
    Timer& tmr = *new Timer) {                   // Таймер для профилирования
    
    // === ЭТАП 1: РАЗБИЕНИЕ НА ПОЛУПЕРИОДЫ И НАХОЖДЕНИЕ НАЧАЛА ===
    int VAC_cnt = static_cast<int>(floor((data[0].size()) / period));

    vector<vector<float>> tmp(3, vector<float>(period, 0));
    vector<vector<vector<float>>> tmp_vancs(VAC_cnt, vector<vector<float>>(3, vector<float>(period, 0)));
    for (int i = 0; i < VAC_cnt; i++) {
        for (int k = 0; k < period; k++) {
            tmp[0][k] = data[0][period * i + k];
            tmp[1][k] = data[1][period * i + k];            
        }
        tmp_vancs[i] = tmp;
    }
    vector<vector<float>>avrg_tmp = avrg_curve(tmp_vancs);
    //reject_bad_vancs(tmp_vancs, avrg_tmp, crit_dev);
    //avrg_tmp = avrg_curve(tmp_vancs);
    //reject_bad_vancs(tmp_vancs, avrg_tmp, crit_dev);
    //avrg_tmp = avrg_curve(tmp_vancs);
    
    // Определяем фазу: находим минимум и максимум в каждом периоде
    // В идеальном треугольном сигнале развертки один из них идет первым
    int min_1pos = local_min(avrg_tmp[0], 0, 9 * period / 10);
    if (min_1pos == 0) min_1pos = period;
    int max_1pos = local_max(avrg_tmp[0], 0, 9 * period / 10);
    if (max_1pos == 0) max_1pos = period;
    //cout << min_1pos << endl << max_1pos << endl;
    int start = min(min_1pos, max_1pos);
    bool fw_first = min_1pos < max_1pos;  // true если минимум идет перед максимумом
    printf(" Start: %i      Start val: %f \n", start, data[0][start]);
    // === ЭТАП 2: РАЗДЕЛЕНИЕ НА ПОЛУПЕРИОДЫ И ФИЛЬТРАЦИЯ ===

    VAC_cnt = static_cast<int> (floor((data[0].size() - start) / period));
    vector<vector<float>> bw_tmp(3, vector<float>(period / 2, 0));
    vector<vector<float>> fw_tmp(3, vector<float>(period / 2, 0));
    vector<vector<vector<float>>> bw_vancs(VAC_cnt, vector<vector<float>>(3, vector<float>(p_num, 0)));
    vector<vector<vector<float>>> fw_vancs(VAC_cnt, vector<vector<float>>(3, vector<float>(p_num, 0)));
    
    // Определяем смещения для FW и BW в зависимости от фазы
    int offset_fw = (fw_first ? 0 : period / 2) + start;
    int offset_bw = (fw_first ? period / 2 : 0) + start;

    // Собираем FW и BW половины каждого периода с фазовыми сдвигами
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
        
        // Переворачиваем обратную VANC для совмещения с прямой
        // (обратная идет от max к min, нам нужно от min к max)
        std::reverse(bw_tmp[0].begin(), bw_tmp[0].end());
        std::reverse(bw_tmp[1].begin(), bw_tmp[1].end());
        std::reverse(bw_tmp[2].begin(), bw_tmp[2].end());
        //cout << i << endl;
        // Децимируем обе половины
        bw_vancs[i] = decimate_vanc(bw_tmp, p_num);
        fw_vancs[i] = decimate_vanc(fw_tmp, p_num);  
    }

    // === ЭТАП 3: ВЫЧИСЛЕНИЕ СРЕДНИХ И ФИЛЬТРАЦИЯ ===
    auto avrg_fw = avrg_curve(fw_vancs);
    auto avrg_bw = avrg_curve(bw_vancs);

    //std::cout <<"Start fw_bw deviation: " << mean_sq_deviation(avrg_bw, avrg_fw) << endl;
    /*reject_bad_vancs(fw_vancs, avrg_fw, 2 * crit_dev);
    reject_bad_vancs(bw_vancs, avrg_bw, 2 * crit_dev);
    avrg_fw = avrg_curve(fw_vancs);
    avrg_bw = avrg_curve(bw_vancs);*/

    reject_bad_vancs(fw_vancs, avrg_fw, crit_dev);
    reject_bad_vancs(bw_vancs, avrg_bw, crit_dev);

    // Пересчитываем средние на отфильтрованных данных
    avrg_fw = avrg_curve(fw_vancs);
    avrg_bw = avrg_curve(bw_vancs);

    // === КОНТРОЛЬ КАЧЕСТВА: ПРОВЕРКА НА ПОЛНОЕ ОТКЛОНЕНИЕ ВСЕХ VANC ===
    if ((avrg_fw.size() == 0) || (avrg_bw.size() == 0)) {
        printf("%i VAC done by thread %i. Time per file %f s; !!!!! ALL VANCS REJECTED !!!!\n", 
               VAC_num, omp_get_thread_num(), tmr.get_loop_interval() / 1000000.0);
        return make_pair(vector<vector<float>> (3, vector<float>(period / 2, EMPTY_VANC)), 
                         vector<vector<float>>(3, vector<float>(period / 2, EMPTY_VANC)));
    }   
    
    // === ФИНАЛЬНЫЙ ОТЧЕТ ===
    std::cout << "Finish FW/BW deviation: " << mean_sq_deviation(avrg_bw, avrg_fw) << endl;
    printf("%i VAC done by thread %i. Time per file %f s   FW count: %i    BW count: %i\n", 
           VAC_num, omp_get_thread_num(), tmr.get_loop_interval() / 1000000.0,
           static_cast<int>(fw_vancs.size()), static_cast<int>(bw_vancs.size()));
    
    return make_pair(avrg_fw, avrg_bw);
}

/// СОХРАНЕНИЕ УСРЕДНЕННЫХ VANC В ФАЙЛЫ
/// Вычисляет финальную интерполированную среднюю VANC и сохраняет отдельно
/// для прямого (FW) и обратного (BW) сканирования.

/// Формат файла: три столбца на строку (V, A, N) для каждой точки
inline void print_avrg_VANCS(
    pair<vector<vector<vector<float>>>, vector<vector<vector<float>>>>& VANCS,
    int p_num,                  // Количество точек выходных VANC
    int tm,                     // Время на точку
    string path) {              // Путь для сохранения
    
    auto avrg_fw = avrg_interpl_curve(VANCS.first, p_num);
    auto avrg_bw = avrg_interpl_curve(VANCS.second, p_num);
    ofstream  fw_vanc, bw_vanc;
    fw_vanc.open(path + "/VANC_FW_P" + to_string(p_num) + "_S" + to_string(tm) + ".dat");
    bw_vanc.open(path + "/VANC_BW_P" + to_string(p_num) + "_S" + to_string(tm) + ".dat");
    for (int i = 0; i < avrg_fw[0].size(); i++) {
        for (int k = 0; k < 3; k++) {
            fw_vanc << avrg_fw[k][i] << "  ";
            bw_vanc << avrg_bw[k][i] << "  ";
        }
        fw_vanc << endl;
        bw_vanc << endl;
    }
    AddFileToSession(path + "/VANC_FW_P" + to_string(p_num) + "_S" + to_string(tm) + ".dat");
}

// ============================================================================
// ГЛАВНАЯ ФУНКЦИЯ ПРОГРАММЫ
// ============================================================================
int main()
{
    // === ИНИЦИАЛИЗАЦИЯ ===
    omp_set_num_threads(THREADS);   // Параллельная обработка с использованием 4 потоков
    setlocale(LC_ALL, "Russian");   // Поддержка русского языка в консоли
    SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS);
    SetThreadPriority(GetCurrentProcess(), THREAD_PRIORITY_TIME_CRITICAL);
    std::cout << GetPriorityClass(GetCurrentProcess()) << endl;
    std::cout << "!!!Чтобы выбрать дефолтные значения вводите 0!!!!" << endl;
    
    // === ПАРАМЕТРЫ ПО УМОЛЧАНИЮ ===
    int start = 5;              // Номер первой VANC для обработки
    int stop =  400;            // Номер последней VANC для обработки
    int phase_shift_cur = 11;   // Фазовый сдвиг канала тока 
    int phase_shift_n = 7;      // Фазовый сдвиг канала шума 
    int p_num =  200;           // Разрешение выходной VANC (количество точек)
    float crit_dev = 0.008f;     // Критерий отклонения VANC (порог качества)
    float freq = 111.111f;      // Частота развертки напряжений (Гц) - определяет период
    string path = ReadLastVANCDiretory();
    std::cout << path << endl;
    float input_ = 0;

    int cnt = 4000000;  // Размер буфера для чтения из файла (макс. 4М отсчетов)
    
    // === ИНТЕРАКТИВНЫЙ ВВОД ПАРАМЕТРОВ ===
    std::cout << "Использовать путь к папке по умолчанию? (да = 1, нет = 0)" << endl;

    cin >> input_;
    if (input_ == 0){
        std::cout << "Введите путь к папке (пример \"11.04.2023/18_35/\"):" << endl;
        std::getline(std::cin, path);
        cin >> path;
        path = "C:/Users/Tunnel Noise/Desktop/STM/scans/"+path;
    }
    std::cout << path << endl;

    std::cout << "Использовать стандартные настройки? (да = 1, нет = 0)" << endl;

    cin >> input_;
    if (input_ == 0)
    {
        std::cout << "Введите номер первой кривой (default = " << start << "):" << endl;
        cin >> input_;
        if (input_ != 0) start = static_cast<int>(input_);

        std::cout << "Введите номер последней кривой (default <= " << stop << "):" << endl;
        cin >> input_;
        if (input_ != 0) stop = static_cast<int>(input_);

        std::cout << "Введите сдвиг фазы напряжения, в точках (default = " << phase_shift_cur << "):" << endl;
        cin >> input_;
        if (input_ != 0) phase_shift_cur = static_cast<int>(input_);

        std::cout << "Введите сдвиг фазы шума, в точках (default = " << phase_shift_n << "):" << endl;
        cin >> input_;
        if (input_ != 0)  phase_shift_n = static_cast<int>(input_);

        std::cout << "Введите желаемое количество точек на  обработанной кривой (default = " << p_num << "):" << endl;
        cin >> input_;
        if (input_ != 0)  p_num = static_cast<int>(input_);

        std::cout << "Введите критическое среднеквадратичное отклонение в единицах (10нА)^2 (default = " << crit_dev << "):" << endl;
        cin >> input_;
        if (input_ != 0)  crit_dev = input_;

        std::cout << "Введите частоту баяса (default = " << freq << "):" << endl;
        cin >> input_;
        if (input_ != 0)  freq = input_;
    }

    // Вычисляем период в отсчетах из частоты (333333 это базовая частота дискретизации)
    int period = static_cast<int>(333333.0f/freq);
    std::cout << "PERIOD = " << period << endl;
    
    // === ПОДГОТОВКА К ОБРАБОТКЕ ===
    string filename = path + "/VANC_";
    string filetype = ".bin";
    vector<vector<vector<float>>> fw_vancs, bw_vancs;  // Накопители для всех VANC из всех файлов
    Timer total_tmr;
    total_tmr.set_to_zero();

    // Проверяем, какие файлы существуют в папке (ограничиваем диапазон обработки)
    for (int i = start; i < stop; i++)
    {
        if (!std::filesystem::exists(filename + to_string(i) + "V" + filetype))
        {
            stop = i - 1;
        }
    }

    // === ПАРАЛЛЕЛЬНАЯ ОБРАБОТКА ФАЙЛОВ ===
    // Каждый файл содержит три канала (V, A, N) с множественными VANC
    // Файлы обрабатываются параллельно 4 потоками OpenMP
#pragma omp parallel for 
    for (int i = start; i < stop; i++)
    {
        Timer tmr;
        tmr.set_to_zero();                     

        FILE* fileV;
        FILE* fileA;
        FILE* fileN;
        
        // Открываем три файла (отдельные каналы) для каждого номера VANC
        std::cout << filename + to_string(i)  + filetype << endl;
        fopen_s(&fileV, (filename + to_string(i) + "V" + filetype).data(), "rb");
        fopen_s(&fileA, (filename + to_string(i) + "A" + filetype).data(), "rb");
        fopen_s(&fileN, (filename + to_string(i) + "N" + filetype).data(), "rb");
    
        // Читаем все три канала в буферы
        std::vector<float> Vbuf(cnt);
        std::vector<float> Abuf(cnt);
        std::vector<float> Nbuf(cnt);
        std::size_t size = std::fread(Vbuf.data(), sizeof Vbuf[0], Vbuf.size(), fileV);
        std::fread(Abuf.data(), sizeof Abuf[0], Abuf.size(), fileA);
        std::fread(Nbuf.data(), sizeof Nbuf[0], Nbuf.size(), fileN);
        
        // Объединяем три канала в единую структуру
        vector<vector<float> > data(3, vector<float>(size, 0));
        printf("!!!!!!!!!!!!!!!!! %i thread.  VANC num %i  !!!!!!!!!!!!!!!!\n file length: %zu \n",
               omp_get_thread_num(), i, size);     
        for (int j = 0; j < size; j++) {
            data[0][j] = Vbuf[j];
            data[1][j] = Abuf[j];
            data[2][j] = Nbuf[j];
        }
        fclose(fileV);
        fclose(fileA);
        fclose(fileN);

        // === ОСНОВНАЯ ОБРАБОТКА: РАЗДЕЛЕНИЕ НА FW/BW, ФИЛЬТРАЦИЯ И УСРЕДНЕНИЕ ===
        auto VANCS = process_vancs(data, crit_dev, phase_shift_cur, phase_shift_n, period, 2 * p_num, i, tmr);
        
        // Сохраняем результаты (если они не помечены как EMPTY_VANC)
        if (VANCS.first[0][0] != EMPTY_VANC) {
            fw_vancs.insert(fw_vancs.end(), VANCS.first);
            bw_vancs.insert(bw_vancs.end(), VANCS.second);
        }
    }

    // === ФИНАЛЬНОЕ УСРЕДНЕНИЕ ВСЕХ ОБРАБОТАННЫХ VANC ===
    std::cout << "FORWARD:" << endl;
    auto avrg_fw = avrg_interpl_curve(fw_vancs, 2 * p_num);
    std::cout << "BACKWARD:" << endl;
    auto avrg_bw = avrg_interpl_curve(bw_vancs, 2 * p_num);
    //auto avrg_fw = avrg_curve(fw_vancs);
    /*reject_bad_vancs(fw_vancs, avrg_fw, 2 * crit_dev);
    avrg_fw = avrg_curve(fw_vancs);*/
    reject_bad_vancs(fw_vancs, avrg_fw, crit_dev/25);
    /*reject_bad_vancs(bw_vancs, avrg_bw, 2 * crit_dev);
    avrg_bw = avrg_curve(bw_vancs);*/
    reject_bad_vancs(bw_vancs, avrg_bw, crit_dev/25);

    // === СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ===
    auto VANCS = make_pair(fw_vancs, bw_vancs);
    print_avrg_VANCS(VANCS, p_num, 2 * (stop - start), path);
    
    std::cout << "Done! Total time, s: " << total_tmr.get_full_interval() / 1000000;
    getchar(); getchar();
    
    return 0;
}
