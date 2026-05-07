#pragma once
#include <iostream>
#include <chrono>
#include <thread>
using namespace std;

/// ============================================================================
/// ИНСТРУМЕНТ ДЛЯ ТОЧНОГО ОТСЧЕТА МИКРОСЕКУНДНЫХ ИТЕРВАЛОВ
/// 
/// Используется для:
/// - Отсчета времени выполнения обработки одного файла VANC
/// - Профилирования производительности параллельной обработки
/// - Контроля общего времени сессии измерений
/// ============================================================================

/// Функция ожидания со спин-локом
/// Активно ожидает заданное время в микросекундах, используя высокоточный таймер.
/// Необходимо для синхронизации при управлении инструментом STM между
/// проведением отдельных VANC.
///
/// Параметр: usec - время ожидания в микросекундах
/// Возвращает: фактическое потраченное время в микросекундах
inline int uwait(double usec){
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    while (std::chrono::duration<double, std::micro>(end - start).count() < usec) {
        end = std::chrono::high_resolution_clock::now();
    }
    return (end - start).count();
}

/// CLASS TIMER - ВЫСОКОТОЧНЫЙ ТАЙМЕР ДЛЯ ПРОФИЛИРОВАНИЯ
/// 
/// Предназначен для отсчета времени выполнения критических операций
/// в программе обработки VANC.
///
/// Особенности:
/// - Использует высокоразрешающий таймер (chrono::high_resolution_clock)
/// - Разрешение микросекунда (1e-6 сек)
/// - Поддерживает отсчет "времени цикла" (loop_interval) и полного времени (full_interval)
/// - Позволяет сравнивать производительность при обработке разных файлов VANC
class Timer {
public:
    chrono::steady_clock::time_point start;  // Начало отсчета (инициализируется в конструкторе)
    chrono::steady_clock::time_point loop;   // Последняя точка для отсчета интервала
    chrono::steady_clock::time_point end;    // Текущий момент времени
    chrono::duration<double, std::micro> interval;  // Интервал в микросекундах
    
    /// Конструктор: инициализирует все метки текущим временем
    /// Интервал устанавливается в 0
    Timer() {
        start = std::chrono::high_resolution_clock::now();
        loop = std::chrono::high_resolution_clock::now();
        end = std::chrono::high_resolution_clock::now();
        interval = (loop - start);
    }
    
    ~Timer() {}
    
    /// Обнуляет таймер: сбрасывает все метки в текущий момент времени
    /// Используется для начала отсчета новой операции
    void set_to_zero() {
        start = loop = end = std::chrono::high_resolution_clock::now();
        interval = (loop - start);
    }
    
    /// Получить время, прошедшее с последнего вызова get_loop_interval()
    /// (или с set_to_zero(), если это первый вызов)
    ///
    /// Возвращает: время интервала в микросекундах (double)
    double get_loop_interval() {
        end = std::chrono::high_resolution_clock::now();
        interval = (end - loop);
        loop = end;
        return interval.count();
    }
    
    /// Получить полное время, прошедшее с момента инициализации таймера
    /// (или с момента последнего вызова set_to_zero())
    ///
    /// Используется для подсчета общего времени обработки сессии всех VANC.
    ///
    /// Возвращает: полное время в микросекундах (double)
    double get_full_interval() {
        end = std::chrono::high_resolution_clock::now();
        interval = (end - start);
        return interval.count();
    }
}; 
