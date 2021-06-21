#ifndef ALLINCLUDES_H
#define ALLINCLUDES_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include<fstream>
#include<iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <omp.h>
#include <limits.h>
#include <thread>

#ifdef _WIN64
    #include <windows.h>
    #include <excpt.h>
    #include <direct.h>
#endif

#ifdef linux
    #include <sys/types.h>
    #include <sys/stat.h>
    #include <unistd.h>
#endif

#include <chrono>
#include <random>
#include <numeric>
#include <functional>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#endif