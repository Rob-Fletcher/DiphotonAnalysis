#ifndef GravitonGGAnalysis_DEBUG_H
#define GravitonGGAnalysis_DEBUG_H

//#define DEBUG_BUILD

#ifdef DEBUG_BUILD
#include <iostream>
#include <ostream>
#include <algorithm>
#endif

#ifdef DEBUG_BUILD
#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)
#else
#define DEBUG(x) do { } while (0)
#endif

#ifdef DEBUG_BUILD
#define DEBUG_ITERABLE(type, x) do { std::copy((x).begin(), (x).end(), std::ostream_iterator<type>(std::cerr, " ")); } while(0)
#else
#define DEBUG_ITERABLE(type, x) do { } while (0)
#endif

#endif //GravitonGGAnalysis_DEBUG_H
