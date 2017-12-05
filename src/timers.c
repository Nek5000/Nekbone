
#include <time.h>  /* for struct timespec */

double fclock_gettime(void) {

  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  double timeval = (double) (ts.tv_sec + ts.tv_nsec*1.0E-9);

  return timeval;
}

double fclock_gettime_( void ) __attribute__((alias("fclock_gettime")));

#ifdef BGQTIMER
#include <hwi/include/bqc/A2_inlines.h>

double ReadTimeBase_Double( void ) {
  return (double) GetTimeBase();
}

double readtimebase_double( void ) __attribute__((alias("ReadTimeBase_Double")));
double readtimebase_double_( void ) __attribute__((alias("ReadTimeBase_Double")));
#endif
