#ifndef MYTIMER_H
#define MYTIMER_H

#include <ctime>
#include <climits>

void wait ( int seconds );

class CTimer
{
 public:
  CTimer();
  double get_time();
  void start();
  void stop();
  void reinit();
 private:
  clock_t t_begin;
  double t_run;
  time_t tt_begin;
  time_t tt_run;
};

#endif
