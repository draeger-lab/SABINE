#include "mytimer.h"

using namespace std;

CTimer::CTimer()
{
  reinit();
}

void CTimer::start()
{
  t_begin = clock();
}

void CTimer::stop()
{
  clock_t t_end = clock();

  //check for overflow - little inexact if not only running process
  t_run += (t_end-t_begin)/(double)CLOCKS_PER_SEC;
}

double CTimer::get_time()
{
  stop();
  return(t_run);
}

void CTimer::reinit()
{
  t_run = 0;
  start();
}

void wait ( int seconds )
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}
