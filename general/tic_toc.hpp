#include <sys/times.h>

class StopWatch
{
   private:
      clock_t real_time, user_time, syst_time;
      clock_t start_rtime, start_utime, start_stime;
      long my_CLK_TCK;
      short Running;
      void Current(clock_t *, clock_t *, clock_t *);
   public:
      StopWatch();  // determines my_CLK_TCK
      void Clear();
      void Start();
      void Stop();
      double RealTime();
      double UserTime();
      double SystTime();
};


extern StopWatch tic_toc;

extern void tic();

extern double toc();

