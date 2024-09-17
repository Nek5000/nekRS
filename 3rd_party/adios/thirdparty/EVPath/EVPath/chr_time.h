#ifndef CHR_TIME_H
#define CHR_TIME_H

#include <stddef.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

  typedef struct chr_time {
    double d1;
    double d2;
    double d3;
  } chr_time;
  
  extern void
  chr_get_time( chr_time *time);
  
  extern void
  chr_timer_diff( chr_time *diff_time, chr_time *src1, chr_time *src2);
  
  extern int
  chr_timer_eq_zero( chr_time *time);
  
  extern void
  chr_timer_sum( chr_time *sum_time, chr_time *src1, chr_time *src2);
  
  extern void
  chr_timer_start (chr_time *timer);
  
  extern void
  chr_timer_stop (chr_time *timer);
  
  extern double
  chr_time_to_nanosecs (chr_time *time);
  
  extern double
  chr_time_to_microsecs (chr_time *time);
  
  extern double
  chr_time_to_millisecs (chr_time *time);
  
  extern double
  chr_time_to_secs (chr_time *time);
  
  extern double 
  chr_approx_resolution();
  
#if defined(__cplusplus) || defined(c_plusplus)
}
#endif
#endif
