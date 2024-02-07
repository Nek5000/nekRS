#if !defined(_CRS_BOX_TIMER_HPP_)
#define _CRS_BOX_TIMER_HPP_

enum BOX_METRIC {
  COPY_RHS = 0,
  CRS_DSAVG1 = 1,
  ASM1 = 2,
  CRS_DSAVG2 = 3,
  MULT_RHS_UPDATE = 4,
  COPY_TO_NEK5000 = 5,
  MAP_VTX_TO_BOX = 6,
  ASM2 = 7,
  MAP_BOX_TO_VTX = 8,
  COPY_FROM_NEK5000 = 9,
  CRS_DSAVG3 = 10,
  COPY_SOLUTION = 11,
  NONE = 100
};

void timer_init();
void timer_tic(const struct comm *c);
void timer_toc(BOX_METRIC m);
void timer_dump(struct comm *c, unsigned interval);
void timer_print(struct comm *c, unsigned interval);

#endif
