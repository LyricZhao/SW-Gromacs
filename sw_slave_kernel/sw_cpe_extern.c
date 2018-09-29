# include "athread.h"

void athread_INIT() {
  athread_init();
  athread_enter64();
}

void athread_LEAVE() {
  athread_leave64();
  athread_halt();
}
