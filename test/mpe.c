# include "athread.h"

extern SLAVE_FUN(cpe_call)();

int main() {
  int answer[64];
  for(int i = 0; i < 64; ++ i) answer[i] = 0;
  athread_init();
  athread_enter64();

  long paras[1] = {(long) answer};
  athread_spawn64(SLAVE_FUN(cpe_call), &paras);
  for(int i = 0; i < 64; ++ i) printf("%d\n", answer[i]);

  athread_leave64();
  athread_halt();
  return 0;
}
