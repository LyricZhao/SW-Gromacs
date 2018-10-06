# include "slave.h"

__thread_kernel("compute") int *answer;

void cpe_call(long *paras) {
  answer = (int *) paras[0];
  int threadID = athread_get_id(-1);
  answer[threadID] = threadID;
  return;
}
