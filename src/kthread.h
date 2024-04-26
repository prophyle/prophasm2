#ifndef KTHREAD_H
#define KTHREAD_H

#ifdef __cplusplus
extern "C" {
#endif

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

#ifdef __cplusplus
}
#endif

#endif
