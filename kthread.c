#include <pthread.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include "kthread.h"

#if (defined(WIN32) || defined(_WIN32)) && defined(_MSC_VER)
#define __sync_fetch_and_add(ptr, addend)     _InterlockedExchangeAdd((void*)ptr, addend)
#endif

/************
 * kt_for() *
 ************/

struct kt_for_t;

typedef struct {
	struct kt_for_t *t;
	long i;
} ktf_worker_t;

typedef struct kt_for_t {
	int n_threads;
	long n;
	ktf_worker_t *w;
	void (*func)(void*,long,int);
	void *data;
} kt_for_t;

static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		if (i >= w->t->n) break;
		w->t->func(w->t->data, i, w - w->t->w);
	}
	while ((i = steal_work(w->t)) >= 0)
		w->t->func(w->t->data, i, w - w->t->w);
	pthread_exit(0);
}

/**
 * 多线程运行，func函数n次，data和i作为func函数的参数，其中i表示第几次运行func函数。
 *
 * @param n_threads
 * @param func
 * @param data
 * @param n
 */
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n)
{
	if (n_threads > 1) {
		int i;
		kt_for_t t;
		pthread_t *tid;
		t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
		t.w = (ktf_worker_t*)calloc(n_threads, sizeof(ktf_worker_t));
		tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
		for (i = 0; i < n_threads; ++i)
			t.w[i].t = &t, t.w[i].i = i;
		for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);
		for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
		free(tid); free(t.w);
	} else {
		long j;
		for (j = 0; j < n; ++j) func(data, j, 0);
	}
}

/*****************
 * kt_pipeline() *
 *****************/

struct ktp_t;

typedef struct {
	struct ktp_t *pl;
	int64_t index;
	int step;
	void *data;
} ktp_worker_t;

typedef struct ktp_t {
	void *shared;
	void *(*func)(void*, int, void*);
	int64_t index;
	int n_workers, n_steps;
	ktp_worker_t *workers;
	pthread_mutex_t mutex;
	pthread_cond_t cv;
} ktp_t;

/**
 * 该函数循环执行ktp_worker_t *data->func函数，func函数的参数中有int参数，表示func函数的所执行的step。
 * 每次执行该函数时，step的值循环增长。并且多个该函数并行执行时，func函数的执行是有限制的，因为不同线程的
 * step之间需要同步，保证不会出现index较小的线程在执行和 小于等于其他index较大的线程的step。
 *
 * @param data		被解析为 ktp_worker_t 结构体指针
 * */
static void *ktp_worker(void *data)
{
	ktp_worker_t *w = (ktp_worker_t*)data;
	ktp_t *p = w->pl;
	while (w->step < p->n_steps) {
		// test whether we can kick off the job with this worker
		pthread_mutex_lock(&p->mutex);
		for (;;) {
			int i;
			// test whether another worker is doing the same step
			for (i = 0; i < p->n_workers; ++i) {
				if (w == &p->workers[i]) continue; // ignore itself
				if (p->workers[i].step <= w->step && p->workers[i].index < w->index)
					break;
			}
			if (i == p->n_workers) break; // no workers with smaller indices are doing w->step or the previous steps
			pthread_cond_wait(&p->cv, &p->mutex);
		}
		pthread_mutex_unlock(&p->mutex);

		// working on w->step
		w->data = p->func(p->shared, w->step, w->step? w->data : 0); // for the first step, input is NULL // 第一步执行前, w->data就是NULL；执行完后w->data才有内容，data一般是func函数中需要保存的中间信息（step_t结构体）。

		// update step and let other workers know
		pthread_mutex_lock(&p->mutex);
		w->step = w->step == p->n_steps - 1 || w->data? (w->step + 1) % p->n_steps : p->n_steps; // 到了最后一步，或者有数据 都会正常更新step；如果没到最后一步就没了数据，就step更新到n_steps。
		if (w->step == 0) w->index = p->index++; // 更新完step之后，如果某线程的step为0，那么就要更新该线程的index。
		pthread_cond_broadcast(&p->cv);
		pthread_mutex_unlock(&p->mutex);
	}
	pthread_exit(0);
}

/**
 * 开n个线程，每个线程执行ktp_worker函数，ktp_worker函数则循环执行func函数。
 * func函数被执行的时候会有一定的同步。
 *
 * @param n_threads 	线程数
 * @param func			存放于aux->func中，在ktp_worker函数中会执行，ktp_worker的参数是ktp_worker_t结构体，该结构体的pl字段指向aux。
 * @param shared_data 	存放于aux->shared
 * @param n_steps		存放于aux->n_steps
 * */
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps)
{
	ktp_t aux;
	pthread_t *tid;
	int i;

	if (n_threads < 1) n_threads = 1;
	aux.n_workers = n_threads;
	aux.n_steps = n_steps;
	aux.func = func;
	aux.shared = shared_data;
	aux.index = 0;
	pthread_mutex_init(&aux.mutex, 0);
	pthread_cond_init(&aux.cv, 0);

	aux.workers = (ktp_worker_t*)calloc(n_threads, sizeof(ktp_worker_t));
	for (i = 0; i < n_threads; ++i) {
		ktp_worker_t *w = &aux.workers[i];
		w->step = 0; w->pl = &aux; w->data = 0;
		w->index = aux.index++;
	}

	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktp_worker, &aux.workers[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	free(tid); free(aux.workers);

	pthread_mutex_destroy(&aux.mutex);
	pthread_cond_destroy(&aux.cv);
}
