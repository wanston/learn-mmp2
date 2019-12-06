//
// Created by tong on 19-12-5.
//

#ifndef MINIMAP2_PROFILE_H
#define MINIMAP2_PROFILE_H

#include <time.h>

#define PROFILE_THREAD_NUM 64

#define PROFILE_INIT0 \
        __thread int pro_tid;

#define PROFILE_SET_TID(i) (pro_tid = (i));

#define PROFILE_INIT(name) \
        struct timespec *name##_start, *name##_end; \
        double name##_cost_sum = 0.0; \
        double *name##_cost = 0;

#define PROFILE_START(name) \
        extern struct timespec *name##_start, *name##_end; \
        extern double *name##_cost; \
        if(! name##_cost) \
            name##_cost = (double*)calloc(PROFILE_THREAD_NUM, sizeof(double)); \
        if(! name##_start) \
            name##_start = (struct timespec*)malloc(PROFILE_THREAD_NUM * sizeof(struct timespec)); \
        if(! name##_end) \
            name##_end = (struct timespec*)malloc(PROFILE_THREAD_NUM * sizeof(struct timespec)); \
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, name##_start + pro_tid);

#define PROFILE_END(name) \
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, name##_end + pro_tid); \
        name##_cost[pro_tid] += name##_end[pro_tid].tv_sec - name##_start[pro_tid].tv_sec + (name##_end[pro_tid].tv_nsec - name##_start[pro_tid].tv_nsec) / 1e9;

#define PROFILE_REPORT(name) \
        extern double name##_cost_sum; \
        int name##_i; \
        for(name##_i=0; name##_i<PROFILE_THREAD_NUM; name##_i++){ \
            name##_cost_sum += name##_cost[name##_i]; \
        } \
        fprintf(stderr, "[WangTong Profile] CPU: %s %.3f seconds \n", #name, name##_cost_sum);

extern __thread int pro_tid;
#endif //MINIMAP2_PROFILE_H
