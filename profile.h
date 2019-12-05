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
        struct timespec name##_start, name##_end; \
        double name##_cost_sum = 0.0; \
        double *name##_cost = 0;

#define PROFILE_START(name) \
        extern struct timespec name##_start, name##_end; \
        extern double name##_cost_sum; \
        extern double *name##_cost; \
        if(! name##_cost) \
            name##_cost = (double*)calloc(PROFILE_THREAD_NUM, sizeof(double)); \
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &name##_start);

#define PROFILE_END(name) \
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &name##_end); \
        name##_cost[pro_tid] += name##_end.tv_sec - name##_start.tv_sec + (name##_end.tv_nsec - name##_start.tv_nsec) / 1e9;

#define PROFILE_REPORT(name) \
        for(int i=0; i<PROFILE_THREAD_NUM; i++){ \
            name##_cost_sum += name##_cost[i]; \
        } \
        fprintf(stderr, "[WangTong Profile] CPU: %s %fsecs \n", #name, name##_cost_sum);

extern __thread int pro_tid;
#endif //MINIMAP2_PROFILE_H
