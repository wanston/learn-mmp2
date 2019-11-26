//
// Created by tong on 19-6-5.
//

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "ksw2.h"

#ifdef __SSE2__
#include <emmintrin.h>

#ifdef KSW_SSE2_ONLY
#undef __SSE4_1__
#endif

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

#ifdef KSW_CPU_DISPATCH
#ifdef __SSE4_1__
void ksw_extd2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#else
void ksw_extd2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#endif
#else

void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#endif // ~KSW_CPU_DISPATCH
{
#define swap(A, B) ({typeof(A) __tmp__ = A; A = B; B = __tmp__})
#define max(A, B) ({typeof(A) __tmpA__ = A, __tmpB__ = B;  __tmpA__ > __tmpB__ ? __tmpA__ : __tmpB__;})
#define min(A, B) ({typeof(A) __tmpA__ = A, __tmpB__ = B;  __tmpA__ < __tmpB__ ? __tmpA__ : __tmpB__;})

    /*************************/
    printf("qlen %d tlen %d\n", qlen, tlen);
    /*************************/

    ksw_reset_extz(ez); // 初始化用于记录返回值的结构体
    if (m <= 1 || qlen <= 0 || tlen <= 0) return;

    if (q2 + e2 < q + e){
        int8_t t;
        t = q, q = q2, q2 = t;
        t = e, e = e2, e2 = t;
    }

	if(w < 0)
        w = tlen > qlen? tlen : qlen;
	int wl = w, wr = w; // wl和wr被用来确定rt坐标系中t的范围

    int n_col = qlen < tlen ? qlen : tlen;
    n_col = (n_col < w + 1 ? n_col : w + 1);

    int max_sc = mat[0], min_sc = mat[1];
    for (int i = 1; i < m * m; ++i) {
        max_sc = max_sc > mat[i] ? max_sc : mat[i];
        min_sc = min_sc < mat[i] ? min_sc : mat[i];
    }
    if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

    int sz = tlen + 2;
    int32_t *mem = (int32_t*)kcalloc(km, (size_t)sz*7* sizeof(int) + sz + tlen + qlen, 1); // tlen_ 和qlen_ 是向上对16取整，然后除以16得到的结果

    int32_t *E = mem;
    int32_t *E2 = E + sz;
    int32_t *F = E2 + sz;
    int32_t *F2 = F + sz;
    int32_t *H = F2 + sz;
    int32_t *H1 = H + sz;
    int32_t *H2 = H1 + sz;
    int8_t *s = (int8_t*)(H2 + sz);
    int8_t *sf = s + sz;
    int8_t *qr = sf + tlen;

    int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);

    int *off = NULL, *off_end = NULL;
    uint8_t *p = NULL;
    if (with_cigar) {
        p = (uint8_t*)kmalloc(km, (size_t)(qlen + tlen - 1) * n_col);
        off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
        off_end = off + qlen + tlen - 1;
    }

    for (int i = 0; i < qlen; ++i) qr[i] = query[qlen - 1 - i];
    memcpy(sf, target, (size_t)tlen);

    int8_t sc_mch = mat[0], sc_mis = mat[1], sc_N = mat[m*m-1] == 0 ? -e2 : mat[m*m-1];

    int H_init = max(-q-e, -q2-e2);
    int E_init = -2*q-e, F_init = -2*q-e, E2_init = -2*q2-e2, F2_init = -2*q2-e2;

    H2[0] = 0; // TODO
    H1[0] = H1[1] = H_init;
    E[0] = E_init, E2[0] = E2_init;
    F[1] = F_init, F2[1] = F2_init;

    int H0 = 0, last_H0_t = 0;

    for(int r=0; r < qlen+tlen-1; ++r){
        int st = 0, en = tlen-1;
        int8_t *qrr = qr + (qlen - 1 - r);

        if (st < r - qlen + 1) st = r - qlen + 1;
        if (en > r) en = r;
        if (st < (r-wr+1)>>1) st = (r-wr+1)>>1;
        if (en > (r+wl)>>1) en = (r+wl)>>1;
        if (st > en) {
            ez->zdropped = 1;
            break;
        }

        if(!(flag & KSW_EZ_GENERIC_SC)){
            for(int i=st; i<=en; ++i){
                s[i+1] = sf[i] == qrr[i] ? sc_mch : sc_mis;
                if(qrr[i] == m-1 || sf[i] == m-1) s[i+1] = sc_N;
            }
        }else{
            for(int i=st; i<=en; ++i){
                s[i+1] = mat[sf[i] + qrr[i]];
            }
        }

        if(!with_cigar){
            for(int t=st+1; t<=en+1; ++t){
                int M = H2[t-1] + s[t];
                int E_rt = max(H1[t-1]-q-e, E[t-1]-e);
                int E2_rt = max(H1[t-1]-q2-e2, E2[t-1]-e2);
                int F_rt = max(H1[t]-q-e, F[t]-e);
                int F2_rt = max(H1[t]-q2-e2, F2[t]-e2);

                H[t] = max(M, max(E_rt, max(E2_rt, max(F_rt, F2_rt))));
                E[t] = E_rt;
                E2[t] = E2_rt;
                F[t] = F_rt;
                F2[t] = F2_rt;
            }
        }else if (!(flag&KSW_EZ_RIGHT)){
            uint8_t *pr = p + r*n_col - st;
            off[r] = st, off_end[r] = en;

            for(int t=en+1; t>=st+1; --t){
                int M = H2[t-1] + s[t];
                int E_rt = max(H1[t-1]-q-e, E[t-1]-e);
                int E2_rt = max(H1[t-1]-q2-e2, E2[t-1]-e2);
                int F_rt = max(H1[t]-q-e, F[t]-e);
                int F2_rt = max(H1[t]-q2-e2, F2[t]-e2);

                uint8_t d = 0;
                H[t] = E_rt > M ? (d=1, E_rt) : (d=0, M);
                H[t] = F_rt > H[t] ? (d=2, F_rt) : H[t];
                H[t] = E2_rt > H[t] ? (d=3, E2_rt) : H[t];
                H[t] = F2_rt > H[t] ? (d=4, F2_rt) : H[t];

//                d = (uint8_t)(E[t-1]  > H1[t-1]-q  ? 1 << 3 : 0) | d;
//                d = (uint8_t)(F[t]    > H1[t]-q    ? 1 << 4 : 0) | d;
//                d = (uint8_t)(E2[t-1] > H1[t-1]-q2 ? 1 << 5 : 0) | d;
//                d = (uint8_t)(F2[t]   > H1[t]-q2   ? 1 << 6 : 0) | d;

                d = (uint8_t)(E_rt  > H[t]-q  ? 1 << 3 : 0) | d;
                d = (uint8_t)(F_rt  > H[t]-q  ? 1 << 4 : 0) | d;
                d = (uint8_t)(E2_rt > H[t]-q2 ? 1 << 5 : 0) | d;
                d = (uint8_t)(F2_rt > H[t]-q2 ? 1 << 6 : 0) | d;

                pr[t-1] = d;

                E[t] = E_rt;
                E2[t] = E2_rt;
                F[t] = F_rt;
                F2[t] = F2_rt;
            }
            /*************************/
            for(int t=st; t<=en; ++t){
                printf("0x%x\t", pr[t]);
            }
            printf("\n");
            /*************************/
        }else{
            uint8_t *pr = p + r*n_col - st;
            off[r] = st, off_end[r] = en;

            for(int t=en+1; t>=st+1; --t){
                int M = H2[t-1] + s[t];
                int E_rt = max(H1[t-1]-q-e, E[t-1]-e);
                int E2_rt = max(H1[t-1]-q2-e2, E2[t-1]-e2);
                int F_rt = max(H1[t]-q-e, F[t]-e);
                int F2_rt = max(H1[t]-q2-e2, F2[t]-e2);

                uint8_t d = 0;
                H[t] = E_rt >= M ? (d=1, E_rt) : (d=0, M);
                H[t] = F_rt >= H[t] ? (d=2, F_rt) : H[t];
                H[t] = E2_rt >= H[t] ? (d=3, E2_rt) : H[t];
                H[t] = F2_rt >= H[t] ? (d=4, F2_rt) : H[t];

                d = (uint8_t)(E_rt  >= H[t]-q  ? 1 << 3 : 0) | d;
                d = (uint8_t)(F_rt  >= H[t]-q  ? 1 << 4 : 0) | d;
                d = (uint8_t)(E2_rt >= H[t]-q2 ? 1 << 5 : 0) | d;
                d = (uint8_t)(F2_rt >= H[t]-q2 ? 1 << 6 : 0) | d;

                pr[t-1] = d;

                E[t] = E_rt;
                E2[t] = E2_rt;
                F[t] = F_rt;
                F2[t] = F2_rt;
            }
            /*************************/
            for(int t=st; t<=en; ++t){
                printf("0x%x\t", pr[t]);
            }
            printf("\n");
            /*************************/
        }
        H[st] = H[en+2] = -min(q+e*(r+2), q2+e2*(r+2));
//        H[st] = H[en+2] = r+1 <= long_thres ? (H_init -= e) : (H_init -= e2); // TODO: BUG所在
        int *pt = H; H = H2; H2 = H1; H1 = pt;

        E[st] = (E_init-=e);
        E2[st] = (E2_init-=e2);
        F[en+2] = (F_init-=e);
        F2[en+2] = (F2_init-=e2);

        /*************************/
        if(with_cigar){
//            printf("st %d en %d\n", st, en);
//            for(int t=0; t<st; ++t){
//                printf("*\t\t\t");
//            }
//            for(int t=st+1; t<=en+1+1; ++t){
//                printf("(%d %d %d %d %d)\t", H1[t], E[t], F[t], E2[t], F2[t]);
//            }
            for(int t=st+1; t<=en+1; ++t){
                printf("%d\t", H1[t]);
            }
            printf("\n");
        }
        /*************************/

        if(!approx_max){
            int32_t max_H = H1[1+st], max_t = st;
            for(int t=st; t<=en; ++t){
                if(H1[1+t] > max_H){
                    max_H = H1[1+t];
                    max_t = t;
                }
            }

            if(en == tlen - 1 && H1[1+en] > ez->mte)
                ez->mte = H1[1+en], ez->mte_q = r - en; // todo: 原来是16对齐的en

            if (r - st == qlen - 1 && H1[1+st] > ez->mqe)
                ez->mqe = H1[1+st], ez->mqe_t = st;
            if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) break;
            if (r == qlen + tlen - 2 && en == tlen - 1)
                ez->score = H1[1 + tlen - 1];
        }else{
            if(r > 0){
                if(last_H0_t >= st && last_H0_t <= en && last_H0_t+1 >= st && last_H0_t+1 <= en){
                    int d0 = H1[1+last_H0_t] - H2[1+last_H0_t];
                    int d1 = H1[1+last_H0_t+1] - H2[1+last_H0_t];
                    if(d0 > d1) H0 += d0;
                    else H0 += d1, ++last_H0_t;
                }else if(last_H0_t >= st && last_H0_t <= en){
                    H0 += H1[1+last_H0_t] - H2[1+last_H0_t];
                }else{
                    ++last_H0_t, H0 += H1[1+last_H0_t] - H2[1+last_H0_t-1];
                }
            }else{
                H0 = H1[1] - H2[1] - q - e;
                last_H0_t = 0;
            }

            if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) break;
            if (r == qlen + tlen - 2 && en == tlen - 1)
                ez->score = H0;
        }
    }

    kfree(km, mem);

    if(with_cigar){
        int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
        if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
            ksw_backtrack(km, 1, rev_cigar, 0, p, off, off_end, n_col, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        } else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
            ez->reach_end = 1;
            ksw_backtrack(km, 1, rev_cigar, 0, p, off, off_end, n_col, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        } else if (ez->max_t >= 0 && ez->max_q >= 0) {
            ksw_backtrack(km, 1, rev_cigar, 0, p, off, off_end, n_col, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        }
        kfree(km, p);
        kfree(km, off);
    }
    /*************************/
    printf("\n");
    print_cigar(ez);
    printf("\n");
    /*************************/
}
#endif // __SSE2__ # 磁盘IO
