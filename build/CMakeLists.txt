cmake_minimum_required(VERSION 2.8)

project(minimap2)


set(SRC ../align.c
        ../bseq.c
        ../chain.c
        ../esterr.c
        ../format.c
        ../hit.c
        ../index.c
        ../kalloc.c
        ../kalloc.h
        ../ksw2_dispatch.c
        ../ksw2_extd2_sse.c
        ../ksw2_exts2_sse.c
        ../ksw2_extz2_sse.c
        ../ksw2_ll_sse.c
        ../kthread.c
        ../kvec.h
        ../main.c
        ../map.c
        ../misc
        ../misc.c
        ../options.c
        ../pe.c
        ../sdust.c
        ../sketch.c
        ../splitidx.c)

add_executable(minimap2 ${SRC})
target_link_libraries(minimap2 pthread m z)