// Copyright (c) 2019-2022 University of Oregon
// Distributed under the BSD Software License
// (See accompanying file LICENSE.txt)

#ifndef _GNU_SOURCE
#define _GNU_SOURCE // needed to define RTLD_DEFAULT
#endif
#include <dlfcn.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pthread.h"
#ifndef PERFSTUBS_USE_TIMERS
#define PERFSTUBS_USE_TIMERS
#endif
#include "perfstubs_api/timer.h"

/* Make sure that the Timer singleton is constructed when the
 * library is loaded.  This will ensure (on linux, anyway) that
 * we can assert that we have m_Initialized on the main thread. */
//static void __attribute__((constructor)) initialize_library(void);

/* Globals for the plugin API */

int perfstubs_initialized = PERFSTUBS_UNKNOWN;
int num_tools_registered = 0;
/* Keep track of whether the thread has been registered */
/* __thread int thread_seen = 0; */
/* Implemented with PGI-friendly implementation, they can't be bothered
 * to implement the thread_local standard like every other compiler... */
static pthread_key_t key;
static pthread_once_t key_once = PTHREAD_ONCE_INIT;

static void make_key(void) {
    (void) pthread_key_create(&key, NULL);
}

/* Function pointers */

ps_initialize_t initialize_function;
ps_finalize_t finalize_function;
ps_pause_measurement_t pause_measurement_function;
ps_resume_measurement_t resume_measurement_function;
ps_register_thread_t register_thread_function;
ps_dump_data_t dump_data_function;
ps_timer_create_t timer_create_function;
ps_timer_start_t timer_start_function;
ps_timer_stop_t timer_stop_function;
ps_start_string_t start_string_function;
ps_stop_string_t stop_string_function;
ps_stop_current_t stop_current_function;
ps_set_parameter_t set_parameter_function;
ps_dynamic_phase_start_t dynamic_phase_start_function;
ps_dynamic_phase_stop_t dynamic_phase_stop_function;
ps_create_counter_t create_counter_function;
ps_sample_counter_t sample_counter_function;
ps_set_metadata_t set_metadata_function;
ps_get_timer_data_t get_timer_data_function;
ps_get_counter_data_t get_counter_data_function;
ps_get_metadata_t get_metadata_function;
ps_free_timer_data_t free_timer_data_function;
ps_free_counter_data_t free_counter_data_function;
ps_free_metadata_t free_metadata_function;

#ifdef PERFSTUBS_USE_STATIC

#if defined(__clang__) && defined(__APPLE__)
#define PS_WEAK_PRE
#define PS_WEAK_POST __attribute__((weak_import))
#define PS_WEAK_POST_NULL __attribute__((weak_import))
#else
#define PS_WEAK_PRE __attribute__((weak))
#define PS_WEAK_POST
#define PS_WEAK_POST_NULL
#endif

PS_WEAK_PRE void ps_tool_initialize(void) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_finalize(void) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_pause_measurement(void) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_resume_measurement(void) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_register_thread(void) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_dump_data(void) PS_WEAK_POST;
PS_WEAK_PRE void* ps_tool_timer_create(const char *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_timer_start(void *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_timer_stop(void *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_start_string(const char *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_stop_string(const char *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_stop_current(void) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_set_parameter(const char *, int64_t) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_dynamic_phase_start(const char *, int) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_dynamic_phase_stop(const char *, int) PS_WEAK_POST;
PS_WEAK_PRE void* ps_tool_create_counter(const char *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_sample_counter(void *, double) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_set_metadata(const char *, const char *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_get_timer_data(ps_tool_timer_data_t *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_get_counter_data(ps_tool_counter_data_t *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_get_metadata(ps_tool_metadata_t *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_free_timer_data(ps_tool_timer_data_t *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_free_counter_data(ps_tool_counter_data_t *) PS_WEAK_POST;
PS_WEAK_PRE void ps_tool_free_metadata(ps_tool_metadata_t *) PS_WEAK_POST;
#endif

void initialize_library(void) {
#ifdef PERFSTUBS_USE_STATIC
    /* The initialization function is the only required one */
    initialize_function = &ps_tool_initialize;
    if (initialize_function == NULL) {
        perfstubs_initialized = PERFSTUBS_FAILURE;
        return;
    }
    // removing printf statement for now, it's too noisy.
    //printf("Found ps_tool_initialize(), registering tool\n");
    finalize_function = &ps_tool_finalize;
    pause_measurement_function = &ps_tool_pause_measurement;
    resume_measurement_function = &ps_tool_resume_measurement;
    register_thread_function = &ps_tool_register_thread;
    dump_data_function = &ps_tool_dump_data;
    timer_create_function = &ps_tool_timer_create;
    timer_start_function = &ps_tool_timer_start;
    timer_stop_function = &ps_tool_timer_stop;
    start_string_function = &ps_tool_start_string;
    stop_string_function = &ps_tool_stop_string;
    stop_current_function = &ps_tool_stop_current;
    set_parameter_function = &ps_tool_set_parameter;
    dynamic_phase_start_function = &ps_tool_dynamic_phase_start;
    dynamic_phase_stop_function = &ps_tool_dynamic_phase_stop;
    create_counter_function = &ps_tool_create_counter;
    sample_counter_function = &ps_tool_sample_counter;
    set_metadata_function = &ps_tool_set_metadata;
    get_timer_data_function = &ps_tool_get_timer_data;
    get_counter_data_function = &ps_tool_get_counter_data;
    get_metadata_function = &ps_tool_get_metadata;
    free_timer_data_function = &ps_tool_free_timer_data;
    free_counter_data_function = &ps_tool_free_counter_data;
    free_metadata_function = &ps_tool_free_metadata;
#else
    initialize_function =
        (ps_initialize_t)dlsym(RTLD_DEFAULT, "ps_tool_initialize");
    if (initialize_function == NULL) {
        perfstubs_initialized = PERFSTUBS_FAILURE;
        return;
    }
    printf("Found ps_tool_initialize(), registering tool\n");
    finalize_function =
        (ps_finalize_t)dlsym(RTLD_DEFAULT, "ps_tool_finalize");
    pause_measurement_function =
        (ps_pause_measurement_t)dlsym(RTLD_DEFAULT, "ps_tool_pause_measurement");
    resume_measurement_function =
        (ps_resume_measurement_t)dlsym(RTLD_DEFAULT, "ps_tool_resume_measurement");
    register_thread_function =
        (ps_register_thread_t)dlsym(RTLD_DEFAULT, "ps_tool_register_thread");
    dump_data_function =
        (ps_dump_data_t)dlsym(RTLD_DEFAULT, "ps_tool_dump_data");
    timer_create_function =
        (ps_timer_create_t)dlsym(RTLD_DEFAULT,
                "ps_tool_timer_create");
    timer_start_function =
        (ps_timer_start_t)dlsym(RTLD_DEFAULT, "ps_tool_timer_start");
    timer_stop_function =
        (ps_timer_stop_t)dlsym(RTLD_DEFAULT, "ps_tool_timer_stop");
    start_string_function =
        (ps_start_string_t)dlsym(RTLD_DEFAULT, "ps_tool_start_string");
    stop_string_function =
        (ps_stop_string_t)dlsym(RTLD_DEFAULT, "ps_tool_stop_string");
    stop_current_function =
        (ps_stop_current_t)dlsym(RTLD_DEFAULT, "ps_tool_stop_current");
    set_parameter_function =
        (ps_set_parameter_t)dlsym(RTLD_DEFAULT, "ps_tool_set_parameter");
    dynamic_phase_start_function = (ps_dynamic_phase_start_t)dlsym(
            RTLD_DEFAULT, "ps_tool_dynamic_phase_start");
    dynamic_phase_stop_function = (ps_dynamic_phase_stop_t)dlsym(
            RTLD_DEFAULT, "ps_tool_dynamic_phase_stop");
    create_counter_function = (ps_create_counter_t)dlsym(
            RTLD_DEFAULT, "ps_tool_create_counter");
    sample_counter_function = (ps_sample_counter_t)dlsym(
            RTLD_DEFAULT, "ps_tool_sample_counter");
    set_metadata_function =
        (ps_set_metadata_t)dlsym(RTLD_DEFAULT, "ps_tool_set_metadata");
    get_timer_data_function = (ps_get_timer_data_t)dlsym(
            RTLD_DEFAULT, "ps_tool_get_timer_data");
    get_counter_data_function = (ps_get_counter_data_t)dlsym(
            RTLD_DEFAULT, "ps_tool_get_counter_data");
    get_metadata_function = (ps_get_metadata_t)dlsym(
            RTLD_DEFAULT, "ps_tool_get_metadata");
    free_timer_data_function = (ps_free_timer_data_t)dlsym(
            RTLD_DEFAULT, "ps_tool_free_timer_data");
    free_counter_data_function = (ps_free_counter_data_t)dlsym(
            RTLD_DEFAULT, "ps_tool_free_counter_data");
    free_metadata_function = (ps_free_metadata_t)dlsym(
            RTLD_DEFAULT, "ps_tool_free_metadata");
#endif
    perfstubs_initialized = PERFSTUBS_SUCCESS;
    /* Increment the number of tools */
    num_tools_registered = 1;
}

char * ps_make_timer_name_(const char * file,
        const char * func, int line) {
    /* The length of the line number as a string is floor(log10(abs(num))) */
    int string_length = (strlen(file) + strlen(func) + floor(log10(abs(line))) + 12);
    char * name = calloc(string_length, sizeof(char));
    sprintf(name, "%s [{%s} {%d,0}]", func, file, line);
    return (name);
}

// used internally to the class
static inline void ps_register_thread_internal(void) {
    if (pthread_getspecific(key) == NULL) {
        if (register_thread_function != NULL) {
            register_thread_function();
            pthread_setspecific(key, (void*)1UL);
        }
    }
}

/* Initialization */
void ps_initialize_(void) {
    /* Only do this once */
    if (perfstubs_initialized != PERFSTUBS_UNKNOWN) {
        return;
    }
    initialize_library();
    if (initialize_function != NULL) {
        initialize_function();
        (void) pthread_once(&key_once, make_key);
        pthread_setspecific(key, (void*)1UL);
    }
}

void ps_finalize_(void) {
    if (finalize_function != NULL)
        finalize_function();
}

void ps_pause_measurement_(void) {
    if (pause_measurement_function != NULL)
        pause_measurement_function();
}

void ps_resume_measurement_(void) {
    if (resume_measurement_function != NULL)
        resume_measurement_function();
}

void ps_register_thread_(void) {
    ps_register_thread_internal();
}

void* ps_timer_create_(const char *timer_name) {
    ps_register_thread_internal();
    void ** objects = (void **)calloc(num_tools_registered, sizeof(void*));
    if (timer_create_function != NULL)
        objects = (void *)timer_create_function(timer_name);
    return (void*)(objects);
}

void ps_timer_create_fortran_(void ** object, const char *timer_name) {
    *object = ps_timer_create_(timer_name);
}

void ps_timer_start_(void *timer) {
    ps_register_thread_internal();
    void ** objects = (void **)timer;
    if (timer_start_function != NULL && objects != NULL)
        timer_start_function(objects);
}

void ps_timer_start_fortran_(void **timer) {
    ps_timer_start_(*timer);
}

void ps_timer_stop_(void *timer) {
    void ** objects = (void **)timer;
    if (timer_stop_function != NULL &&
            objects != NULL)
        timer_stop_function(objects);
}

void ps_timer_stop_fortran_(void **timer) {
    ps_timer_stop_(*timer);
}

void ps_start_string_(const char *timer_name) {
    ps_register_thread_internal();
    if (start_string_function != NULL)
        start_string_function(timer_name);
}

void ps_stop_string_(const char *timer_name) {
    if (stop_string_function != NULL)
        stop_string_function(timer_name);
}

void ps_stop_current_(void) {
    if (stop_current_function != NULL)
        stop_current_function();
}

void ps_set_parameter_(const char * parameter_name, int64_t parameter_value) {
    if (set_parameter_function != NULL)
        set_parameter_function(parameter_name, parameter_value);
}

void ps_dynamic_phase_start_(const char *phase_prefix, int iteration_index) {
    if (dynamic_phase_start_function != NULL)
        dynamic_phase_start_function(phase_prefix, iteration_index);
}

void ps_dynamic_phase_stop_(const char *phase_prefix, int iteration_index) {
    if (dynamic_phase_stop_function != NULL)
        dynamic_phase_stop_function(phase_prefix, iteration_index);
}

void* ps_create_counter_(const char *name) {
    ps_register_thread_internal();
    void ** objects = (void **)calloc(num_tools_registered, sizeof(void*));
    if (create_counter_function != NULL)
        objects = (void*)create_counter_function(name);
    return (void*)(objects);
}

void ps_create_counter_fortran_(void ** object, const char *name) {
    *object = ps_create_counter_(name);
}

void ps_sample_counter_(void *counter, const double value) {
    void ** objects = (void **)counter;
    if (sample_counter_function != NULL &&
            objects != NULL)
        sample_counter_function(objects, value);
}

void ps_sample_counter_fortran_(void **counter, const double value) {
    ps_sample_counter_(*counter, value);
}

void ps_set_metadata_(const char *name, const char *value) {
    ps_register_thread_internal();
    if (set_metadata_function != NULL)
        set_metadata_function(name, value);
}

void ps_dump_data_(void) {
    if (dump_data_function != NULL)
        dump_data_function();
}

void ps_get_timer_data_(ps_tool_timer_data_t *timer_data) {
    if (get_timer_data_function != NULL)
        get_timer_data_function(timer_data);
}

void ps_get_counter_data_(ps_tool_counter_data_t *counter_data) {
    if (get_counter_data_function != NULL)
        get_counter_data_function(counter_data);
}

void ps_get_metadata_(ps_tool_metadata_t *metadata) {
    if (get_metadata_function != NULL)
        get_metadata_function(metadata);
}

void ps_free_timer_data_(ps_tool_timer_data_t *timer_data) {
    if (free_timer_data_function != NULL)
        free_timer_data_function(timer_data);
}

void ps_free_counter_data_(ps_tool_counter_data_t *counter_data) {
    if (free_counter_data_function != NULL)
        free_counter_data_function(counter_data);
}

void ps_free_metadata_(ps_tool_metadata_t *metadata) {
    if (free_metadata_function != NULL)
        free_metadata_function(metadata);
}

