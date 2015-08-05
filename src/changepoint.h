/*
 * changepoint.h
 *
 * Bayesian changepoint detection library.
 * C wrapper around Scott Niekum's CHAMP ROS package
 * (written in C++).
 */

#ifndef CHANGEPOINT_H_
#define CHANGEPOINT_H_

#include "common/zarray.h"

#ifdef __cplusplus
#include "cp_detector.h"
extern "C" {
#else
  typedef struct CPDetector CPDetector;
#endif

typedef void  (*cp_data_free_fcn)(void*);
typedef void* (*cp_data_copy_fcn)(const void*);

typedef struct {
    double len_mean;
    double len_sigma;
    int    min_seg_len;
    int    max_particles;
    int    resamp_particles;
} cp_params_t;

void cp_params_t_print(const cp_params_t *params);

typedef struct {
    const char *model_name;
    zarray_t *param_names;  // <const char*>
    zarray_t *model_params; // <double>
    zarray_t *seg_stats;    // <zarray_t<double>*>
    int first_point;
    int last_point;
} cp_model_segment_t;

cp_model_segment_t *cp_model_segment_t_copy(const cp_model_segment_t *s);
void cp_model_segment_t_print(const cp_model_segment_t *s);
void cp_model_segment_t_destroy(cp_model_segment_t *s);

double changepoint_get_seg_param_by_name(const cp_model_segment_t *s, const char *param_name);
//const char *changepoint_get_segment_policy(const cp_model_segment_t *s);

// Performs Bayesian change point detection on time-series data (batch mode)
// Input:
//     model_class_name: class of models to be fitted to data
//     cp_params: detection parameters
//     data_pts: zarray_t* of cp_data_point_t*
// Returns: zarray_t* of cp_model_segment_t*
// Caller must handle cleanup (use cp_model_segment_t_destroy())
zarray_t *changepoint_detect_batch(const char *model_class_name,
                                   const cp_params_t *cp_params,
                                   const zarray_t *data_pts,
                                   cp_data_free_fcn data_free_fcn,
                                   cp_data_copy_fcn data_copy_fcn,
                                   void *user);

// Fits all models of a class to a data series.
// Returns: zarray_t* of cp_model_segment_t*, with one
// segment per model in the class
zarray_t *changepoint_fit_class(const char *model_class_name,
                                const cp_params_t *cp_params,
                                const zarray_t *data_pts,
                                cp_data_free_fcn data_free_fcn,
                                cp_data_copy_fcn data_copy_fcn,
                                void *user);

CPDetector *changepoint_detector_create(const char *model_class_name,
                                        const cp_params_t *cp_params,
                                        void *initial_point,
                                        cp_data_free_fcn data_free_fcn,
                                        cp_data_copy_fcn data_copy_fcn,
                                        void *user);

// Performs Bayesian change point detection on time-series data (online/incremental mode)
// Input:
//     model_class_name: class of models to be fitted to data
//     cp_params: detection parameters
//     new_point: new data point to be considered in the detection
// Returns: zarray_t* of cp_model_segment_t*
// Caller must handle cleanup (use cp_model_segment_t_destroy())
zarray_t *changepoint_detector_update(CPDetector* d,
                                      const char *model_class_name,
                                      const cp_params_t *cp_params,
                                      void *new_point,
                                      int max_data_points,
                                      void *user);

// Returns a shallow copy (void*'s) of the detector's data
// Returned pointers might be invalidated after a call to _update or _destroy
// Caller must call zarray_destroy on it
zarray_t *changepoint_detector_get_data(const CPDetector* d);

void changepoint_detector_destroy(CPDetector* d);

#ifdef __cplusplus
}
#endif

#endif /* CHANGEPOINT_H_ */
