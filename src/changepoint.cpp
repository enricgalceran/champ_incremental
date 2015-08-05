#include "changepoint/changepoint.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

void cp_params_t_print(const cp_params_t * params)
{
    printf("len_mean: %f\n", params->len_mean);
    printf("len_sigma: %f\n", params->len_sigma);
    printf("min_seg_len: %d\n", params->min_seg_len);
    printf("max_particles: %d\n", params->max_particles);
    printf("resamp_particles: %d\n", params->resamp_particles);
}

cp_model_segment_t *cp_model_segment_t_copy(const cp_model_segment_t *s)
{
    assert(s != NULL);

    cp_model_segment_t *sc = (cp_model_segment_t*)calloc(1, sizeof(cp_model_segment_t));

    sc->model_name = strdup(s->model_name); // FIXME: possible leak?
//    std::cout << "Copied model segment name: [" << sc->model_name << "], ptr=" << (void*)(sc->model_name) << std::endl;

    sc->param_names = zarray_create(sizeof(const char*));
    zarray_ensure_capacity(sc->param_names, zarray_size(s->param_names));
    for (int i = 0; i < zarray_size(s->param_names); i++) {
        const char *param_name = NULL; zarray_get(s->param_names, i, &param_name);
        param_name = strdup(param_name);
        zarray_add(sc->param_names, &param_name);
    }

    sc->model_params = zarray_create(sizeof(double));
    zarray_ensure_capacity(sc->model_params, zarray_size(s->model_params));
    for (int i = 0; i < zarray_size(s->model_params); i++) {
        double model_param; zarray_get(s->model_params, i, &model_param);
        zarray_add(sc->model_params, &model_param);
    }

    sc->seg_stats = zarray_create(sizeof(zarray_t*));
    zarray_ensure_capacity(sc->seg_stats, zarray_size(s->seg_stats));
    for (int i = 0; i < zarray_size(s->seg_stats); i++) {
        zarray_t *dp = NULL; zarray_get(s->seg_stats, i, &dp);
        zarray_t *dpc = zarray_copy(dp);
        zarray_add(sc->seg_stats, &dpc);
    }

    sc->first_point = s->first_point;
    sc->last_point = s->last_point;

    return sc;
}

void cp_model_segment_t_print(const cp_model_segment_t *s)
{
    assert(s != NULL);
    assert(s->model_name != NULL);
    assert(s->model_params != NULL);
    assert(s->param_names != NULL);
    assert(s->seg_stats != NULL);
    printf("Model: %s    Length: %d\n", s->model_name, s->last_point - s->first_point + 1);
    printf("Start: %d\n", s->first_point);
    printf("End: %d\n", s->last_point);
    for (int j = 0; j < zarray_size(s->model_params); j++) {
        char *pname = NULL; zarray_get(s->param_names, j, &pname);
        double p; zarray_get(s->model_params, j, &p);
        printf("    %s: %f\n", pname, p);
    }
}

void cp_model_segment_t_destroy(cp_model_segment_t *s)
{
    assert(s != NULL);
    assert(s->model_name != NULL);
    assert(s->param_names != NULL);
    assert(s->model_params != NULL);
    assert(s->seg_stats != NULL);

    // Ugly casts to void* required by C++ compiler
//    printf("cp_modle_segment_t_destroy, freeing s->model_name=[%s], ptr=%p\n", s->model_name, ((void*)s->model_name));
    free((void*)s->model_name);
    zarray_vmap(s->param_names, (void (*)())free);
    zarray_destroy(s->param_names);
    zarray_destroy(s->model_params);
    zarray_vmap(s->seg_stats, (void (*)())zarray_destroy);
    zarray_destroy(s->seg_stats);
    free(s);
}

double changepoint_get_seg_param_by_name(const cp_model_segment_t *s, const char *param_name)
{
    // sanity checks
    assert(NULL != s);
    assert(NULL != s->param_names);
    assert(NULL != s->model_params);
    assert(zarray_size(s->param_names) == zarray_size(s->model_params));
    assert(NULL != param_name);

    for (int i = 0; i < zarray_size(s->param_names); i++) {
        const char *cpn = NULL; zarray_get(s->param_names, i, &cpn);
        if (std::string(cpn) == std::string(param_name)) {
            double pv; zarray_get(s->model_params, i, &pv);
            return pv;
        }
    }

    assert(0 && "param name not found");
    return 42;
}

static
CPDataVector pod2cpp(const zarray_t *data_pts)
{
    // Convert from zarray to std::vector
    CPDataVector cpp_data_pts(zarray_size(data_pts));
    for (int i = 0; i < (int)(cpp_data_pts.size()); i++) {
        void *pod_pt; zarray_get(data_pts, i, &pod_pt);
        cpp_data_pts[i] = pod_pt;
    }
    return cpp_data_pts;
}

static
zarray_t *cpp2pod(const std::vector<ModelSegment>& cpp_segments)
{
    zarray_t *segments = zarray_create(sizeof(cp_model_segment_t*));
    for (std::vector<ModelSegment>::const_iterator i = cpp_segments.begin(); i != cpp_segments.end(); ++i) {
        cp_model_segment_t *s = (cp_model_segment_t*)calloc(1, sizeof(cp_model_segment_t));

        s->model_name = strdup(i->model_name.c_str()); // FIXME: possible leak?
//        std::cout << "Added model segment name: [" << s->model_name << "], ptr=" << (void*)(s->model_name) << std::endl;

        s->param_names = zarray_create(sizeof(const char*));
        zarray_ensure_capacity(s->param_names, i->param_names.size());
        for (std::vector<std::string>::const_iterator j = i->param_names.begin(); j != i->param_names.end(); ++j) {
            const char *param_name = strdup(j->c_str());
            zarray_add(s->param_names, &param_name);
        }

        s->model_params = zarray_create(sizeof(double));
        zarray_ensure_capacity(s->model_params, i->model_params.size());
        for (std::vector<double>::const_iterator j = i->model_params.begin(); j != i->model_params.end(); ++j) {
            double model_param = *j;
            zarray_add(s->model_params, &model_param);
        }

        s->seg_stats = zarray_create(sizeof(zarray_t*));
        zarray_ensure_capacity(s->seg_stats, i->seg_stats.size());
        for (std::vector< std::vector<double> >::const_iterator j = i->seg_stats.begin(); j != i->seg_stats.end(); ++j) {
            zarray_t *pt = zarray_create(sizeof(double));
            for (std::vector<double>::const_iterator k = j->begin(); k != j->end(); ++k) {
                double elm = *k;
                zarray_add(pt, &elm);
            }
            zarray_add(s->seg_stats, &pt);
        }

        s->first_point = i->first_point;
        s->last_point = i->last_point;

        zarray_add(segments, &s);
    }
    return segments;
}

// Returns: zarray_t* of <cp_model_segment_t*>
zarray_t *changepoint_detect_batch(const char *model_class_name,
                                   const cp_params_t *cp_params,
                                   const zarray_t *data_pts,
                                   cp_data_free_fcn data_free_fcn,
                                   cp_data_copy_fcn data_copy_fcn,
                                   void *user)
{
    if (NULL == model_class_name ||
        NULL == cp_params ||
        NULL == data_pts ||
        zarray_size(data_pts) < cp_params->min_seg_len) {
        return NULL;
    }

    // Convert POD input to C++ input
    CPDataVector cpp_data_pts = pod2cpp(data_pts);
    CPParams cpp_cp_params;
    memcpy(&cpp_cp_params, cp_params, sizeof(cp_params_t)); // sizeof(CPParams) == sizeof(cp_params_t)

    // Get C++ output
        CPDetector c(cpp_data_pts, cpp_cp_params, model_class_name, data_free_fcn, data_copy_fcn, user);
        std::vector<ModelSegment> cpp_segments = c.detectChangepoints(model_class_name, user);

    // Test incremental interface
//        std::vector<DataPoint> initial_data;
//        initial_data.push_back(cpp_data_pts[0]);
//        CPDetector c(initial_data, cpp_cp_params, model_class_name, user);
//        std::vector<ModelSegment> cpp_segments;
//        for (size_t i = 1; i < cpp_data_pts.size(); ++i) {
//            cpp_segments = c.detectChangepointsIncrementally(model_class_name,
//                                                             cpp_data_pts[i].point,
//                                                             user);
//        }

    // Convert C++ output to POD
    return cpp2pod(cpp_segments); // FIXME: possible leak?
}

zarray_t *changepoint_fit_class(const char *model_class_name,
                                const cp_params_t *cp_params,
                                const zarray_t *data_pts,
                                cp_data_free_fcn data_free_fcn,
                                cp_data_copy_fcn data_copy_fcn,
                                void *user)
{
    printf("Performing class fitting...\n");
    if (NULL == model_class_name ||
        NULL == cp_params ||
        NULL == data_pts ||
        zarray_size(data_pts) < cp_params->min_seg_len) {
        return NULL;
    }

    // Convert POD input to C++ input
    CPDataVector cpp_data_pts = pod2cpp(data_pts);
    CPParams cpp_cp_params;
    memcpy(&cpp_cp_params, cp_params, sizeof(cp_params_t)); // sizeof(CPParams) == sizeof(cp_params_t)

//    printf("C params:\n");
//    cp_params_t_print(cp_params);
//    std::cout << "C++ params: " << std::endl <<
//                 "len_mean" << cpp_cp_params.len_mean << std::endl <<
//                 "len_sigma" << cpp_cp_params.len_sigma << std::endl <<
//                 "min_seg_len" << cpp_cp_params.min_seg_len << std::endl <<
//                 "max_particles" << cpp_cp_params.max_particles << std::endl <<
//                 "resamp_particles" << cpp_cp_params.resamp_particles << std::endl;

    // Get C++ output
    CPDetector c(cpp_data_pts, cpp_cp_params, model_class_name, data_free_fcn, data_copy_fcn, user);
    std::vector<ModelSegment> cpp_segments = c.fitClass(model_class_name, user);

    // Convert C++ output to POD
    return cpp2pod(cpp_segments); // FIXME: possible leak?
}


CPDetector *changepoint_detector_create(const char *model_class_name,
                                        const cp_params_t *cp_params,
                                        void *initial_point,
                                        cp_data_free_fcn data_free_fcn,
                                        cp_data_copy_fcn data_copy_fcn,
                                        void *user)
{
    if (NULL == model_class_name ||
        NULL == cp_params ||
        NULL == initial_point ) {
        return NULL;
    }

    // Convert POD input to C++ input
    CPDataVector cpp_data_pts;
    cpp_data_pts.push_back(initial_point);

    CPParams cpp_cp_params;
    memcpy(&cpp_cp_params, cp_params, sizeof(cp_params_t)); // sizeof(CPParams) == sizeof(cp_params_t)

    CPDetector *d = new CPDetector(cpp_data_pts, cpp_cp_params, model_class_name, data_free_fcn, data_copy_fcn, user);
    return d;
}

// Returns: zarray_t* of <cp_model_segment_t*>
zarray_t *changepoint_detector_update(CPDetector* d,
                                      const char *model_class_name,
                                      const cp_params_t *cp_params,
                                      void *new_point,
                                      int max_data_points,
                                      void *user)
{
    std::vector<ModelSegment> cpp_segments = d->detectChangepointsInc(model_class_name, new_point, max_data_points, user);
    return cpp2pod(cpp_segments);
}

zarray_t *changepoint_detector_get_data(const CPDetector* d)
{
    CPDataVector data = d->getData();
    zarray_t *pod = zarray_create(sizeof(void*));
    zarray_ensure_capacity(pod, data.size());
    for (CPDataVector::const_iterator it = data.begin(); it != data.end(); ++it) {
        void *dp = *it;
        zarray_add(pod, &dp);
    }
    return pod;
}

void changepoint_detector_destroy(CPDetector* d)
{
    delete d;
}
