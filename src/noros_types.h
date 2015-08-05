#ifndef NOROS_TYPES_H_
#define NOROS_TYPES_H_

#include <vector>
#include <string>

// Workaround for ROS macros
#define ROS_ERROR printf

typedef struct {
    double len_mean;
    double len_sigma;
    int    min_seg_len;
    int    max_particles;
    int    resamp_particles;
} CPParams;

typedef struct {
    std::string model_name;
    std::vector<std::string> param_names;
    std::vector<double> model_params;
    std::vector< std::vector<double> > seg_stats;
    int first_point;
    int last_point;
} ModelSegment;

#endif /* NOROS_TYPES_H_ */
