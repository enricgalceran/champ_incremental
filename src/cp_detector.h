/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2014, Scott Niekum
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Willow Garage nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************/

/**
  * \author Scott Niekum
  */

// Main changes (Enric Galceran): removed ROS dependency, removed "changepoint" namespace

#ifndef CP_DETECTOR_H_
#define CP_DETECTOR_H_

#include <math.h>
#include <iostream>
#include <string>
#include <typeinfo>
#include <queue>

#include "noros_types.h"

//namespace changepoint{

typedef std::vector<void*> CPDataVector; // TODO: move all into namespace

typedef void  (*CPDataFreeFcn)(void*);
typedef void* (*CPDataCopyFcn)(const void*);

class ModelParams
{
public:
    ModelParams():modelEvidence(-INFINITY), logLikelihood(-INFINITY) {};
    virtual ~ModelParams(){};

    virtual void printParams() = 0;
    virtual void fillParams(ModelSegment &seg) = 0;
    virtual std::string getModelName() = 0;

    double getModelEvidence(){return modelEvidence;}
    double getLogLikelihood(){return logLikelihood;}

    double modelEvidence;  //If params integrated, then Bayes model evidence, otherwise probably BIC
    double logLikelihood;
};


class ModelFitter{
public:
    ModelFitter(){mp=NULL;}
    virtual ~ModelFitter(){};

    virtual void initialize(int model_id) = 0;
    virtual void copyFrom(ModelFitter *rhs) = 0;
    // Fit a model to the segment of start+1 to end. Params+BIC+LL in mp should be filled by fit.
    virtual bool fitSegment(const CPDataVector& data/*double **data*/, const int start, const int end, void *user) = 0;
    virtual int nModels() = 0;  // Returns number of model types for this set of models
    // Calculates arbitrary model-specific stats for a specified data segment under the current model
    virtual std::vector< std::vector<double> > calcFinalSegStats(const CPDataVector& data/*double **data*/, const int start, const int end) = 0;

    ModelParams *mp;
};


class Particle
{
    Particle();
//    void *get_user_data_copy(void *original)
//    {
//        // Deep-copies `original'
//        const policy_model_io_t *original_mio = static_cast<const policy_model_io_t*>(original);
//
//    }

public:
    Particle(double prev_MAP, int pos, ModelFitter &mf)
    {
        this->MAP = -INFINITY;
        this->nMAP = -INFINITY;
        this->prev_MAP = prev_MAP;
        this->pos = pos;
        this->fitter = &mf;
    }
    // Default copy ctor. and assignment op. are fine
    ~Particle(){};

    double MAP;
    double nMAP;
    double prev_MAP;
    int pos;
    ModelFitter *fitter;
};


class CPDetector
{
private:
    CPDetector();

public:
    CPDetector(const CPDataVector& data_pts,
               const /*changepoint::*/CPParams& cp_params,
               const std::string& class_type,
               CPDataFreeFcn data_free_fcn,
               CPDataCopyFcn data_copy_fcn,
               void *user);
    ~CPDetector();

    std::string getModelClass() const
    {
        return this->model_class;
    }

    CPDataVector getData() const
    {
        return this->data;
    }

    std::vector<ModelSegment> detectChangepoints(const std::string& class_type, void *user);
    std::vector<ModelSegment> detectChangepointsInc(const std::string& class_type, void *new_point, int max_data_points, void *user);
    std::vector<ModelSegment> fitClass(const std::string& class_type, void *user);

private:
    double gaussCDF(double t);
    double logLenPDF(int t);
    double logOneMinLenCDF(int t);
    void normalizeMAP();
    double calculateAlpha(int M);
    void resampleParticles(int max_particles, int resamp_particles);

    void detectionStep(int t, const std::string& class_type, void *user);
    std::vector<ModelSegment> getCurrentViterbiPath();

    bool once;
    std::string model_class;
    int d_len;
    CPDataVector data;
    int    size_Q;
    double max_MAP;
    double l_prior;
    std::vector<Particle> particles;
    std::vector<int> max_path_indices;
    std::vector<ModelFitter*> max_path_models;
    std::queue<double> prev_max_MAP;

    // For proper cleanup
    std::vector<const ModelFitter*> model_fitter_ptrs;

    // For managing data points memory
    CPDataFreeFcn data_free_fcn;
    CPDataCopyFcn data_copy_fcn;

    // Detection params
    double LEN_MEAN;                // Mean of segment length gaussian
    double LEN_SIG;                 // Sigma of segment length gaussian
    int MIN_SEG_LEN;                // The minimum length of a segment for model fitting purposes
    int MAX_PARTICLES;              // The most particles to ever keep in the filter
    int RESAMP_PARTICLES;           // The number of particles to resample back down to when resampling
};

//} // end namespace

#endif /* CP_DETECTOR_H_ */
