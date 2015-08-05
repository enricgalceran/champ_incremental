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

// Enric:
// fixed memory leaks (should add smart pointers)
// added

#include "cp_detector.h"

#include "in_lane.h"
#include "gauss1D.h"
#include "policies.h"

#define SQRT2PI 2.50662827463

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <iterator>
using namespace std;

#define CP_DETECTOR_VERBOSE

//namespace changepoint{

static
bool createHelper(const string& class_type, ModelFitter* &c, void *user)
{
    try {
        //c = c_loader.createClassInstance(class_type);
        //c = c_loader.createUnmanagedInstance(class_type);

        // FIXME: likely memory leak in block below
        if (class_type == "changepoint/Gauss1DFitter")
            c = new Gauss1DFitter;
        else if (class_type == "changepoint/InLaneFitter")
            c = new InLaneFitter;
        else if (class_type == "changepoint/PoliciesFitter")
            c = new PoliciesFitter(user);
        else
            throw std::runtime_error(std::string("Invalid model class: [") + class_type + "]");
    }
    //catch (pluginlib::PluginlibException& ex) {
    catch (std::exception& ex){
//        ROS_ERROR("Classifer plugin failed to load! Error: %s", ex.what());
        ROS_ERROR("Model class failed to load! Error: %s\n", ex.what());
        exit(EXIT_FAILURE);
    }

    return true;
}


CPDetector::CPDetector(const CPDataVector& data_pts,
                       const /*changepoint::*/CPParams& cp_params,
                       const std::string& class_type,
                       CPDataFreeFcn data_free_fcn,
                       CPDataCopyFcn data_copy_fcn,
                       void *user)
:
    once(false),
    model_class(class_type),
    d_len(0),
    data(),
    size_Q(0),
    max_MAP(0),
    l_prior(0),
    particles(),
    max_path_indices(),
    max_path_models(),
    prev_max_MAP(),
    model_fitter_ptrs()
{
    assert(!data_pts.empty());

    // Load data
    d_len = data_pts.size();
//    if (d_len > 0) {
//        d_dim =  data_pts[0].point.size();
//        assert(d_dim > 0); // E Galceran, if no elements in point
//    }

//    data = new double*[d_len];
//    for(int i=0; i<d_len; i++){
//        data[i] = new double[d_dim];
//        for(int j=0; j<d_dim; j++){
//            data[i][j] = data_pts[i].point[j];
//        }
//    }

    data = CPDataVector(d_len);
    for (int i=0; i<d_len; i++) {
        data[i] = data_copy_fcn(data_pts[i]);
//        for (int j=0; j<d_dim; j++) {
//            data[i][j] = data_pts[i].point[j];
//        }
    }

    this->data_free_fcn = data_free_fcn;
    this->data_copy_fcn = data_copy_fcn;

    // Set up parameters
    LEN_MEAN = cp_params.len_mean;
    LEN_SIG = cp_params.len_sigma;
    MIN_SEG_LEN = cp_params.min_seg_len;
    MAX_PARTICLES = cp_params.max_particles;
    RESAMP_PARTICLES = cp_params.resamp_particles;

    // Initialize prior
    srand(time(0));
    ModelFitter *temp;
    createHelper(class_type, temp, user);
    size_Q = temp->nModels();      // The number of models
    delete temp;

    max_MAP = log(1.0 / size_Q);
    l_prior = max_MAP;         // A uniform prior for all models

    //Initialize particle filter and Viterbi stats
//    particles.clear();
    for (int i=0; i<size_Q; i++) {
        ModelFitter *mf;
        createHelper(class_type, mf, user); this->model_fitter_ptrs.push_back(mf);

        mf->initialize(i);

        //particles.push_back(*(new Particle(max_MAP,-1,*mf))); // E Galceran: original line, leak
        Particle particle(max_MAP,-1,*mf);
        particles.push_back(particle);
    }
//    vector<int> max_path_indices;
//    vector<ModelFitter*> max_path_models;
//    prev_max_MAP = queue<double>();
}

CPDetector::~CPDetector()
{
    for (std::vector<ModelFitter*>::const_iterator it = max_path_models.begin(); it != max_path_models.end(); ++it)
        if (*it != NULL) delete *it;

//    for(int i=0; i<d_len; i++)
//        delete[] data[i];
//    delete[] data;

    for (CPDataVector::const_iterator it = data.begin(); it != data.end(); ++it)
        data_free_fcn(*it);

    for (std::vector<const ModelFitter*>::const_iterator it = model_fitter_ptrs.begin(); it != model_fitter_ptrs.end(); ++it)
        delete *it;
}


// Gaussian cdf calculation
double CPDetector::gaussCDF(double t){
    return 0.5 * erfc(-(t-LEN_MEAN)/(LEN_SIG*M_SQRT2));
}


// Calc log length distribution with discretized truncated Gaussian distribution
// NOTE: This only is accurate up to a normalizing constant introduced by truncation.
double CPDetector::logLenPDF(int t){
    // Integrate between t-1 and t, so that the CDF is consistent later for cdf(t) = pdf(t) + pdf(t-1) + ... + pdf(MIN_SEG_LEN)
    return log( gaussCDF(t) - gaussCDF(t-1.0) );
}


// Calc cdf of log len dist by truncating normal gauss cdf: cdf(t) - cdf(a)
// NOTE: This only is accurate up to a normalizing constant introduced by truncation.
double CPDetector::logOneMinLenCDF(int t){
    // Truncate at MIN_SEG_LEN
    double p = gaussCDF(t) - gaussCDF(MIN_SEG_LEN-1);
    return log(1-p);
}


// Helper fxn: Sorting fxn for Particles to be used by std::sort
bool MAPSortingFxn(const Particle &i, const Particle &j) { return i.MAP < j.MAP; }


// Calculate the alpha value for SOR resampling.  M is num particles to reasmple down to
double CPDetector::calculateAlpha(int M)
{
    //for(size_t i=1; i<particles.size(); i++)
    //    printf("Part %i:  pos: %i  Norm: %f  MAP: %f\n", i, particles[i].pos, particles[i].nMAP, particles[i].MAP);

    // Implementation of algo from Fearnhead & Clifford (2003)
    // Look for smallest kappa = w_i that satisfies sum[min(w_i/kappa , 1)] <= M
    int N = particles.size();
    for(int i=N-M-1; i<N; i++){   // i is the cutoff index
        int A_i = N-i-1;          // Number of elements > w_i
        double B_i = 0;           // Sum of elements <= w_i
        for(int j=0; j<=i; j++)
            B_i += particles[j].nMAP;

        double kappa = particles[i].nMAP;
        if(kappa == 0) continue;                  // Account for underflow
        double stat = (1/kappa)*B_i + A_i;        // Check if sum[min(w_i/kappa , 1)] <= M
        if(stat <= M){
             double alpha = B_i / (M - A_i);      // Prob mass of w_i <= kappa div by size of set
             //printf("i %i  kappa %f  A_i %i  B_i %f  Stat: %f  M: %i  Alpha: %f\n", i, kappa, A_i, B_i, stat, M, alpha);
             return alpha;
        }
    }

    printf("ERROR: calculateAlpha(): No suitible alpha value found\n");
    return -1;
}


// Normalize the MAP values (into nMAP) of the vector of particles. These are NOT log.
void CPDetector::normalizeMAP()
{
    // Find the max log
    double max_log = particles[0].MAP;
    for(size_t i=1; i<particles.size(); i++){
        if(particles[i].MAP > max_log)
            max_log = particles[i].MAP;
    }

    // Factor out (subtract) the max log and un-log them. Some may still underflow, but they are basically zero anyway.
    for(size_t i=0; i<particles.size(); i++)
        particles[i].nMAP = exp(particles[i].MAP - max_log);

    // Normalize
    double total = 0;
    for(size_t i=0; i<particles.size(); i++){
        double curr = particles[i].nMAP;
        if(!isnan(curr) && curr != -INFINITY)
            total += curr;
    }
    for(size_t i=0; i<particles.size(); i++){
        double norm = particles[i].nMAP / total;
        particles[i].nMAP = norm;
    }

}


// Resample down to resamp_particles if there are more than max_particles, using Stratified Optimal Resampling.
void CPDetector::resampleParticles(int max_particles, int resamp_particles)
{
    // Get ready to perform resampling if there are too many particles
    if((int)particles.size() > max_particles){
        // Normalize MAP values of all particles
        normalizeMAP();

        // Throw out particles with a p->nMAP value of -INFINITY or NAN
        int nerased = 0;
        for(std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); ) {
            if(isnan(iter->nMAP) || iter->nMAP == -INFINITY) {
                iter = particles.erase(iter);
                nerased++;
            }
            else ++iter;
        }
        cout << "$$$ Normalization erased " << nerased << " particles" << endl;
    }

    // Only continue with resampling if the renormalizing step didn't get rid of enough
    if((int)particles.size() > max_particles){
        // Sort the particles by normalized p->MAP values, smallest to largest
        std::sort(particles.begin(), particles.end(), MAPSortingFxn);

        // Calculate alpha
        double alpha = calculateAlpha(resamp_particles);

        int initial_np = particles.size();

        // Keep particles that have weight >= alpha
        int A=0;
        std::vector<Particle> new_particles;
        for(std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); ) {
            //cout << "\t\titer->nMAP=" << iter->nMAP << endl;
            if(iter->nMAP >= alpha){
                new_particles.push_back(*iter);
                iter = particles.erase(iter);
                A++;
            }
            else ++iter;
        }

        // Choose random u from uniform dist [0, alpha]
        double u = (rand() / double(RAND_MAX)) * alpha;

        // Resample using SOR(3)
        for(std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); ) {
            u -= iter->nMAP;
            if(u <= 0){
                new_particles.push_back(*iter);
                iter = particles.erase(iter);
                u = u + alpha;
                A++;
            }
            else ++iter;
        }

        cout << "$$$ Resampling killed " << initial_np - new_particles.size() << " particles, "
                "alpha=" << alpha << endl;

        // If underflow caused too few particles to be resampled, choose randomly
        while(A < resamp_particles){
            int ind = (rand() / double(RAND_MAX)) * particles.size();
            new_particles.push_back(particles[ind]);
            particles.erase(particles.begin() + ind);
            A++;
        }

        // Replace old particles with new resampled set
        particles = new_particles;
    }
}


void CPDetector::detectionStep(int t, const string& class_type, void *user)
{
    // Only create new particle for first time when a CP there would divide data in 2 halves, each of MIN_SEG_LEN
    if(t >= (2*MIN_SEG_LEN)-1){
        double prev = prev_max_MAP.front();
        prev_max_MAP.pop();
        //printf("Adding particles for p=%i with prev: %f\n",t-MIN_SEG_LEN,prev);

        // Create new particles for a changepoint at time t-MIN_SEG_LEN, one for each model in Q
        for(int i=0; i<size_Q; i++){
            ModelFitter *mf;
            createHelper(class_type, mf, user); this->model_fitter_ptrs.push_back(mf);
            mf->initialize(i);
            // Particle *new_p = new Particle(prev,t-MIN_SEG_LEN,*mf); // E Galceran: original lines, leak
            // particles.push_back(*new_p);
            Particle new_p(prev,t-MIN_SEG_LEN,*mf);
            particles.push_back(new_p);
        }
    }

    // Compute fit probs for all particles
    for(size_t i=0; i<particles.size(); i++){
        Particle *p = &(particles[i]);
        int seg_len = t-p->pos;

        // Fit the model and calc the data likelihood
        p->fitter->fitSegment(data, p->pos, t, user);

        // p_tjq is the prob of the CP **PRIOR** to time t occuring at j (p_pos)
        // p->MAP is the prob of the MAP CP occuring at time j (p_pos)
        double p_tjq;
        double p_ME = p->fitter->mp->getModelEvidence();
        if(isnan(p_ME) || p_ME == -INFINITY)
            p_tjq = -INFINITY;
        else
            p_tjq = logOneMinLenCDF(seg_len-1) + p_ME + l_prior + p->prev_MAP;

        if(isnan(p_tjq) || p_tjq == -INFINITY)
            p->MAP = -INFINITY;
        else
            p->MAP = p_tjq + logLenPDF(seg_len) - logOneMinLenCDF(seg_len-1);

        cout << "p_tjq = (" << logOneMinLenCDF(seg_len-1) << ") + p_ME=" << p_ME << " + l_prior=" << l_prior << " + prev_MAP=" << p->prev_MAP << endl;
        cout << "pMAP  = p_tjq=" << p_tjq << " + (" << logLenPDF(seg_len)  << ") - (" << logOneMinLenCDF(seg_len-1) << ")" << endl << endl;

        //printf("t %i   pos %i   model %s   ll %f   evi %f   tjq %f   map %f\n",t,p->pos,p->fitter->mp->getModelName().c_str(),p->fitter->mp->logLikelihood,p->fitter->mp->modelEvidence,p_tjq,p->MAP);

    }

    // Update global stats and Viterbi path
    double max = -INFINITY;
    Particle *max_particle = NULL;
    for(size_t i=0; i<particles.size(); i++){
        if(particles[i].MAP > max){
            max_particle = &(particles[i]);
            max = max_particle->MAP;
        }
    }
    if(max_particle != NULL){
        max_MAP = max_particle->MAP;
        max_path_indices.push_back(max_particle->pos);
        ModelFitter *copy_model;
        createHelper(class_type, copy_model, user);
        copy_model->copyFrom(max_particle->fitter);
        max_path_models.push_back(copy_model);
#ifdef CP_DETECTOR_VERBOSE
        cout << "MAX " << t << "  pos: " << max_particle->pos
             << "  model: " << max_particle->fitter->mp->getModelName()
             << "  map: " << max_particle->MAP << "\n\n";
#endif
    }
    else{
        max_path_indices.push_back(-1);
        max_path_models.push_back(NULL);
    }

    // If there are more than MAX_PARTICLES, resample down to RESAMP_PARTICLES
    resampleParticles(MAX_PARTICLES, RESAMP_PARTICLES);

    // Keep track of MIN_SEG_LEN many of the previous max_MAP values to create particles with later.
    prev_max_MAP.push(max_MAP);

    cout << "MAX_PATH_INDICES:" << endl;
    std::copy(max_path_indices.begin(), max_path_indices.end(), std::ostream_iterator<int>(std::cout, " "));
    cout << endl << endl;
}

vector<ModelSegment> CPDetector::getCurrentViterbiPath()
{
    // Now max_path contains the final Viterbi path, so trace it
    int curr_cp = d_len - 1;
    int path_index = curr_cp - MIN_SEG_LEN + 1;  // This isn't a CP number, but an index into the max_path
    vector<ModelSegment> segments;

//    cout << "size of max_path_indices: " << max_path_indices.size() << endl;
//    cout << "size of max_path_models:  " << max_path_models.size() << endl;
//    cout << "curr_cp=" << curr_cp << endl;
//    cout << "path_index=" << path_index << endl;

#ifdef CP_DETECTOR_VERBOSE
    printf("\nFINAL CHANGEPOINTS:\n");
#endif
    while (curr_cp > -1) {
//        cout << "in loop, path_index=" << path_index << endl;

        if (path_index < 0 || path_index >= static_cast<int>(max_path_indices.size())) break;
        if (max_path_indices[path_index] < 0) break;

        // Print CP info
#ifdef CP_DETECTOR_VERBOSE
        cout << "start: " << (max_path_indices[path_index]+1) <<
                "   model: " << (max_path_models[path_index]->mp->getModelName()) <<
                "   len: " << (curr_cp - (max_path_indices[path_index])) << "\n"; // E. Galceran: comma operator (,) originally present here
        max_path_models[path_index]->mp->printParams();
#endif

        // Compute stats on each final segment
        std::vector< std::vector<double> > stats = max_path_models[path_index]->calcFinalSegStats(data,max_path_indices[path_index],curr_cp);

        // Add a ModelSegment for the ROS message response
        ModelSegment temp;
        temp.model_name = max_path_models[path_index]->mp->getModelName();
        temp.first_point = max_path_indices[path_index]+1;
        temp.last_point = curr_cp;
        max_path_models[path_index]->mp->fillParams(temp);
        temp.seg_stats = stats;
        segments.insert(segments.begin(), temp);

        // Go to previous CP in chain
        curr_cp = max_path_indices[path_index];
        path_index = curr_cp - MIN_SEG_LEN + 1;
    }
    return segments;
}


vector<ModelSegment> CPDetector::detectChangepoints(const string& class_type, void *user)
{
    assert(once == false);

#ifdef CP_DETECTOR_VERBOSE
    cout << "detectChangepoints(), class_type=" << class_type << endl;
#endif

    // Process each time step, approximating the distribution of the most recent changepoint *BEFORE* time t
    // NOTE: By starting at MIN_SEG_LEN-1, we ensure the first segment will be >= MIN_SEG_LEN
    for(int t=MIN_SEG_LEN-1; t<d_len; t++)
        detectionStep(t, class_type, user);

    once = true;

    cout << "NPARTICLES ==== " << particles.size() << endl;

    return getCurrentViterbiPath();
}

std::vector<ModelSegment> CPDetector::detectChangepointsInc(const std::string& class_type,
                                                            void *new_point,
                                                            int max_data_points,
                                                            void *user)
{
    // Add new data point (account for max_data_points)
    data.push_back(this->data_copy_fcn(new_point));
    if (data.size() > static_cast<size_t>(max_data_points)) {
        // Relocate data points and drop old particles
        this->data_free_fcn(data[0]);
        data = std::vector<void*>(data.begin()+1, data.end());
        std::vector<Particle> new_particles;
        new_particles.reserve(this->particles.size());
        typedef std::vector<Particle>::const_iterator It;
        for (It it = this->particles.begin(); it != this->particles.end(); ++it) {
            Particle np = *it;
            np.pos--;
            if (np.pos >= -1)
                new_particles.push_back(np);
        }
        this->particles = new_particles;
    }
    d_len = data.size();

    int t = d_len-1;
    detectionStep(t, class_type, user);
    return getCurrentViterbiPath();
}

vector<ModelSegment> CPDetector::fitClass(const string& class_type, void *user)
{
#ifdef CP_DETECTOR_VERBOSE
    cout << "fitClass(), class_type=" << class_type << endl;
#endif

    ModelFitter *temp;
    createHelper(class_type, temp, user);
    int size_Q = temp->nModels();      // The number of models
    delete temp;

    vector<ModelSegment> fits;

    for(int i=0; i<size_Q; i++) {
        ModelFitter *mf;
        createHelper(class_type, mf, user);
        mf->initialize(i);
        int start = -1;
        int end = this->d_len-1;
        mf->fitSegment(this->data, start, end, user);
        ModelSegment s;
        s.first_point = 0;
        s.last_point = this->d_len-1;
        s.model_name = mf->mp->getModelName();
        mf->mp->fillParams(s);
        s.seg_stats = mf->calcFinalSegStats(data,start,end);
        delete mf;
        fits.push_back(s);
    }

    return fits;
}

//} // end namespace
