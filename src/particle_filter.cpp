/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

// ============================================================================
// initialization
// ============================================================================
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set the number of particles
  num_particles = 100;
  
  // Gaussian random noise generator
  std::default_random_engine gen;
  
  // initialize the standard distribution and theta
  std::normal_distribution<double> init_x(0, std[0]);
  std::normal_distribution<double> init_y(0, std[1]);
  std::normal_distribution<double> init_theta(0, std[2]);
  
  // Generate particles
  for (int i=0; i<num_particles; i++){
    Particle particle;
    particle.id = i;
    particle.x = x; 
    particle.y = y;
    particle.theta = theta;
    particle.weight = 1.0;
   
    // assign randon gaussian noise
    particle.x += init_x(gen); 
    particle.y += init_y(gen);
    particle.theta += init_theta(gen);
 
    // push the particles and weights to the vectors
    particles.push_back(particle);
    weights.push_back(particle.weight);
  }
  
  // flag that the particles are initialized
  is_initialized = true;
}

// ============================================================================
// Prediction step
// ============================================================================
void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  // Gaussian random noise generator
  std::default_random_engine gen;

  // define normal distributions
  std::normal_distribution<double> norm_x(0, std_pos[0]);
  std::normal_distribution<double> norm_y(0, std_pos[1]);
  std::normal_distribution<double> norm_theta(0, std_pos[2]);
  
  // loop over the particles
  for (int i = 0; i < num_particles; i++) {

    if (fabs(yaw_rate) < 0.00001) {
      
      // driving straight
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
      
    } else {
      
      // vehicle is turning
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    
    // update with some random gaussian noise
    particles[i].x += norm_x(gen);
    particles[i].y += norm_y(gen);
    particles[i].theta += norm_theta(gen);
  }
}

// ============================================================================
// Map landmark positions
// ============================================================================
void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  // loop over all observations
  for (unsigned int i = 0; i< observations.size(); i++){
    
    // init distance to maximum possible
    double minDistance = std::numeric_limits<double>::max();
    int index_map = -1;
    
    // loop over all predictions
    for (unsigned int j=0; j<predicted.size(); j++){
      
      // calculate the distance between each observation and prediction
      double x_distance = observations[i].x - predicted[j].x;
      double y_distance = observations[i].y - predicted[j].y;
      double distance = x_distance * x_distance + y_distance * y_distance;
      
      // find the nearest neighbour
      if (distance < minDistance){
        minDistance = distance;
        index_map =  predicted[j].id;
      }
    }
    
    // add the nearest prediction index
    observations[i].id = index_map;
  }
}

// ============================================================================
// helper: calculate the multi-variate Gaussian distribution
// ============================================================================
double multiv_prob_gaussian(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  
  double gaussian_norm;
  double exponent;
  double weight;
  
  // calculate normalization term
  gaussian_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2))) + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

  // calculate weight using normalization terms and exponent
  weight = gaussian_norm * exp(-exponent);
    
  return weight;
}

// ============================================================================
// Update weight
// ============================================================================
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  
  // Read map data
  std::vector<Map::single_landmark_s> landmark_list = map_landmarks.landmark_list;
  double land_x;
  double land_y;
  double max_val = 2 * sensor_range;
  
  // loop over all particles
  for (unsigned int i = 0; i < particles.size(); ++i) {
     
    // get the particle object
    Particle particle = particles[i];
    
    // initialize the particle probability to 1.0
    double prob = 1.0;

    // loop over all observations
    for (unsigned int j = 0; j < observations.size(); j++) {
      
      // transforms observations from vehicle's coordinate system to map's coordinate system
      double x_map = particle.x + (cos(particle.theta) * observations[j].x) - (sin(particle.theta) * observations[j].y);
      double y_map = particle.y + (sin(particle.theta) * observations[j].x) + (cos(particle.theta) * observations[j].y);

      // loop over the landmarks
      for (unsigned int k = 0; k < landmark_list.size(); k++) {
                
        // Calculate distance between particle and landmarks
        double local_land_x = landmark_list[k].x_f;
        double local_land_y = landmark_list[k].y_f;

        // dist function
        double distance = dist(x_map, y_map, local_land_x, local_land_y);
                
        if ((distance <= sensor_range) && (distance <= max_val)){ 
          
          // Calculate multivariate Gaussian normal distribution
          land_x = local_land_x;
          land_y = local_land_y;
          max_val = distance;
          prob = multiv_prob_gaussian(std_landmark[0], std_landmark[1], x_map, y_map, land_x, land_y);
          
          // update the particle weights with the gaussian distribution
          particles[i].weight = prob;
          weights[i] = prob;       
        }     
      }   
    }  
  }
}

// ============================================================================
// resample
// ============================================================================
void ParticleFilter::resample() {

  // Gaussian random noise generator
  std::default_random_engine gen;

  // Generate discrete distribution based on the weights
  std::discrete_distribution<> d(weights.begin(), weights.end());
  
  // collect the resamples particles here
  std::vector<Particle> resampled_particles;

  // loop over all particles
  for (int n = 0; n < num_particles; ++n) {

    // generate new particles based on the discrete weights
    Particle particle = particles[d(gen)];
    
    // collect the new particles in the dedicated vector
    resampled_particles.push_back(particle);
  }
  
  // overwrite the old vector with the new one
  particles = resampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}