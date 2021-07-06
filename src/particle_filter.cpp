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
  std::normal_distribution<double> norm_x(x, std[0]);
  std::normal_distribution<double> norm_y(y, std[1]);
  std::normal_distribution<double> norm_theta(theta, std[2]);
  
  // Generate particles
  for (int i=0; i < num_particles; i++) {

    // create the particle
    Particle particle;
    particle.id = i;
    particle.x = norm_x(gen); 
    particle.y = norm_y(gen);
    particle.theta = norm_theta(gen);
    particle.weight = 1.0;
 
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
  std::normal_distribution<double> norm_x(0.0, std_pos[0]);
  std::normal_distribution<double> norm_y(0.0, std_pos[1]);
  std::normal_distribution<double> norm_theta(0.0, std_pos[2]);
  
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
  for (unsigned int i = 0; i < observations.size(); i++) {
    
    // init distance to maximum possible
    double minDistance = std::numeric_limits<double>::max();
    
    // loop over all predictions
    for (unsigned int j=0; j < predicted.size(); j++){
      
      // calculate the distance between each observation and prediction
      double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      
      // find the nearest neighbour
      if (distance < minDistance){
        minDistance = distance;
        observations[i].id = predicted[j].id;
      }
    }
  }
}

// ============================================================================
// helper: calculate the multi-variate Gaussian distribution
// ============================================================================
double multiv_prob_gaussian(double sig_x, double sig_y, double x_obs, double y_obs, double mu_x, double mu_y) {
  
  // calculate normalization term
  double gaussian_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2))) + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

  // calculate weight using normalization terms and exponent
  double weight = gaussian_norm * exp(-exponent);
    
  return weight;
}

// ============================================================================
// Update weight
// ============================================================================
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  
  // Loop over all particles
  for (unsigned int i = 0; i < particles.size(); ++i) {
    
    // Loop over all observations and transform it's perspective
    std::vector<LandmarkObs> observations_transformed;
    for (unsigned int j = 0; j< observations.size(); j++) {
      LandmarkObs temp;
      temp.x = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
      temp.y = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
      observations_transformed.push_back(temp);
    }

    // Loop over the landmarks and find within the sensor range only
    std::vector<LandmarkObs> predicted;
    std::vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;
    for (unsigned int j = 0; j < landmarks.size(); j++) {
      double distance = dist(landmarks[j].x_f, landmarks[j].y_f, particles[i].x, particles[i].y);
      if(distance <= sensor_range){
        LandmarkObs pred;
        pred.id = landmarks[j].id_i;
        pred.x = landmarks[j].x_f;
        pred.y = landmarks[j].y_f;
        predicted.push_back(pred);
      }
    }

    // Assosiate data
    dataAssociation(predicted, observations_transformed);

    // Calculate the weight
    double weight = 1.0;
    for (unsigned int j = 0; j< observations_transformed.size(); j++) {
      double mu_x = 0;
      double mu_y = 0;
      for (unsigned int k=0; k < predicted.size(); k++){
        if(predicted[k].id == observations_transformed[j].id){
          mu_x = predicted[k].x;
          mu_y = predicted[k].y;
        }
      }
      weight *= multiv_prob_gaussian(std_landmark[0], std_landmark[1], observations_transformed[j].x, observations_transformed[j].y, mu_x, mu_y);
    }

    // update the particle + weights, weight
    particles[i].weight = weight;
    weights[i] = weight;
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