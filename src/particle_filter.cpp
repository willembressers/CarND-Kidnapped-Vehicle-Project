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
 
	// push the particles and weights to the vector
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

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

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