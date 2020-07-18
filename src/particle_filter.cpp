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
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  // Set standard deviations for x, y, and theta
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  
  std::default_random_engine gen;
  // This line creates a normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  for (int i = 0; i < num_particles; ++i) {
    // Declare an object of particle
    Particle particle;
    
    // Sample from these normal distributions like this: 
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    //   where "gen" is the random engine initialized earlier.

    // Print your samples to the terminal.
    particles.push_back(particle);
  }
  
  is_initialized = true;
  
//   for (auto i : particles){
//     std::cout << i.x << ' ' << i.y << ' ' << i.theta << ' ' << std::endl;
//   }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // Set standard deviations for x, y, and theta
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];
  
  std::default_random_engine gen;
  

  // Update each particle
  for(unsigned int i=0; i< particles.size(); ++i){
    
    double x0 = particles[i].x;
    double y0 = particles[i].y;
    double theta0 = particles[i].theta;
    if (fabs(yaw_rate) < 0.0001) {
      particles[i].x = x0 + velocity * delta_t * cos(theta0);
      particles[i].y = y0 + velocity * delta_t * sin(theta0);
    }
    else{
      particles[i].x = x0 + velocity / yaw_rate * (sin(theta0 + yaw_rate * delta_t) - sin(theta0));
      particles[i].y = y0 + velocity / yaw_rate * (cos(theta0) - cos(theta0 + yaw_rate * delta_t));
      particles[i].theta = theta0 + yaw_rate * delta_t;
    }
    
    // Sample from these normal distributions:
    // This line creates a normal (Gaussian) distribution for x, y and theta
    normal_distribution<double> dist_x(particles[i].x, std_x);
    normal_distribution<double> dist_y(particles[i].y, std_y);
    normal_distribution<double> dist_theta(particles[i].theta, std_theta);
    
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
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
  
  for (unsigned int i = 0; i < observations.size(); ++i){
    double x1 = observations[i].x;
    double y1 = observations[i].y;
    double min_dis = std::numeric_limits<double>::max();
    
    for (unsigned int j = 0; j < predicted.size(); ++j){    
      double x2 = predicted[j].x;
      double y2 = predicted[j].y;
      double distance = dist(x1, y1, x2, y2);
      if (distance < min_dis){
        min_dis = distance;
        observations[i].id = predicted[j].id;
      }
      
    }
  }

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
  double std_landmark_x = std_landmark[0];
  double std_landmark_y = std_landmark[1];
  
  
  
  for (unsigned int i = 0; i < particles.size(); ++i){
    // Create a vector containing all of the observation map  global coordinate and landmark map global coordinate
    vector<LandmarkObs> observations_map;
    vector<LandmarkObs> landmark_map;
    
    // particle map coordinate and heading
    double x_part = particles[i].x;
    double y_part = particles[i].y;
    double theta = particles[i].theta;
    
    // get the landmark within sensor range
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j){
      double landmark_x = (double)map_landmarks.landmark_list[j].x_f;
      double landmark_y = (double)map_landmarks.landmark_list[j].y_f;
      
      if(fabs(landmark_x - x_part) <= sensor_range && fabs(landmark_y - y_part) <= sensor_range){
        LandmarkObs landmark;
        landmark.x = landmark_x;
        landmark.y = landmark_y;
        landmark.id = map_landmarks.landmark_list[j].id_i;
        landmark_map.push_back(landmark);
      }
    }
    
    // convert the observation coordinate to map coordinate and push_back to vector "observations_map"
    for (unsigned int j = 0; j < observations.size(); ++j){
      LandmarkObs obs_map;
      double x_obs = observations[j].x;
      double y_obs = observations[j].y;
      obs_map.x = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
      obs_map.y = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
      obs_map.id = observations[j].id;
      observations_map.push_back(obs_map);
    }
    
    dataAssociation(landmark_map, observations_map);
    
    double weight = 1.0;
    for (unsigned int k = 0; k < observations_map.size(); ++k){
      int id = observations_map[k].id;
      double nearst_landmark_x, nearst_landmark_y;
      for (unsigned int i = 0; i < landmark_map.size(); ++i){
        if(landmark_map[i].id == id){
          nearst_landmark_x = landmark_map[i].x;
          nearst_landmark_y = landmark_map[i].y;
          break;
        }
      }
      double obs_x = observations_map[k].x;
      double obs_y = observations_map[k].y;
      weight *= multiv_prob(std_landmark_x, std_landmark_y, obs_x, obs_y,
                   nearst_landmark_x, nearst_landmark_y);
    }
    
    particles[i].weight = weight;
 
  }
  

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> particles_2;
  vector<double> weight;
  std::default_random_engine gen;
  
  // Save weight into a vector
  for (int i = 0; i < num_particles; i++) {
    weight.push_back(particles[i].weight);
  }
  // Create a discrete distribution
  std::discrete_distribution<> dist(weight.begin(), weight.end());
  
  // "dist(gen)" will pick a number from 0 to num_particles with weight.
  for(int i = 0; i < num_particles; i++) {
    particles_2.push_back(particles[dist(gen)]);
  }
  particles=particles_2;

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

