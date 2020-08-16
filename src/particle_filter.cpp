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
using std::uniform_real_distribution;
using std::uniform_int_distribution;

std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles    = 20;  // TODO: Set the number of particles
  
    
  //create normal distirbutions for positions and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  for (int i=0 ; i<num_particles; i++) {
    
    //initialize new particle
    Particle pi;
    
    //assign values to particle    
    pi.id = i;
    pi.x = dist_x(gen);
    pi.y = dist_y(gen);
    pi.theta = dist_theta(gen);
    pi.weight = 1;
    
    //add to vector of particles
    particles.push_back(pi);
      
  }
  
  is_initialized=true; 

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
  
  //create normal distirbutions for positions and theta
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  for (unsigned int i=0; i<particles.size(); i++){
    
    
    if(fabs(yaw_rate) < 0.0001){
         
      particles[i].x += velocity * cos(particles[i].theta) * delta_t;
      particles[i].y += velocity * sin(particles[i].theta) * delta_t;   
    }
    
    else{
      //move each particle according to bicicle model
      particles[i].x = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y = particles[i].y + velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t));
      particles[i].theta = particles[i].theta + yaw_rate*delta_t;
    }
      
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
    
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
  
  
  //for all observations
  for (unsigned int k=0 ; k< observations.size();k++) {
    
    double near_distance = std::numeric_limits<double>::max();
    int nearest = -1;
    
    for(unsigned int l=0 ; l<predicted.size(); l++) {
      
      double distance = dist(observations[k].x, observations[k].y, predicted[l].x, predicted[l].y); 
      
      if (distance < near_distance){
        
      nearest = predicted[l].id;
      near_distance = distance;
      }      
      
    }    
    
    observations[k].id =nearest;
    
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
  
  
  /*
  * predict measurement to map landmarks
  */
  
  //for all particles
  for (unsigned int i=0 ; i<particles.size(); i++) {
    
    particles[i].weight=1;
    
    vector<LandmarkObs> landmarks_local;
    
    //for all landmarks check if is in sensor range and transform to local coordinates
    for (unsigned int j=0; j<map_landmarks.landmark_list.size(); j++){
      
      
      double x_l_global = map_landmarks.landmark_list[j].x_f;
      double y_l_global = map_landmarks.landmark_list[j].y_f;
      
      //distance between particle and landmakr
      double dist_part_land = dist(particles[i].x, particles[i].y , x_l_global, y_l_global);
      
      // if landmark is in sensor range
      if ( dist_part_land <sensor_range){
        
        //tanslation of landmark to paricle coordiantes
        double x_l_trans = x_l_global-particles[i].x;
        double y_l_trans = y_l_global-particles[i].y;
        
        //rotation of landmark to particle coordinates
        double x_l_local = x_l_trans*cos(particles[i].theta) + y_l_trans*sin(particles[i].theta);
        double y_l_local = -x_l_trans*sin(particles[i].theta) + y_l_trans*cos(particles[i].theta);
        
        LandmarkObs landmark_local;
        
        landmark_local.id = map_landmarks.landmark_list[j].id_i;
        landmark_local.x = x_l_local;
        landmark_local.y = y_l_local;  
          
        landmarks_local.push_back(landmark_local);
        
      }      
      
    }
    
    vector<LandmarkObs> observations_assigned = observations; 
    dataAssociation(landmarks_local, observations_assigned);
    
    /*
    * update weights of each particle using milti-variate Gaussian pdf 
    */
    
    float x_obs;
    float y_obs;
    float mu_x;
    float mu_y;
    float sigma_x = std_landmark[0];
    float sigma_y = std_landmark[1];
    
    double weight;
    double exponent;
    
    // pre-calculate normalization term
    double gauss_norm = 1 / (2 * M_PI * sigma_x * sigma_y);    

    
    for (unsigned int ii=0; ii<observations_assigned.size(); ii++){      
      
      x_obs = observations_assigned[ii].x;
      y_obs = observations_assigned[ii].y;      
      
      //for all landmarks in range compare landmark id with assigned id
      for (unsigned int jj=0; jj<landmarks_local.size(); jj++){
        
        if (observations_assigned[ii].id == landmarks_local[jj].id){
          
         mu_x= landmarks_local[jj].x;
         mu_y= landmarks_local[jj].y;
          
        }        
      }
      
      // calculate exponent
      exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sigma_x, 2))) + (pow(y_obs - mu_y, 2) / (2 * pow(sigma_y, 2)));      
      
      // calculate weight using normalization terms and exponent
      weight = gauss_norm * exp(-exponent);
      particles[i].weight*= weight;      
    }      
  
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  int num_particless = particles.size();
  
  // find max_weight of all particles
  double max_weight = 0.0;
  
  for (int i=0 ; i < num_particless; i++){
    if (particles[i].weight > max_weight){
      max_weight =  particles[i].weight;
    }
  }  
  
  double beta = 0;
  vector<Particle> particles_resampled;
  
  uniform_real_distribution<double> dist_double_weight(0, 2*max_weight);
  uniform_int_distribution<int> dist_index(0, num_particles - 1);
  
  //set random initial index
  int index = dist_index(gen);
  
  
  
  //repeat num_particles times
  for (int j=0; j<num_particless; j++){
    
    beta += dist_double_weight(gen);
    
    while (particles[index].weight < beta){      
      beta -= particles[index].weight;
      index = (index+1) % num_particles;      
    }
    
    particles_resampled.push_back(particles[index]);      
    
  }
  
 // std::cout << particles[0].theta << std::endl; 
  
  
  particles = particles_resampled;  

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