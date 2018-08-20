/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "map.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 1000;

	for(int i = 0; i < num_particles; i++){
		weights.push_back(1.0);
		Particle part;
		part.id = i;
		part.x = random.gauss(x, std[0]);
		part.y = random.gauss(y, std[1]);
		part.theta = random.gauss(theta, std[2]);
		part.weight = 1.0;
		particles.push_back(part);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	for(unsigned int i = 0; i < particles.size(); i++){
		double xf = particles[i].x + velocity*(sin(particels[i].theta + yaw_rate * delta_t) - sin(particles[i].theta))/yaw_rate;
		double yf = particles[i].y + velocity*(cos(particles[i].theta) - cos(particels[i].theta + yaw_rate * delta_t))/yaw_rate;
		double thetaf = partiicles[i].theta + yaw_rate * delta_t;
		particles[i].x = random.gauss(xf, std_pos[0]);
		particles[i].y = random.gauss(yf, std_pos[1]);
		particles[i].theta = random.gauss(thetaf, std_pos[2]);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for(unsigned int j = 0; i < observations.size(); j++){

		double min_distance = 999999999;

		for(unsigned int i = 0; i < predicted.size(); i++){

			double dx = predicted.x - *observations[j].x;
			double dy = predicted.y - *observations[j].y;
			double distance = sqrt(dx*dx + dy*dy);

			if(min_distance > distance){
				min_distance = distance;
				*observations[j].id = predicted.id;
			}
		}
	}
	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	for(unsigned int j = 0; j < particles.size(); j++){

		vector<LandmarkObs> trans_observations;
		
		for(unsigned int i = 0; i < *observations.size(); i++){
			LandmarkObs tempt;
			tempt.x = particles[j].x + *observations[i].x * cos(particles[j].theta) - *observations[i].y * sin(particles[j].theta);
			tempt.y = particles[j].y + *observations[i].x * sin(particles[j].theta) - *observations[i].y * cos(particles[j].theta);
			if((abs(tempt.x) < sensor_range) && (abs(tempt.y) < sensor_range)){
				trans_observations.push_back(tempt);
			}
		}
		dataAssociation(*map_landmarks.landmark_list, &trans_observations);

		weights[j] = 1.0;
		for(unsigned int i = 0; i < trans_observations.size(); i++){
			double gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
			double dx = trans_observations[i].x - map_landmarks.landmark_list[trans_observations.id].x;
			double dy = trans_observations[i].y - map_landmarks.landmark_list[trans_observations.id].y;
			double exponent = (dx * dx /(std_landmark[0] * std_landmark[0] * 2) + dy * dy / (std_landmark[1] * std_landmark[1] * 2));
			weights[j] *= (gauss_norm * exp(-exponent));
		}
	}

	sum = accumulate(weights.begin(), weights.end());
	for(unsigned int i = 0; i < weights.size(); i++){
		particles[i].weight = weights[j] / sum;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;

	vector<Particle> tempt;
	double belta = 0.0;
	int index = (int)(rand() % num_particles);
	vector<double>::iterator p = max_element(weights);
	for(unsigned int i = 0; i < particles.size(); i++){
		discrete_distribution<int> wei(0,1000);
		belta += (double)wei(gen)/1000.0 * 2 * (*p);
		while(weights[index] < belta){
			belta -= weights[index];
			tempt.push_back(particles[index]);
		}
	}
	particles = tempt;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
