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

using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	num_particles = 100;

	for(int i = 0; i < num_particles; i++){

		weights.push_back(1.0);
		Particle part;
		part.id = i;
		part.x = dist_x(gen);
		part.y = dist_y(gen);
		part.theta = dist_theta(gen);
		part.weight = 1.0 / num_particles;
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
		double xf, yf, thetaf;
		if(abs(yaw_rate)<0.0001){
			xf = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			yf = particles[i].y + velocity * delta_t * sin(particles[i].theta);
		}
		else{
			xf = particles[i].x + velocity*(sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta))/yaw_rate;
			yf = particles[i].y + velocity*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t))/yaw_rate;			
		}
		thetaf = particles[i].theta + yaw_rate * delta_t;

		normal_distribution<double> dist_x(xf, std_pos[0]);
		normal_distribution<double> dist_y(yf, std_pos[1]);
		normal_distribution<double> dist_tehta(thetaf, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_tehta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for(unsigned int j = 0; j < observations.size(); j++){

		double min_distance = numeric_limits<double>::max();

		for(unsigned int i = 0; i < predicted.size(); i++){

			double dx = predicted[i].x - observations[j].x;
			double dy = predicted[i].y - observations[j].y;
			double distance = sqrt(dx*dx + dy*dy);

			if(min_distance > distance){
				min_distance = distance;
				observations[j].id = predicted[i].id;
			}
		}
		//cout<<"min distance: "<<min_distance<<endl;
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

	for(int j = 0; j < num_particles; j++){

		vector<LandmarkObs> trans_observations, predictions;

		for(unsigned i = 0; i < map_landmarks.landmark_list.size(); i++){

			if((abs(map_landmarks.landmark_list[i].x_f - particles[j].x) < sensor_range) && 
			   (abs(map_landmarks.landmark_list[i].y_f - particles[j].y) < sensor_range)){

				LandmarkObs tempt;
				tempt.x = map_landmarks.landmark_list[i].x_f;
				tempt.y = map_landmarks.landmark_list[i].y_f;
				tempt.id = map_landmarks.landmark_list[i].id_i;
				predictions.push_back(tempt);
				}
			}
		
		for(unsigned int i = 0; i < observations.size(); i++){

			LandmarkObs tempt;
			tempt.x = particles[j].x + observations[i].x * cos(particles[j].theta) - observations[i].y * sin(particles[j].theta);
			tempt.y = particles[j].y + observations[i].x * sin(particles[j].theta) + observations[i].y * cos(particles[j].theta);
			tempt.id = -1;
			trans_observations.push_back(tempt);
		}

		dataAssociation(predictions, trans_observations);

		weights[j] = 1.0;
		for(unsigned int i = 0; i < trans_observations.size(); i++){

			double predict_x, predict_y;

			for(vector<LandmarkObs>::iterator iter = predictions.begin(); iter != predictions.end(); iter++){
				if(trans_observations[i].id == (*iter).id){
					predict_x = (*iter).x;
					predict_y = (*iter).y;
				}
			}

			double gauss_norm = 1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]);
			double dx = trans_observations[i].x - predict_x;
			double dy = trans_observations[i].y - predict_y;
			//cout<<"dx: "<<dx<<", dy: "<<dy<<endl;
			double exponent = (dx * dx /(std_landmark[0] * std_landmark[0] * 2) + dy * dy / (std_landmark[1] * std_landmark[1] * 2));
			weights[j] *= (gauss_norm * exp(-exponent));
			if(weights[j] < 0.000001){
				weights[j] = 0.000001;
			}
			//cout<<i<<"th"<<" exponent: "<<exponent<<endl;
		}			
	}

	double sum = 0.0;
	for(vector<double>::iterator iter = weights.begin(); iter != weights.end(); iter++){
		sum += *iter;
	}
	//cout<<"summary: "<<sum<<endl;

	for(unsigned int i = 0; i < weights.size(); i++){
		
		weights[i] /= sum;
		particles[i].weight = weights[i];
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> tempt;

	double belta = 0.0;
	discrete_distribution<> d(weights.begin(), weights.end());
	//cout<<d(gen)<<endl;
	/*int index = (int)(rand() % num_particles);
	vector<double>::iterator p = max_element(weights.begin(), weights.end());
	discrete_distribution<int> wei(0,1000);
	uniform_real_distribution<double> add_belta(0.0, 2 * (*p));*/

	int index = d(gen);

	for(int i = 0; i < num_particles; i++){

		//belta += add_belta(gen);
		belta += weights[d(gen)] * 2;

		while(weights[index] < belta){

			belta -= weights[index];
			index = (index + 1) % num_particles;
		}
		tempt.push_back(particles[index]);
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
