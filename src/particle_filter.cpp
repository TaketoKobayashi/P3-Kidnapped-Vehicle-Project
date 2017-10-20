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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> yaw(theta, std[2]);
	
	num_particles = 36;
	
	for(int i = 0; i < num_particles; i++) {
		Particle pf;
		pf.id = i;
		pf.x = dist_x(gen);
		pf.y = dist_y(gen);
		pf.theta = yaw(gen);
		pf.weight = 1.0;
		particles.push_back(pf);
		weights.push_back(pf.weight);
	}
	
	is_initialized = true;

	return;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> yaw(0, std_pos[2]);
	
	double theta = yaw_rate * delta_t;

	for(int i = 0; i < num_particles; i++) {

		if( fabs(yaw_rate) > 0.001) {
			double new_theta = particles[i].theta + theta;
			particles[i].x = particles[i].x + velocity / yaw_rate * (sin(new_theta) - sin(particles[i].theta));
		
		    particles[i].y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(new_theta));
		
		    particles[i].theta = new_theta;
		}
		else {
			particles[i].x = particles[i].x + velocity * cos(particles[i].theta) * delta_t;
			particles[i].y = particles[i].y + velocity * sin(particles[i].theta) * delta_t;
		}
		
		// TODO: add noise or not? in init function, it is already add noise.
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += yaw(gen);
		
		//particles[i].id = i;

	}
	
	return;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	vector<LandmarkObs> particular_landmark;
	
	LandmarkObs landmark;
	
	for(int i = 0; i < observations.size(); i++ ) {
		double x0 = observations[i].x;
		double y0 = observations[i].y;
		std::vector< std::pair<double, int> > dist_id;
		
		for(int j = 0; j < predicted.size(); j++ ) {
			double x1 = predicted[j].x;
			double y1 = predicted[j].y;
			
			// calculate min distance between the two points
			double distance = dist(x0, y0, x1, y1);  //sqrt((x0-x1) * (x0-x1) + (y0-y1) * (y0-y1));
			dist_id.push_back(make_pair(distance, j));
		}
		
		//get the min distance and matched id
		double min_dist = dist_id[0].first;
		int min_id = dist_id[0].second;
		
		for(int k = 0; k < dist_id.size(); k++) {
			if(min_dist > dist_id[k].first) {
				min_dist = dist_id[k].first;
				min_id = dist_id[k].second;
			}
		}
		
		landmark.id = predicted[min_id].id;
		landmark.x = x0;
		landmark.y = y0;
		
		particular_landmark.push_back(landmark);
		
	}
	
	observations.assign(particular_landmark.begin(), particular_landmark.end() );
	
	return;
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
	
	weights.clear();
	double x_w = 2.0 * std_landmark[0] * std_landmark[0];
	double y_w = 2.0 * std_landmark[1] * std_landmark[1];
	double xy_w = 2.0 * M_PI * std_landmark[0] * std_landmark[1];
	
	for(int i = 0; i < num_particles; i++) {
		vector<LandmarkObs> filter_obs_landmark;
		
		//filter the map landmark which is out of range of sensor_range
		for(int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			double distance = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
			if(distance < sensor_range) {
				LandmarkObs landobs;
				landobs.id = map_landmarks.landmark_list[j].id_i;
				landobs.x = map_landmarks.landmark_list[j].x_f;
				landobs.y = map_landmarks.landmark_list[j].y_f;
				filter_obs_landmark.push_back(landobs);
			}
		}
		
		//transform each observation marker from the vehicle's coordinates to the map's coordinates, with respect to our particle
		vector<LandmarkObs> landmark_map;
		LandmarkObs landmark_new;
		for(int k = 0; k < observations.size(); k++) {
			LandmarkObs obs = observations[k];
			landmark_new = TransformToMap(obs, particles[i]);
			landmark_map.push_back(landmark_new);
		}

		dataAssociation(filter_obs_landmark, landmark_map);
		
		//set initial value
		particles[i].weight = 1.0;
		
		for(int m = 0; m < landmark_map.size(); m++) {
			Map::single_landmark_s landmark = map_landmarks.landmark_list.at(landmark_map[m].id);
			double x_w1 = (landmark_map[m].x - landmark.x_f) * (landmark_map[m].x - landmark.x_f) / x_w;
			double y_w1 = (landmark_map[m].y - landmark.y_f) * (landmark_map[m].y - landmark.y_f) / y_w;
			double weigh = 1.0 / xy_w * exp(-1.0 * (x_w1 + y_w1));
			
			//double ParticleFilter::CalGaussian(LandmarkObs predicted, LandmarkObs sensor, double std[]);
			
			particles[i].weight *= weigh;
		}
		
		weights.push_back(particles[i].weight);
	}
	
	return;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::default_random_engine gen;
	std::discrete_distribution<int> dist(weights.begin(), weights.end());
	int index = int(dist(gen) * num_particles);
	
	std::vector<Particle> p_w;
	double beta = 0.0;
	vector<double>::iterator biggest = std::max_element(std::begin(weights), std::end(weights));
	double mw = *biggest;
	for(int i = 0; i < num_particles; i++) {
		beta += dist(gen) * 2.0 * mw;
		while(beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		p_w.push_back(particles[index]);
	}
	particles.assign(p_w.begin(), p_w.end() );
	return;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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

double ParticleFilter::CalGaussian(LandmarkObs predicted, LandmarkObs sensor, double std[]) {
	double pow_x = (sensor.x - predicted.x) / std[0];
	pow_x = pow_x * pow_x;
	double pow_y = (sensor.y - predicted.y) / std[1];
	pow_y = pow_y * pow_y;
	
	return exp(-1.0 * (pow_x + pow_y) / 2) / (2.0 * M_PI * std[0] * std[1]);
}

LandmarkObs ParticleFilter::TransformToMap(LandmarkObs& landmark, Particle& particle) {
	LandmarkObs landmark_map;
	double theta = particle.theta;
	landmark_map.x = particle.x + landmark.x * cos(theta) - landmark.y * sin(theta);
	landmark_map.y = particle.y + landmark.x * sin(theta) + landmark.y * cos(theta);
	return landmark_map;
}
