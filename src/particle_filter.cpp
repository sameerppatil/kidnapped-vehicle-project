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

default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 70;

    double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

    // TODO: Set standard deviations for x, y, and theta
     std_x = std[0];
     std_y = std[1];
     std_theta = std[2];

     particles.resize(num_particles);
     weights.resize(num_particles);

    // This line creates a normal (Gaussian) distribution for x
    normal_distribution<double> dist_x(x, std_x);

    // TODO: Create normal distributions for y and theta
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);


    for (int i = 0; i < 3; ++i) {
        double sample_x, sample_y, sample_theta;

        // TODO: Sample  and from these normal distrubtions like this:
        //   sample_x = dist_x(gen);
        //   where "gen" is the random engine initialized earlier.

            particles[i].id = i;
            particles[i].x = dist_x(gen);
            particles[i].y = dist_y(gen);
            particles[i].theta = dist_theta(gen);
            particles[i].weight = 1.0;

         // Print your samples to the terminal.
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

  // TODO: Set standard deviations for x, y, and theta
   std_x = std_pos[0];
   std_y = std_pos[1];
   std_theta = std_pos[2];

  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(0, std_x);

  // TODO: Create normal distributions for y and theta
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);

  for (int i = 0; i < num_particles; i++)
  {
    if (fabs(yaw_rate) < 0.0001)
    {
      particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
      particles[i].y = particles[i].y + velocity*delta_t*cos(particles[i].theta);
    }
    else
    {
      particles[i].x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate
        * delta_t) - sin(particles[i].theta));
      particles[i].y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) -
        cos(particles[i].theta + yaw_rate * delta_t));

      particles[i].theta = particles[i].theta + yaw_rate * delta_t;
    }

    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
  return;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
  double dist1, dist2, min_distance;
  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];
  double mu_x, mu_y;

  for (int i = 0; i < num_particles; i++)
  {
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;

    particles[i].weight = 1;

    vector<LandmarkObs> nearby_landmarks;
    for (int m = 0; m < map_landmarks.landmark_list.size(); ++m)
    {
      /* code */
      dist1 = dist(p_x, p_y, map_landmarks.landmark_list[m].x_f,
        map_landmarks.landmark_list[m].y_f);
      if (dist1 <=  sensor_range)
      {
        nearby_landmarks.push_back(LandmarkObs{map_landmarks.landmark_list[m].id_i,
          map_landmarks.landmark_list[m].x_f, map_landmarks.landmark_list[m].y_f});
      }
    }
    mu_x = 0;
    mu_y = 0;

    for (LandmarkObs &obs : observations)
    {
      /* code */
      double x_map = cos(p_theta) * obs.x - sin(p_theta) * obs.y + p_x;
      double y_map = sin(p_theta) * obs.x + cos(p_theta) * obs.y + p_y;

      min_distance = dist(x_map, y_map, nearby_landmarks[0].x, nearby_landmarks[0].y);

      for(LandmarkObs &local_landmarks : nearby_landmarks)
      {
        dist2 = dist(x_map, y_map, local_landmarks.x, local_landmarks.y);
        if (dist2 <= min_distance)
        {
          min_distance = dist2;
          mu_x = local_landmarks.x;
          mu_y = local_landmarks.y;
        }
      }
      double exponent = (pow(x_map - mu_x, 2))/(2 * pow(sig_x, 2)) +
        (pow(y_map - mu_y, 2))/(2 * pow(sig_y, 2));

      double partial_weight = (1/(2 * M_PI * sig_x * sig_y))
        * exp(-exponent);

      particles[i].weight *= partial_weight;
      weights[i] = particles[i].weight;
    }
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  vector<Particle> resampled_particle_list(num_particles);

  discrete_distribution<int> distribution(weights.begin(), weights.end());
  for (Particle &sample : resampled_particle_list)
  {
    sample = particles[distribution(gen)];
  }

  particles = resampled_particle_list;
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
