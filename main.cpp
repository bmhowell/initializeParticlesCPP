#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <chrono>
#include <random>



// particle parameters
const double dparticleL{2 * (3.985714285714286e-5)};            // diameter of large particles
const double dparticleS{dparticleL / 10};                       // diameter of small particles
const double dparticleM{(6*dparticleL + 4*dparticleS) / 10};      // median diameter
const double dNozzle{8.37e-4};                                  // diamter of nozzle
const double zHeight{6.2e-4};                                   // zHeight of nozzle
const int scale = 1;                                            // scaling factor
const double xmin{(dparticleL / 2) * scale};                    // xmin for location
const double xmax{(dNozzle - xmin) * scale};                    // xmax for location
const double ymin{(zHeight / 10) * scale};                      // ymin for location
const double ymax{2 * dNozzle * scale};                         // ymax for location

const double areaNozzle{(xmax - xmin) * (ymax - ymin)};
const int nParticles{1000};
double randPosX;                                                // new random position
double randPosY;                                                // new random position
double randDiam;                                                // new random diameter
double rstate[nParticles][3];                                   // state matrix [x, y, diameter]

double packingDensity;                                          // total packing density
double sampleD;                                                 // random diameter from normal distribution

// data outputs
std::ofstream printParticles;
std::ofstream printDiameters;

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "xmin = " << xmin << " xmax = " << xmax << " ymin = " << ymin << " ymax = " << ymax << std::endl;

//    srand(time(NULL));      // seed for random number
    std::random_device rd{};
    std::mt19937 gen{rd()};

    // initialize random positions in state matrix
    rstate[0][0] = ((double) rand() / RAND_MAX) * (xmax - xmin) + xmin;
    rstate[0][1] = ((double) rand() / RAND_MAX) * (ymax - ymin) + ymin;
    rstate[0][2] = ((double) rand() / RAND_MAX) * (dparticleL - dparticleS) + dparticleS;

    // start populating canvas with random, non-overlapping particles
    int particle{1};
    while (particle < nParticles){
        std::normal_distribution<> d{dparticleM, 1e-5};
        sampleD = d(gen);

//        std::cout << "random diameter: " << sampleD << std::endl;

        bool test = true;
        randDiam = sampleD;
        while (randDiam < 0){
            std::cout << "HEERREEE" << std::endl;
            randDiam = d(gen);          // ensure diameter is positive
            if (randDiam > 0){
                std::cout << "   randDiam now positive: " << randDiam << std::endl;
            }
        }
        randPosX = ((double) rand() / RAND_MAX) * (xmax - xmin - sampleD / 2) + xmin + sampleD / 2; // ensure diameter
        randPosY = ((double) rand() / RAND_MAX) * (ymax - ymin - sampleD / 2) + ymin + sampleD / 2; // is considered

        // Loop through all proceeding particles
        for (int i=0; i<particle; ++i){

            double posX_ = rstate[i][0];       // past particle x Position
            double posY_ = rstate[i][1];       // past particle y Position
            double posD_ = rstate[i][2];       // past particle diameter
            double dist = sqrt(pow(posX_ - randPosX, 2.0) + pow(posY_ - randPosY, 2.0));

            // Check if current random particle overlaps with past particles
            if (dist < (randDiam + posD_) / 2){
                test = false;
                break;
            }
        }

        // If there are no overlaps, accept new random particle into array
        if (test == true){
            rstate[particle][0] = randPosX;
            rstate[particle][1] = randPosY;
            rstate[particle][2] = randDiam;
            std::cout << particle << std::endl;
            particle += 1;
        }
    }


    // Compute packing density
    for (int particle=0; particle<nParticles; ++particle){
        packingDensity += M_PI * pow((rstate[particle][2] / 2), 2);
    }
    packingDensity = packingDensity / areaNozzle;
    std::cout << "packing density: " << packingDensity << std::endl;

    // computational time
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "computational time: " << duration << "s" << std::endl;

    // write to file
    printParticles.open("/Users/bhopro/Desktop/Berkeley/MSOL/Projects/initializeParticlesCPP/particlePosition.dat");
    printParticles << "X" << "," << "Y" << "," << "D" << std::endl;

    for (int particle=0; particle<nParticles; ++particle){
        printParticles << rstate[particle][0] << "," << rstate[particle][1] << "," << rstate[particle][2] << std::endl;
    }
    printParticles.close();

    return 0;
}
