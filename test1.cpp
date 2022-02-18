/*
My initial goal is to recreate this simple test script in C++:

  import galsim

  print("Creating object")
  obj = galsim.Gaussian(flux=10000, sigma=0.3)
  print("Creating image")
  im1 = galsim.ImageD(1024, 1024, scale=0.3)
  print("Creating RNG")
  rng1 = galsim.BaseDeviate(5678)
  print("Creating sensor")
  silicon = galsim.SiliconSensor(rng=rng1, diffusion_factor=0.0)
  print("Drawing image")
  obj.drawImage(im1, method='phot', poisson_flux=False, sensor=silicon, rng=rng1)
  print("Done!")
*/

#include "galsim/GSParams.h"
#include "galsim/SBGaussian.h"
#include "galsim/Image.h"
#include "galsim/Random.h"
#include "galsim/Silicon.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main(int argc, char *argv[])
{
    std::string sensor_data_filename = "/global/homes/j/jamesp/gpu/GalSim/share/sensors/lsst_itl_50_8.dat";
    std::string abs_length_filename = "/global/homes/j/jamesp/gpu/GalSim/share/sensors/abs_length.dat";
    
    // Create object
    std::cout << "Creating object" << std::endl;    
    galsim::GSParams gsparams;
    galsim::SBGaussian gaussian(0.3, 10000.0, gsparams);

    // Create image
    std::cout << "Creating image" << std::endl;
    galsim::ImageAlloc<double> im1(1024, 1024);
    // scale??
    
    // Create RNG
    std::cout << "Creating RNG" << std::endl;
    galsim::BaseDeviate rng1(5678);

    // Create sensor
    /*
        Silicon(int numVertices, double numElec, int nx, int ny, int qDist, double nrecalc,
                double diffStep, double pixelSize, double sensorThickness, double* vertex_data,
                const Table& tr_radial_table, Position<double> treeRingCenter,
                const Table& abs_length_table, bool transpose);
    */
    std::cout << "Creating sensor" << std::endl;

    // need to read from vertex data file
    std::vector<double> vertex_data_v;
    std::ifstream vertex_data_file(sensor_data_filename.c_str());
    if (vertex_data_file.fail()) {
	std::cerr << "Error opening data file " << sensor_data_filename << std::endl;
	return 1;
    }
    // skip header line
    std::string header;
    std::getline(vertex_data_file, header);
    while (!vertex_data_file.eof()) {
	double val;
	vertex_data_file >> val;
	if (!vertex_data_file.eof()) vertex_data_v.push_back(val);
    }
    vertex_data_file.close();
    std::cout << "Read " << vertex_data_v.size() << " values from sensor data file" << std::endl;
    
    double* vertex_data = new double[vertex_data_v.size()];
    for (int i = 0; i < vertex_data_v.size(); i++) {
	vertex_data[i] = vertex_data_v[i];
    }

    double tr_radial_table_args[2] = { 0.0, 1.0 };
    double tr_radial_table_vals[2] = { 0.0, 0.0 };
    int tr_radial_table_N = 2;
    galsim::Table tr_radial_table(tr_radial_table_args, tr_radial_table_vals, tr_radial_table_N,
				  galsim::Table::interpolant::linear);
    
    galsim::Position<double> treeRingCenter(0.0, 0.0);

    // need to generate this table
    std::vector<double> abs_length_args_v, abs_length_vals_v;
    std::ifstream abs_length_file(abs_length_filename.c_str());
    if (abs_length_file.fail()) {
	std::cerr << "Error opening data file " << abs_length_filename << std::endl;
	return 1;
    }
    // skip header line
    std::getline(abs_length_file, header);
    while (!abs_length_file.eof()) {
	double arg, val;
	abs_length_file >> arg;
	abs_length_file >> val;
	if (!abs_length_file.eof()) {
	    abs_length_args_v.push_back(arg);
	    abs_length_vals_v.push_back(val);
	}
    }
    abs_length_file.close();
    std::cout << "Read " << abs_length_args_v.size() << " entries from abs_length table" << std::endl;
    
    double* abs_length_table_args = new double[abs_length_args_v.size()];
    double* abs_length_table_vals = new double[abs_length_vals_v.size()];
    for (int i = 0; i < abs_length_args_v.size(); i++) {
	abs_length_table_args[i] = abs_length_args_v[i];
	abs_length_table_vals[i] = abs_length_vals_v[i];
    }
    int abs_length_table_N = abs_length_args_v.size();
    galsim::Table abs_length_table(abs_length_table_args, abs_length_table_vals, abs_length_table_N,
				   galsim::Table::interpolant::linear);

    std::cout << "Creating silicon sensor" << std::endl;
    //galsim::Silicon silicon(8, 100000.0, 9, 9, 3, 10000, 0.0, 10.0, 100.0, vertex_data,
    //  tr_radial_table, treeRingCenter, abs_length_table, false);
    galsim::Silicon silicon(8, 100000.0, 9, 9, 3, 0.0, 10.0, 100.0, vertex_data,
      tr_radial_table, treeRingCenter, abs_length_table, false);
    
    delete[] abs_length_table_vals;
    delete[] abs_length_table_args;
    delete[] vertex_data;

    // Draw image
    std::cout << "Shooting photons" << std::endl;
    // Get photons first
    galsim::PhotonArray photons(1000000);
    gaussian.shoot(photons, rng1);

    // Now fire them at sensor
    galsim::Position<int> origCentre(0, 0);
    silicon.initialize(im1.view(), origCentre);

    silicon.accumulateGPU(photons, 0, 1000000, rng1, im1.view());  // GPU version
    //silicon.accumulate(photons, 0, 1000000, rng1, im1.view());   // CPU version

    // Add delta image to actual image
    silicon.addDelta(im1.view());

    // Save image data in raw binary
    std::cout << "Saving result" << std::endl;
    double* imageData = im1.getData();
    int imageXMin = im1.getXMin();
    int imageYMin = im1.getYMin();
    int imageXMax = im1.getXMax();
    int imageYMax = im1.getYMax();
    int imageStep = im1.getStep();
    int imageStride = im1.getStride();

    std::ofstream out("output.bin", std::ios::binary | std::ios::out);
    if (out.fail()) {
	std::cerr << "Unable to open output.bin for writing" << std::endl;
	return 1;
    }
    for (int y = 0; y < (imageYMax - imageYMin); y++) {
	for (int x = 0; x < (imageXMax - imageXMin); x++) {
	    int idx = (x * imageStep) + (y * imageStride);
	    out.write((const char *)&imageData[idx], sizeof(double));
	}
    }
    out.close();
    
    std::cout << "Done" << std::endl;
    
    return 0;
}
