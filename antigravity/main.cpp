//
//  main.cpp
//  antigravity
//
//  Created by Tom Andersen on 2019-08-21.
//  Copyright © 2019 Tom Andersen. All rights reserved.
//

//#include <iostream>
//
//int main(int argc, const char * argv[]) {
//    // insert code here...
//    std::cout << "Hello, Worldjjjfff!\n";
//    return 0;
//}
//
#include <cstring>

#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/program_options.hpp>
#include <exception>

namespace po = boost::program_options;
// Make struct/class for settings
class Settings
{
    public:
        double lattice;
    
        double xi;
        double yi;
        double zi;

        double ahbar;
    
        long nX;
        long nY;
        long nZ;
};
// make functions to evaluate Kerr metric at given x, y, z

// util to read an item
static double readItem(po::variables_map &vm, const char* itemName, double defaultVal)
{
    using std::cout;
    if (vm.count(itemName)) {
        cout << itemName << " set to " << vm[itemName].as<double>() << " \n";
        return vm[itemName].as<double>();
    }
    cout << itemName << " is default " << defaultVal << " .\n";
    return defaultVal;
}
// read commands
static int inputVariables(int ac, char** av, Settings& outSettings)
{
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::exception;
    try {
        
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "This program calculates the metric underneath an array of spinning atoms, assuming that each atom is some kind of spinning Kerr like system")
        ("lattice", po::value<double>(), "atomic lattice spacing, nanometers")
        ("xi", po::value<double>(), "x distance from crystal face, nanometers")
        ("yi", po::value<double>(), "y distance from crystal face, nanometers")
        ("zi", po::value<double>(), "y distance from crystal face, nanometers")
        ("ahbar", po::value<double>(), "spin per atom, units of hbar - eg 0.5")
        ("nX", po::value<double>(), "number of atoms in sample, X eg 1000")
        ("nY", po::value<double>(), "number of atoms in sample, Y eg 1000")
        ("nZ", po::value<double>(), "number of atoms in sample, Z eg 20");

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);
        
        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }
        
        outSettings.lattice = readItem(vm, "lattice", 2.0);
        outSettings.xi = readItem(vm, "xi", 0.0);
        outSettings.yi = readItem(vm, "yi", 0.0);
        outSettings.zi = readItem(vm, "zi", 40.0);
        outSettings.ahbar = readItem(vm, "ahbar", 0.5);
        outSettings.nX = readItem(vm, "nX", 100);
        outSettings.nY = readItem(vm, "nY", 100);
        outSettings.nZ = readItem(vm, "nZ", 20);
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        return 1;
    }
    return 0;
}

int main(int ac, char** av)
{
    using namespace boost::lambda;
    Settings settings;
    inputVariables(ac, av, settings);
    
    
    
//    typedef std::istream_iterator<double> in;
//
//    std::for_each(in(std::cin), in(),
//        std::cout << "*3\n----\n"  << (_1 * 3) << "\n\n"
//    );
    
}
