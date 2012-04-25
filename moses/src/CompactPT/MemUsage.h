#ifndef MEMUSAGE_H__
#define MEMUSAGE_H__

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace Moses {

void get_mem_usage(double& vm_usage, double& rss_usage);
double get_vm_usage();
double get_rss_usage();

template <typename OStream>
void print_mem_usage(OStream &out) {
    double vm_usage;
    double rss_usage;
    get_mem_usage(vm_usage, rss_usage);
    
    std::stringstream ss;
    ss << std::fixed;
    ss << std::setprecision(3);
    ss << "Memory Usage:\t" << "VIRT: " << vm_usage << "m\tRSS: " << rss_usage << "m" << std::endl;
    
    out << ss.str();
}

void print_mem_usage();

}
#endif
