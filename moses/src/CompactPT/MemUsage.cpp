#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "MemUsage.h"

namespace Moses {

void get_mem_usage(double& vm_usage, double& rss_usage) {
    std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);
    std::string dummy;
    long vsize;
    long rss;
    
    size_t i = 0;
    while(++i < 23)
	stat_stream >> dummy;
    stat_stream >> vsize;
    stat_stream >> rss;
    
    stat_stream.close();
    
    vm_usage = vsize/(1024.0*1024.0);
    rss_usage = rss * (sysconf(_SC_PAGE_SIZE)/(1024.0*1024.0));
}

double get_vm_usage() {
    double vm_usage;
    double rss_usage;
    get_mem_usage(vm_usage, rss_usage);
    return vm_usage;
}

double get_rss_usage() {
    double vm_usage;
    double rss_usage;
    get_mem_usage(vm_usage, rss_usage);
    return rss_usage;
}

void print_mem_usage() {
    print_mem_usage(std::cerr);
}

}
