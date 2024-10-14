#include "lbm.hpp"

void write_restart(int &step, LBM *lbm);
LBM read_restart(const std::string& filename);
size_t getFileSize(int fd);