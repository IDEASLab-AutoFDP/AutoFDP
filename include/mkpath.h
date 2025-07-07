#ifndef _MKPATH_H_
#define _MKPATH_H_
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <errno.h>

#include <dirent.h>
#include <unistd.h>

int IsFolderExist(const char* path);
namespace light
{
    int mkpath(std::string s,mode_t mode=0755);
}
#endif