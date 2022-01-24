#ifndef READ_DATA_H
#define READ_DATA_H
#include "TFile.h"
#include "TTree.h"
#include "AngDisXiXi.hh"
#include <iostream>
using namespace std;

int readData(const char *filename, AngDisXiXi *ngdis, double **para, int flag, const char *type, const int index, const int MM);

#endif // READ_DATA_H
