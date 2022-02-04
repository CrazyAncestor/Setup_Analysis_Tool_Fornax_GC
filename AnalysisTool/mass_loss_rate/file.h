#ifndef FILE_H
#define FILE_H
#include<fstream>
#include<string.h>
#include<iostream>
using namespace std;
int file_size();
int read_file(int size,double pos[][3],double vel[][3]);
#endif // #ifndef FILE_H