#include<file.h>
int file_size() {
    fstream file;
    file.open("Attributes", ios::in);        

    if (!file)     
    {
        cerr << "Can't open file!\n";
        exit(1);     
    }
    string a;
    int id = 0;
    while (file >> a) {     
        double s = atof(a.c_str());
        id++;
    }

    file.close();
    return id;
}//FUNCTION: int file_size()

int read_file(int size,double pos[][3],double vel[][3]) {
    fstream file;
    file.open("Attributes", ios::in);  

    if (!file)    
    {
        cerr << "Can't open file!\n";
        exit(1);    
    }
    string a;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            file >> a;
            pos[i][j] = atof(a.c_str());
        }
        for (int j = 0; j < 3; j++) {
            file >> a;
            vel[i][j] = atof(a.c_str());
        }
    }
    file.close();
}//FUNCTION: int read_file