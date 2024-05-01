#include <iostream>
#include <fstream>
#include<string>
#include<omp.h>
#include<iomanip>
#include<cmath>

//#define CHUNK 2
double Coefficient(double F_Prev, double F_Cur, double F_Next, double dx, double sigma, double K_Prev, double K_Cur, double K_Next, double dt, double constant)  //K_i
{
  double var = 0.0;
  double temp = 0.0;
  if (constant != 0)
    var = dt / constant;
  else throw "Error";

  temp = F_Prev + K_Prev * var - 2 * F_Cur - 2 * K_Cur * var + F_Next + K_Next * var;
  return sigma / (dx * dx) * temp;
}

void Insert_In_File(double* new_array, size_t count, std::fstream& stream)
{
  if (stream.is_open())
  {
    for (size_t i = 0; i < count; i++)
      stream << new_array[i] << std::setprecision(16)<< " ";
    stream << "\n";
    stream << "\n";
  }
  else
    throw "The file isn't exist";
}

int main()
{
  double tMax = 0.0;
  double sigma = 0.0;
  double xMax = 0.0;
  double deltaX = 0.0;
  double deltaT = 0.0;
  double time = 0.0;
  double pi = 3.141592653589793;
  double t = 0.0;
  size_t i = 1;
  double t1 = 0, t2 = 0, dt = 0;

  std::fstream write;
  write.open("output.txt", std::fstream::out);

  std::string path = "const_initial.txt";
  std::ifstream read;
  read.open(path);


  if (!read.is_open())
  {
    std::cout << "file not found";

  }
  else {
    read >> xMax;
    read >> deltaX;
    read >> tMax;
    read >> deltaT;
    read >> sigma;
    //cout << tMax;
  }
  read.close();

  int count = static_cast<int>((xMax / deltaX) + 1);

  double* temporary = nullptr;
  temporary = new double[count];

  double* current = nullptr;
  current = new double[count];

  double* temp_array = nullptr;
  temp_array = new double[count];

  double* K1 = new double[count];
  double* K2 = new double[count];
  double* K3 = new double[count];
  double* K4 = new double[count];


  for (size_t j = 0; j < count; j++)
  {
    double x = j * deltaX;
    temporary[j] = sin(pi * deltaX * j);
   // std::cout << temporary[j] << std::endl;

  }

  temporary[count - 1] = 0.0;
  //std::cout << "start" << std::endl;


  Insert_In_File(temporary, count, write);

  omp_set_num_threads(8); //задает число потоков для выполнения следующего параллельного региона
  //std::cout << omp_get_num_threads() << std::endl;

  t1 = omp_get_wtime();

  while (t <= tMax)
  {
    t = i * deltaT;


    current[0] = 0.0;
    current[count - 1] = 0.0;

    K1[0] = 0.0;
    K1[count - 1] = 0.0;

    K2[0] = 0.0;
    K2[count - 1] = 0.0;

    K3[0] = 0.0;
    K3[count - 1] = 0.0;

    K4[0] = 0.0;
    K4[count - 1] = 0.0;


//#pragma omp parallel for shared(K1) schedule (static)
    for (int n = 1; n < count - 1; n++) {
      //      if (n == 1) {
      //       std::cout << omp_get_num_threads() << std::endl;
      //        std::cout.flush();
      //      }
      K1[n] = Coefficient(temporary[n - 1], temporary[n], temporary[n + 1], deltaX, 1, 0, 0, 0, 0, 1);
    }


//#pragma omp parallel for shared(K1,K2) schedule (static)

    for (int n = 1; n < count - 1; n++)
      K2[n] = Coefficient(temporary[n - 1], temporary[n], temporary[n + 1], deltaX, 1, K1[n - 1], K1[n], K1[n + 1], deltaT, 2);



//#pragma omp parallel for shared(K2,K3) schedule (static)

    for (int n = 1; n < count - 1; n++)
      K3[n] = Coefficient(temporary[n - 1], temporary[n], temporary[n + 1], deltaX, 1, K2[n - 1], K2[n], K2[n + 1], deltaT, 2);



//#pragma omp parallel for shared(K3,K4) schedule (static)

    for (int n = 1; n < count - 1; n++)
      K4[n] = Coefficient(temporary[n - 1], temporary[n], temporary[n + 1], deltaX, 1, K3[n - 1], K3[n], K3[n + 1], deltaT, 1);


//#pragma omp parallel for shared(K1,K2,K3,K4) schedule (static)

    for (int n = 1; n < count - 1; n++)
      current[n] = temporary[n] + (deltaT / 6) * (K1[n] + (2 * K2[n]) + (2 * K3[n]) + K4[n]);


    temp_array = temporary;
    temporary = current;
    current = temp_array;


   /* for (int s = 0; s < count; s++)
    {
      std::cout << current[s] << ' ';
    }
    std::cout << std::endl;*/


    i++;
  

  }

  /*t2 = omp_get_wtime();*/
  dt = t2 - t1;
  std::cout << "Time  " << dt << std::endl;
  Insert_In_File(current, count, write);
  

  std::cout << "end" << std::endl;
  write.close();

  delete[] temporary;
  //delete[] current;
  delete[] temp_array;


  delete[] K1;
  delete[] K2;
  delete[] K3;
  delete[] K4;

  return 0;

}