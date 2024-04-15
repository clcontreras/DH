#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

// This function returns the Debye-Hueckel g(r) between
// two ions with charges z1 and z2 and a separation distance r
// The function receives the valences, Bjerrum length, and Debye-Hueckel
// parameter as arguments
double dh(double z1, double z2, double r, double lB, double kappa) {
  return 1 - z1 * z2 * lB * exp(-kappa * r) / r;
}

// This function saves the cation-cation, anion-anion, cation-anion, and
// anion-cation g(r) functions to files
void save_gr(double z1, double z2, double lB_red, double kappa_red) {
  // We open the file g_11.dat to save the cation-cation g(r)
  std::ofstream g11("gr_11.dat");
  // We add a header to the file g_11.dat
  g11 << "# r g_11" << std::endl;
  g11 << "#" << std::endl;
  // This loop calculates the cation-cation g(r) and saves it to the file
  // g_11.dat with 12 decimal places of precision for g(r)
  for (double r = 0.0008; r < 5; r += 0.0005) {
    g11 << r << " " << std::setprecision(12)
        << dh(z1, z1, r, lB_red, kappa_red) << std::endl;
  }
  // We close the file g_11.dat
  g11.close();
  // We open the file g_22.dat to save the anion-anion g(r)
  std::ofstream g22("gr_22.dat");
  // We add a header to the file g_22.dat
  g22 << "# r g_22" << std::endl;
  g22 << "#" << std::endl;
  // This loop calculates the anion-anion g(r) and saves it to the file
  // g_22.dat

  for (double r = 0.0008; r < 5; r += 0.0005) {
    g22 << r << " " << std::setprecision(12)
        << dh(z2, z2, r, lB_red, kappa_red) << std::endl;
  }
  // We close the file g_22.dat
  g22.close();
  // We open the file g_12.dat to save the cation-anion g(r)
  std::ofstream g12("gr_12.dat");
  // We add a header to the file g_12.dat
  g12 << "# r g_12" << std::endl;
  g12 << "#" << std::endl;
  // This loop calculates the cation-anion g(r) and saves it to the file
  // g_12.dat
  for (double r = 0.0008; r < 5; r += 0.0005) {
    g12 << r << " " << std::setprecision(12)
        << dh(z1, z2, r, lB_red, kappa_red) << std::endl;
  }
  // We close the file g_12.dat
  g12.close();
  // We open the file g_21.dat to save the anion-cation g(r)
  std::ofstream g21("gr_21.dat");
  // We add a header to the file g_21.dat
  g21 << "# r g_21" << std::endl;
  g21 << "#" << std::endl;
  // This loop calculates the anion-cation g(r) and saves it to the file
  // g_21.dat
  for (double r = 0.0008; r < 5; r += 0.0005) {
    g21 << r << " " << std::setprecision(12)
        << dh(z2, z1, r, lB_red, kappa_red) << std::endl;
  }
  // We close the file g_21.dat
  g21.close();
}

// This function calculates the cation-cation, anion-anion, cation-anion,
// and anion-cation S(q) functions, and saves them to files.
void save_Sq(double kappa_red) {
  // We open the file S_11.dat to save the cation-cation S(q)
  std::ofstream S11("Sq_11.dat");
  // We add a header to the file S_11.dat
  S11 << "# q S_11" << std::endl;
  S11 << "#" << std::endl;
  // This loop calculates the cation-cation S(q) and saves it to the file
  // S_11.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    S11 << q << " " << 1 - (kappa_red * kappa_red / 2) /
                           (q * q + kappa_red * kappa_red) << std::endl;
  }
  // We close the file S_11.dat
  S11.close();
  // We open the file S_22.dat to save the anion-anion S(q)
  std::ofstream S22("Sq_22.dat");
  // We add a header to the file S_22.dat
  S22 << "# q S_22" << std::endl;
  S22 << "#" << std::endl;
  // This loop calculates the anion-anion S(q) and saves it to the file
  // S_22.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    S22 << q << " " << 1 - (kappa_red * kappa_red / 2) /
                           (q * q + kappa_red * kappa_red) << std::endl;
  }
  // We close the file S_22.dat
  S22.close();
  // We open the file S_12.dat to save the cation-anion S(q)
  std::ofstream S12("Sq_12.dat");
  // We add a header to the file S_12.dat
  S12 << "# q S_12" << std::endl;
  S12 << "#" << std::endl;
  // This loop calculates the cation-anion S(q) and saves it to the file
  // S_12.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    S12 << q << " "
        << (kappa_red * kappa_red / 2) / (q * q + kappa_red * kappa_red)
        << std::endl;
  }
  // We close the file S_12.dat
  S12.close();
  // We open the file S_21.dat to save the anion-cation S(q)
  std::ofstream S21("Sq_21.dat");
  // We add a header to the file S_21.dat
  S21 << "# q S_21" << std::endl;
  S21 << "#" << std::endl;
  // This loop calculates the anion-cation S(q) and saves it to the file
  // S_21.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    S21 << q << " "
        << (kappa_red * kappa_red / 2) / (q * q + kappa_red * kappa_red)
        << std::endl;
  }
}

//This function calculates the cation-cation, anion-anion, cation-anion,
//and anion-cation C(q) functions, and saves them to files.
void save_Cq(double kappa_red) {
  // We open the file C_11.dat to save the cation-cation C(q)
  std::ofstream C11("Cq_11.dat");
  // We add a header to the file C_11.dat
  C11 << "# q C_11" << std::endl;
  C11 << "#" << std::endl;
  // This loop calculates the cation-cation C(q) and saves it to the file
  // C_11.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    C11 << q << " " << -(kappa_red * kappa_red / 2) / (q * q)
        << std::endl;
  }
  // We close the file C_11.dat
  C11.close();
  // We open the file C_22.dat to save the anion-anion C(q)
  std::ofstream C22("Cq_22.dat");
  // We add a header to the file C_22.dat
  C22 << "# q C_22" << std::endl;
  C22 << "#" << std::endl;
  // This loop calculates the anion-anion C(q) and saves it to the file
  // C_22.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    C22 << q << " " << -(kappa_red * kappa_red / 2) / (q * q)
        << std::endl;
  }
  // We close the file C_22.dat
  C22.close();
  // We open the file C_12.dat to save the cation-anion C(q)
  std::ofstream C12("Cq_12.dat");
  // We add a header to the file C_12.dat
  C12 << "# q C_12" << std::endl;
  C12 << "#" << std::endl;
  // This loop calculates the cation-anion C(q) and saves it to the file
  // C_12.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    C12 << q << " " << (kappa_red * kappa_red / 2) / (q * q)
        << std::endl;
  }
  // We close the file C_12.dat
  C12.close();
  // We open the file C_21.dat to save the anion-cation C(q)
  std::ofstream C21("Cq_21.dat");
  // We add a header to the file C_21.dat
  C21 << "# q C_21" << std::endl;
  C21 << "#" << std::endl;
  // This loop calculates the anion-cation C(q) and saves it to the file
  // C_21.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    C21 << q << " " << (kappa_red * kappa_red / 2) / (q * q)
        << std::endl;
  }
  // We close the file C_21.dat
  C21.close();
}

// This is the Oseen function for hydrodynamic interactions
double FO(double x) {
  return 3 * (x + (x * x - 1) * atan(x)) / (4 * (x * x * x));
}

// This function calculates the cation-cation, anion-anion, cation-anion,
// and anion-cation H(q) functions, and saves them to files.
void save_Hq(double kappa_red, double Dkappa, double D0) {
  // We open the file H_11.dat to save the cation-cation H(q) and a file
  // Hred_11.dat to save the reduced cation-cation H(q)
  std::ofstream H11("Hq_11.dat");
  std::ofstream Hred11("Hqred_11.dat");
  // We add a header to the file H_11.dat
  H11 << "# q H_11" << std::endl;
  H11 << "#" << std::endl;
  // We add a header to the file Hred_11.dat
  Hred11 << "# q Hred_11" << std::endl;
  Hred11 << "#" << std::endl;
  // This loop calculates the cation-cation H(q) and saves it to the file
  // H_11.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    H11 << q << " " << D0 - (Dkappa / 2) * FO(q / kappa_red)
        << std::endl;
    Hred11 << q << " " << 1 - (Dkappa / (2 * D0)) * FO(q / kappa_red)
           << std::endl;
  }
  // We close the file H_11.dat and Hred_11.dat
  H11.close();
  Hred11.close();
  // We open the file H_22.dat to save the anion-anion H(q) and a file
  // Hred_22.dat to save the reduced anion-anion H(q)
  std::ofstream H22("Hq_22.dat");
  std::ofstream Hred22("Hqred_22.dat");
  // We add a header to the file H_22.dat
  H22 << "# q H_22" << std::endl;
  H22 << "#" << std::endl;
  // We add a header to the file Hred_22.dat
  Hred22 << "# q Hred_22" << std::endl;
  Hred22 << "#" << std::endl;
  // This loop calculates the anion-anion H(q) and saves it to the file
  // H_22.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    H22 << q << " " << D0 - (Dkappa / 2) * FO(q / kappa_red)
        << std::endl;
    Hred22 << q << " " << 1 - (Dkappa / (2 * D0)) * FO(q / kappa_red)
           << std::endl;
  }
  // We close the file H_22.dat and Hred_22.dat
  H22.close();
  Hred22.close();
  // We open the file H_12.dat to save the cation-anion H(q) and a file
  // Hred_12.dat to save the reduced cation-anion H(q)
  std::ofstream H12("Hq_12.dat");
  std::ofstream Hred12("Hqred_12.dat");
  // We add a header to the file H_12.dat
  H12 << "# q H_12" << std::endl;
  H12 << "#" << std::endl;
  // We add a header to the file Hred_12.dat
  Hred12 << "# q Hred_12" << std::endl;
  Hred12 << "#" << std::endl;
  // This loop calculates the cation-anion H(q) and saves it to the file
  // H_12.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    H12 << q << " " << (Dkappa / 2) * FO(q / kappa_red) << std::endl;
    Hred12 << q << " " << (Dkappa / (2 * D0)) * FO(q / kappa_red)
           << std::endl;
  }
  // We close the file H_12.dat and Hred_12.dat
  H12.close();
  Hred12.close();
  // We open the file H_21.dat to save the anion-cation H(q) and a file
  // Hred_21.dat to save the reduced anion-cation H(q)
  std::ofstream H21("Hq_21.dat");
  std::ofstream Hred21("Hqred_21.dat");
  // We add a header to the file H_21.dat
  H21 << "# q H_21" << std::endl;
  H21 << "#" << std::endl;
  // We add a header to the file Hred_21.dat
  Hred21 << "# q Hred_21" << std::endl;
  Hred21 << "#" << std::endl;
  // This loop calculates the anion-cation H(q) and saves it to the file
  // H_21.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    H21 << q << " " << (Dkappa / 2) * FO(q / kappa_red) << std::endl;
    Hred21 << q << " " << (Dkappa / (2 * D0)) * FO(q / kappa_red)
           << std::endl;
  }
  // We close the file H_21.dat and Hred_21.dat
  H21.close();
  Hred21.close();
}

double lambdaEL(double y, double b, double Dkappa, double D0) {
  return (y * y + b * b) * (1 - (Dkappa / D0) * FO(y / b));
}

double expD0(double y, double t) {
  return exp(-y * y * t);
}

double expEL(double y, double b, double Dkappa, double D0, double t) {
  return (y * y * exp(-lambdaEL(y, b, Dkappa, D0) * t) -
          (y * y + b * b) * expD0(y, t)) / (2 * (y * y + b * b));
}


// This function calculates the cation-cation, anion-anion, cation-anion,
// and anion-cation FS(q,t) functions, and saves them to files.

void save_FSqt(double kappa_red, double Dkappa, double D0) {
  // We open the file FS_11.dat to save the cation-cation FS(q,t)
  std::ofstream FS11("FSqt_11.dat");
  // We add a header to the file FS_11.dat
  FS11 << "# q FS_11" << std::endl;
  FS11 << "#" << std::endl;
  // This loop calculates the cation-cation FS(q,t) and saves it to the file
  // FS_11.dat
  double t = 1.0;
  double fsaux;
  for (double q = 0.001; q < 24; q += 0.02) {
    fsaux = expD0(q, t) + expEL(q, kappa_red, Dkappa, D0, t);
    FS11 << q << " " << fsaux << std::endl;
  }
  // We close the file FS_11.dat
  FS11.close();

  // We open the file FS_12.dat to save the cation-anion FS(q,t)
  std::ofstream FS12("FSqt_12.dat");
  // We add a header to the file FS_22.dat
  FS12 << "# q FS_12" << std::endl;
  FS12 << "#" << std::endl;
  // This loop calculates the cation-anion FS(q,t) and saves it to the file
  // FS_12.dat
  for (double q = 0.001; q < 24; q += 0.02) {
    fsaux = -expEL(q, kappa_red, Dkappa, D0, t);
    FS11 << q << " " << fsaux << std::endl;
  }
  // We close the file FS_12.dat
  FS12.close();
}


int main() {
  // Define the valences of the ions
  double z1 = 1;
  double z2 = -1;
  // Define the Bjerrum length
  double lB = 7.13997e-10;
  // Define the molar concentration of each ion in mol/l
  double n1M = 0.005;
  double n2M = 0.005;
  // Convert the molar concentration to number density
  double n1 = n1M * 6.022e23 * 1e3;
  double n2 = n2M * 6.022e23 * 1e3;
  // Define the Debye-Hueckel parameter
  double kappa = sqrt(4 * M_PI * lB * (n1 * z1 * z1 + n2 * z2 * z2));
  // Average ion separation
  double d = 1 / pow(n1 + n2, 1. / 3.);
  // Reduce the Debye-Hueckel parameter to the average ion separation
  double kappa_red = kappa * d;
  // Reduce the Bjerrum length to the average ion separation
  double lB_red = lB / d;
  // Save the g(r) functions to files
  save_gr(z1, z2, lB_red, kappa_red);
  // Save the S(q) functions to files
  save_Sq(kappa_red);
  // Save the C(q) functions to files
  save_Cq(kappa_red);
  double K_B = 1.381e-23;
  double T = 298.15;
  double eta = 0.00085;
  double Dkappa = K_B * T * kappa / (6 * M_PI * eta);
  double D0 = 2.032e-9;
  // Save the H(q) functions to files
  save_Hq(kappa_red, Dkappa, D0);
  // Save the FS(q,t) functions to files
  save_FSqt(kappa_red, Dkappa, D0);

  std::cout << "********************************************************"
            << std::endl;
  std::cout << "Reduced short-time cation-cation partial mobility is "
            << 1 - (Dkappa / (2 * D0)) << std::endl;
  std::cout << "Reduced short-time cation-anion partial mobility is "
            << (Dkappa / (2 * D0)) << std::endl;

  std::cout << "********************************************************"
            << std::endl;

  std::cout << "Reduced electrophoretic mobility of cations is "
            << 1 - (Dkappa / D0) << std::endl;
  
  std::cout
      << "Mcirr(q=0,t) reducida integrada en el tiempo es para cation-cation "
      << ((2 - sqrt(2)) / 12) * kappa * lB << std::endl;

  return 0;
}
