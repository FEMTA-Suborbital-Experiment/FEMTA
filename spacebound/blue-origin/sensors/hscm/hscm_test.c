#define HSCM_ADDRESS 0x28

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pigpio.h>


int read_hscm(int handle) {
  /* 
   * Binary address 00101000; hex address 28
   *
   * To read, send address followed by read bit. Then the pressure sensor will
   * return 4 bytes -- the first 2 are MSB and LSB of pressure, respectively.
   * The 3rd is the MSB of temperature, and the last one has 3 LSB.
   */
  
  char read_raws[4];
  int p_counts;
  int t_counts;
  //double pressure;
  double temperature;
  
  if (!i2cReadDevice(handle, read_raws, 4)) return 0;
  
  uint8_t p_upper = read_raws[0];
  uint8_t p_lower = read_raws[1];
  uint8_t t_upper = read_raws[2];
  uint8_t t_lower = read_raws[3] >> 5;

  p_counts = (p_upper << 8) + p_lower;
  t_counts = (t_upper << 3) + t_lower;
  temperature = (((double) t_counts) / 2047 * 200) - 50;  // conversion taken from datasheet
  
  printf("pressure: %d, temperature %f\n", p_counts, temperature);
  
  return 0;
}

int main() {

  if (gpioInitialise() < 0) return 1;

  int handle = i2cOpen(1, HSCM_ADDRESS, 0);
  

  for (int i = 0; i < 10; i++) {
    if (read_hscm(handle)) break;
    sleep(1);
  }

  return 0;

  i2cClose(handle);
  gpioTerminate();
}
