#pragma once

#include "../include/headers.h"

typedef struct spi_device {    // mechanism for spi communication
  
  Sensor * sensor;             // the sensor bound to this device

  FILE * log;                  // the log file

  uint8 channel;               // the SPI channel (0 or 1 unless we use aux peripheral)

  int interval;                // time span between reads in ms
  int hertz_numerator;         // schedule frequency in hertz
  int hertz_denominator;       // engenders fractional frequency through deferrals
  int count;                   // counts since last read

  int handle;                  // the device's pigpio handle

  i2c_reader read;

} spi_device;
