

#include <stdlib.h>
#include <stdio.h>

#include "ad15.h"
#include "ds32.h"

#include "../system/i2c.h"
#include "../system/color.h"


// information regarding this can be found pages 18-19 of the datasheet
// more options exist in the datasheet, but this should cover most needs

#define AD15_POINTER_REG            0x00     // pointer register
#define AD15_CONFIG_REG             0x01     // configuration register

#define AD15_OS_NO_EFFECT           0b0      // do nothing
#define AD15_OS_BEGIN_CONVERSION    0b1      // begin single conversion

#define AD15_MUX_DIFFERENTIAL_AIN01 0b000    // differential between A0 and A1
#define AD15_MUX_DIFFERENTIAL_AIN23 0b011    // differential between A2 and A3
#define AD15_MUX_SINGLE_AIN0        0b100    // single channel muxed on A0
#define AD15_MUX_SINGLE_AIN1        0b101    // single channel muxed on A1
#define AD15_MUX_SINGLE_AIN2        0b110    // single channel muxed on A2
#define AD15_MUX_SINGLE_AIN3        0b111    // single channel muxed on A3

#define AD15_PGA_6144V              0b000    // analog range +/- 6.144V
#define AD15_PGA_4096V              0b001    // analog range +/- 4.096V
#define AD15_PGA_2048V              0b010    // analog range +/- 2.048V
#define AD15_PGA_1024V              0b011    // analog range +/- 1.024V
#define AD15_PGA_0512V              0b100    // analog range +/- 0.512V
#define AD15_PGA_0256V              0b101    // analog range +/- 0.256V

#define AD15_MODE_CONTINUOUS        0b0      // continually sample
#define AD15_MODE_LOW_POWER         0b1      // sample only when asked

#define AD15_RATE_8HZ               0b000    // sample at   8Hz
#define AD15_RATE_16HZ              0b001    // sample at  16Hz
#define AD15_RATE_32HZ              0b010    // sample at  32Hz
#define AD15_RATE_64HZ              0b011    // sample at  64Hz
#define AD15_RATE_128HZ             0b100    // sample at 128Hz
#define AD15_RATE_250HZ             0b101    // sample at 250Hz
#define AD15_RATE_475HZ             0b110    // sample at 475Hz
#define AD15_RATE_860HZ             0b111    // sample at 860Hz

#define AD15_COMP_TRADITIONAL       0b0      // use traditional comparator with hysteresis
#define AD15_COMP_WINDOW            0b1      // use window comparator

#define AD15_POL_ACTIVE_LOW         0b0      // set comparator to active low
#define AD15_POL_ACTIVE_HIGH        0b1      // set comparator to active high

#define AD15_LAT_DONT_LATCH         0b0      // use a non-latching comparator
#define AD15_LAT_DO_LATCH           0b1      // use a latching comparator

#define AD15_QUE_ASSERT_ONE         0b00     // assert comparator after a single conversion
#define AD15_QUE_ASSERT_TWO         0b01     // assert comparator after two conversions
#define AD15_QUE_ASSERT_FOUR        0b10     // assert comparator after four conversions
#define AD15_QUE_DISABLE            0b11     // disable comparator

void free_ad15(Sensor * ad15);
bool read_ad15(i2c_device * ad15_i2c);
void configure_ad15(Sensor * ad15);

Sensor * init_ad15(uint8 address, char * title, List * modes, List * names) {
  
  Sensor * ad15 = malloc(sizeof(Sensor));
  
  ad15 -> name = "ADS1115";
  ad15 -> free = free_ad15;
  
  ad15 -> i2c = create_i2c_device(ad15, address, read_ad15, 10);    // 10ms between reads
  
  char file_name[32];
  
  sprintf(file_name, "logs/ad15-%x.log", address);
  
  FILE * file = fopen(file_name, "a");
  
  ad15 -> i2c -> file = file;

  fprintf(file, RED "\n\n");
  fprintf(file, "ADS115 - %s\n", title);
  fprintf(file, "Start time %s\n", formatted_time);
  
  for (iterate(names, char *, column_name))
    fprintf(file, "%s\t", column_name);
  
  fprintf(file, "\n" RESET);
  
  // set up the configuration  (page 19)

  AD15_Config * sensor_config = malloc(sizeof(AD15_Config));
  
  sensor_config -> COMP_QUE  = AD15_QUE_DISABLE;
  sensor_config -> COMP_LAT  = AD15_LAT_DONT_LATCH;
  sensor_config -> COMP_POL  = AD15_POL_ACTIVE_LOW;
  sensor_config -> COMP_MODE = AD15_COMP_TRADITIONAL;
  sensor_config -> DATA_RATE = AD15_RATE_128HZ;
  
  sensor_config -> MODE      = AD15_MODE_CONTINUOUS;
  sensor_config -> PGA       = AD15_PGA_6144V;
  sensor_config -> MUX       = AD15_MUX_SINGLE_AIN0;
  sensor_config -> OS        = AD15_OS_BEGIN_CONVERSION;
  
  sensor_config -> modes     = modes;
  
  ad15 -> data = sensor_config;
  
  configure_ad15(ad15);
  
  return ad15;
}

void configure_ad15(Sensor * ad15) {
  
  AD15_Config * sensor_config = ad15 -> data;
  
  uint8 request[3] = {
    AD15_CONFIG_REG,               // point to config register
    sensor_config -> high_byte,    // fill in the request according to the config
    sensor_config -> low_byte,     // -------------------------------------------
  };
  
  i2c_raw_write(ad15 -> i2c, request, 3);    
}

bool read_ad15(i2c_device * ad15_i2c) {
  
  Sensor * ad15 = ad15_i2c -> sensor;
  
  AD15_Config * config = ad15 -> data;
  
  for (iterate(config -> modes, uint8, mode)) {
    
    uint8 ad15_raws[2];
    
    config -> MUX = mode << 4;    // gotta shift to match the MUX field
    
    configure_ad15(ad15);
    
    i2c_read_bytes(ad15_i2c, 0x00, ad15_raws, 2);
    
    uint16 counts = (ad15_raws[0] << 8) | ad15_raws[1];
    
    fprintf(ad15_i2c -> file, "%d\t", counts);
  }
  
  fprintf(ad15_i2c -> file, "\n");
  
  return true;
}

void free_ad15(Sensor * ad15) {
  
}
