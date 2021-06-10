
#include "../include/program.h"

local bool read_hscm(i2c_device * hscm_i2c);

Sensor * init_hscm(Sensor * hscm) {
  
  hscm -> name = "HSCMNNN015PA2A3";
  hscm -> i2c = create_i2c_device(hscm, read_hscm, "An ambient pressure sensor");
  
  sensor_log_header(hscm, RED);
  return hscm;
}

bool read_hscm(i2c_device * hscm_i2c) {
  /* 
   * Binary address 00101000; hex address 28
   *
   * To read, send address followed by read bit. Then the pressure sensor will
   * return 4 bytes -- the first 2 are MSB and LSB of pressure, respectively.
   * The 3rd is the MSB of temperature, and the last one has 3 LSB.
   */
  
  Sensor * hscm = hscm_i2c -> sensor;
  
  uint8 read_raws[4];
  int p_counts;
  int t_counts;
  //double pressure;
  double temperature;
  
  if (!i2c_raw_read(hscm_i2c, read_raws, 4)) return false;
  
  uint8 p_upper = read_raws[0];
  uint8 p_lower = read_raws[1];
  uint8 t_upper = read_raws[2];
  uint8 t_lower = read_raws[3] >> 5;

  p_counts = (p_upper << 8) + p_lower;
  t_counts = (t_upper << 3) + t_lower;

  temperature = (((double) t_counts) / 2047 * 200) - 50;  // conversion taken from datasheet
  
  bind_stream(hscm, p_counts, HSCM_MEASURE_PRESSURE);
  bind_stream(hscm, temperature, HSCM_MEASURE_TEMPERATURE);
  
  sensor_log_outputs(hscm, hscm_i2c -> log, NULL);
  sensor_process_triggers(hscm);
  return true;
}
