
#include "../include/program.h"

local bool read_mcp9(i2c_device * mcp9_i2c);

Sensor * init_mcp9(Sensor * mcp9) {

  mcp9 -> name = "MCP9808";
  mcp9 -> i2c = create_i2c_device(mcp9, read_mcp9, "An ambient thermometer");

  sensor_log_header(mcp9, RED);
  return mcp9;
}

bool read_mcp9(i2c_device * mcp9_i2c) {
  /*
   * Binary address 0011111; hex address 1F
   * Can read max every t_conv, or 250ms w/ 0.0625 *C accuracy
   * Can read max every 30ms w/ 0.5 *C accuracy
   *
   * A high-to-low transition of the SDA line (while SCL is high) is the Start
   * condition. All data transfers must be preceded by a Start condition from
   * the master. A low-to-high transition of the SDA line (while SCL is high)
   * signifies a Stop condition.
   */

  Sensor * mcp9 = mcp9_i2c -> sensor;

  // Initialize array of past values
  static double T_array[MCP9_SMOOTHING_SAMPLES] = {0};

  memmove(&T_array[1], &T_array[0], sizeof(double)*(MCP9_SMOOTHING_SAMPLES-1));

  uint8 read_raws[2];
  double temperature;

  if (!i2c_read_bytes(mcp9_i2c, 0x05, read_raws, 2)) return false;

  uint8 upper = read_raws[0];
  uint8 lower = read_raws[1];

  if (upper & 0x10) {
    upper = upper & 0x0F;  // Mask last 4 bits
    temperature = 256 - ((upper << 4) + (lower >> 4));  // Get ambient temp. (-)
    temperature -= (lower & 0x0F) * 0.0625f;            // Get decimal value
  } else {
    upper = upper & 0x0F;  // Mask last 4 bits
    temperature = (upper << 4) + (lower >> 4);          // Get ambient temp. (+)
    temperature += (lower & 0x0F) * 0.0625f;            // Get decimal value
  }

  double total = 0;
  double filled = 0;
  T_array[0] = temperature;

  for (int i = 0; i < MCP9_SMOOTHING_SAMPLES; i++) {
    if (T_array[i] > 0) {
      filled += (MCP9_SMOOTHING_SAMPLES - i);
      total += ((MCP9_SMOOTHING_SAMPLES - i) * T_array[i]);
    }
  }

  bind_stream(mcp9, temperature, MCP9_MEASURE_TEMPERATURE);
  bind_stream(mcp9, total/filled, MCP9_MEASURE_SMOOTH);

  sensor_log_outputs(mcp9, mcp9_i2c -> log, NULL);
  sensor_process_triggers(mcp9);
  return true;
}
