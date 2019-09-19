#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mcp9.h"

#include "../system/color.h"
#include "../system/gpio.h"
#include "../system/i2c.h"

void free_mcp9(Sensor * sensor);
bool read_mcp9(i2c_device * mcp9_i2c);

Sensor * init_mcp9(ProtoSensor * proto) {

  Sensor * mcp9 = sensor_from_proto(proto);

  mcp9 -> name = "MCP9808";
  mcp9 -> free = free_mcp9;

  mcp9 -> i2c = create_i2c_device(mcp9, proto -> address, read_mcp9, proto -> hertz);
  mcp9 -> i2c -> log = fopen("logs/mcp9.log", "a");
  
  mcp9 -> print = proto -> print;

  setlinebuf(mcp9 -> i2c -> log);  // write out every read

  fprintf(mcp9 -> i2c -> log, GREEN "\nMCP9808\n" RESET);
  printf("Started " GREEN "%s " RESET "at " YELLOW "%dHz " RESET "on " BLUE "0x%x " RESET,
	 mcp9 -> name, proto -> hertz, proto -> address);

  printf("logged in logs/mcp9.log\n");
  
  return mcp9;
}

bool read_mcp9(i2c_device * mcp9_i2c) {
  /*
   * Binary address 0011000; hex address 18
   * 0th bit: 1 is read, 0 is write
   * Can read max every t_conv, or 250ms w/ 0.0625 *C accuracy
   * Can read max every 30ms w/ 0.5 *C accuracy
   */
  
  Sensor * mcp9 = mcp9_i2c -> sensor;

  uint8 read_raws[2];
  int sign = 1;
  double temp;
  
  if (!i2c_read_bytes(mcp_i2c, 0x05, read_raws, 2)) return false;

  upper = read_raws[0];
  lower = read_raws[1];

  if (upper & 0x10) {
    upper = upper & 0x0F;  // Mask last 4 bits
    temp = 256 - (upper << 4 + lower >> 4);  // Get ambient temp. (-)
  }
  else {
    upper = upper & 0x0F;  // Mask last 4 bits
    temp = upper << 4 + lower >> 4;  // Get ambient temp. (+)
  }

  fprintf(mcp9 -> log, "%d\n", temp);

  return true;
}

void free_mcp9(Sensor * ds32) {
  // TODO: Does the MCP9808 have special release conditions?
}
