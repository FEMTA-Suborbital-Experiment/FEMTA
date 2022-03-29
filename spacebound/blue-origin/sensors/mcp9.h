#pragma once

/**
 * Sensor:
 *   0.5C Maximum Accuracy Digital Temperature Sensor
 *   MCP9808
 *   Connected via I2C
 *
 * Datasheet:
 *   http://ww1.microchip.com/downloads/en/DeviceDoc/25095A.pdf
 *
 * Author:
 *   Mark Hartigan
 */

#define MCP9_ADDRESS 0x18

#define MCP9_MEASURE_TEMPERATURE 0
#define MCP9_SMOOTH_TEMPERATURE 0
