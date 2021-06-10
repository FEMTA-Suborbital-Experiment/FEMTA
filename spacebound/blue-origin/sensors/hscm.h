#pragma once

/**
 * Sensor:
 *   SMT Board Mount Pressure Sensor 0-15 psi
 *   HSCMNNN015PA2A3
 *   Connected via I2C
 *
 * Datasheet:
 *   https://www.mouser.com/datasheet/2/187/honeywell-sensing-trustability-hsc-series-high-acc-1228679.pdf
 *   https://prod-edam.honeywell.com/content/dam/honeywell-edam/sps/siot/en-us/products/sensors/pressure-sensors/board-mount-pressure-sensors/common/documents/sps-siot-i2c-comms-digital-output-pressure-sensors-tn-008201-3-en-ciid-45841.pdf
 *
 * Author:
 *   Mark Hartigan
 */

#define HSCM_ADDRESS 0x28

#define HSCM_MEASURE_PRESSURE    0
#define HSCM_MEASURE_TEMPERATURE 1
