#pragma once

/**
 * Sensor:
 *   Ethernet Module
 *   Microchip ENC28J60
 *   Mounted on controller module
 *   Connected via SPI
 * 
 * Datasheet:
 *   https://ww1.microchip.com/downloads/en/DeviceDoc/39662e.pdf
 *   https://www.openimpulse.com/blog/products-page/product-category/enc28j60-ethernet-module/
 * 
 * Author:
 *   Evan Rittner
 */

#include "../include/headers.h"

#define ETH0_MEASURE_STATE 0
#define ETH0_CHANNEL       0
