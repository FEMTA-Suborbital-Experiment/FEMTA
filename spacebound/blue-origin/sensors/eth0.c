
#include "../include/program.h"

local bool read_eth0(spi_device * eth0_spi);

Sensor * init_eth0(Sensor * eth0) {

    eth0 -> name = "ENC28J60";
    eth0 -> spi = create_spi_device(eth0, read_eth0, "The Ethernet interface to New Shepard");

    sensor_log_header(eth0, BLUE);

    micro_sleep(300); // Allow Ethernet module oscillator to stabilize (datasheet sec 2.2)
    /* Send any initialization commands to the module here */

    return eth0;
}

bool read_eth0(spi_device * eth0_spi) {

    Sensor * eth0 = eth0_spi -> sensor;

    int flight_state = 0; // Figure out how to get this from the sensor, obviously

    bind_stream(eth0, flight_state, ETH0_MEASURE_STATE);

    sensor_log_outputs(eth0, eth0_spi -> log, NULL);
    sensor_process_triggers(eth0);
    return true;
}
