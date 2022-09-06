

#include <stdlib.h>
#include <stdio.h>

#include "adxl.h"
#include "ds32.h"
#include "fram.h"
#include "list.h"
#include "sensor.h"


void init_sensors() {
  
  sensors = create_list(SLL, NULL);
  
  list_insert(sensors, init_ds32());    // first so we can get the time
  list_insert(sensors, init_adxl());
  list_insert(sensors, init_fram());
  
  list_insert(sensors, init_ads(ADS_DEFAULT_GND,
                                init_noise(),
                                init_noise(),
                                init_noise(),
                                init_noise()
                                )
              );

  list_insert(sensors, init_ads(ADS_DEFAULT_VDD,
                                init_noise(),
                                init_noise(),
                                init_noise(),
                                init_noise()
                                )
              );

  list_insert(sensors, init_ads(ADS_DEFAUT_SDA,
                                init_noise(),
                                init_noise(),
                                init_noise(),
                                init_noise()
                                )
              );

  list_insert(sensors, init_ads(ADS_DEFAUT_SDL,
                                init_noise(),
                                init_noise(),
                                init_noise(),
                                init_noise()
                                )
              );
}

Sensor * create_sensor(char * name, sensor_free freer) {

  Sensor * sensor = malloc(sizeof(Sensor));
  
  sensor -> name = name;
  sensor -> free = freer;

  return sensor;
}


void destroy_sensors() {
  
  
}

