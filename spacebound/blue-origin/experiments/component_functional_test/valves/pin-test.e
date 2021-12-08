/*
**  
*/

define enter test;
define leave finished;

Sensor ds32 1Hz {
    [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
    [conversions | Time, raw, s       |                  ];
	
	if (State test | Time > 0s) {
		// Start with all pins off
		set pin 12 neg; //large 1 open
		set pin 16 neg; //large 1 close
		set pin 20 neg; //large 2 open
		set pin 21 neg; //large 2 close
		set pin 14 neg; //small 1 forward
		set pin 15 neg; //small 1 reverse
		set pin 18 neg; //small 2 forward
		set pin 23 neg; //small 2 reverse
		set pin 24 neg; //small 3 forward
		set pin 25 neg; //small 3 reverse
		
		// Then go through each and pulse it on for one second
		set pin 12 pos after 1s;
		set pin 12 neg after 2s;
		
		set pin 16 pos after 3s;
		set pin 16 neg after 4s;
		
		set pin 20 pos after 5s;
		set pin 20 neg after 6s;
		
		set pin 21 pos after 7s;
		set pin 21 neg after 8s;
		
		set pin 14 pos after 9s;
		set pin 14 neg after 10s;
		
		set pin 15 pos after 11s;
		set pin 15 neg after 12s;
		
		set pin 18 pos after 13s;
		set pin 18 neg after 14s;
		
		set pin 23 pos after 15s;
		set pin 23 neg after 16s;
		
		set pin 24 pos after 17s;
		set pin 24 neg after 18s;
		
		set pin 25 pos after 19s;
		set pin 25 neg after 20s;
		
		leave test after 250ms;
		enter finished after 500ms;
	}
}
