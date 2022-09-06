
ds18 1Hz {
  
  //  Temperature less_than   73F trigger {};             // should result in error  
  //  Temperature less_than 32.5C trigger {} singular;    // ----------------------
  
  Temperature less_than 36C trigger {
      +12,
      -13 [pulse :: 300ms],
      +14 [pulse :: 250ms]
  } singular;
	
  Temperature more_than 216K trigger {
    -8, +32
  } forever, reverses;
  
  [calibrate : Temperature, hart : 48.0, 24.5, 0.1];

  [print];
  
[smooth :: 0.2];


  [print  :: 5Hz];
};
