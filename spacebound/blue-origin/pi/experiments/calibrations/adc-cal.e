define enter start;

Sensor ad15_sda 4Hz {
  
  where A0 is numa
  where A1 is numb;
  where A2 is numc;
  where A3 is numd; 

  [print | green, numa  | 4];
  [print | green, A0    | 4];
  [print | red, numb    | 4];
  [print | red, A1      | 4];
  [print | blue, numc   | 4];
  [print | blue, A2     | 4];
  [print | orange, numd | 4];
  [print | orange, A3   | 4];
}
