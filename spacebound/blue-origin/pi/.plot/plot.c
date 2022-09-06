
#include <stdlib.h>
#include <stdio.h>

#include <curses.h>
#include <string.h>

#include "../structures/list.h"


char * plot_title;

List * data[3];
char * axis_names[3];
char axes_enabled;

void update_plot() {

  mvprintw(0, COLS / 2 - 16, "%s         ", plot_title);
  
  double min_value = *(double *) data[0] -> head -> value;
  double max_value = *(double *) data[0] -> head -> value;
  
  for (uchar l = 0; l < axes_enabled; l++) {

    bool first = false;
    
    for (Node * node = data[l] -> head; node && (!first || node != data[l] -> head); node = node -> next) {

      first = true;
      
      double value = *(double *) node -> value;
      
      if (value < min_value) min_value = value;
      if (value > max_value) max_value = value;
    }
  }
  
  if (min_value == max_value) {
    // Uniform data, bail out
    
    refresh();
    return;
  }
  
  // Draw y labels
  for (uchar l = 2; l < LINES - 5; l++) {
    
    double left_label = (double) (l - 2) * (max_value - min_value) / (LINES - 5 - 2) + min_value;
    
    mvprintw(LINES - 4 - l, 2 + (left_label >= 0), "%.2lf", left_label);
  }
  
  for (uchar y = 2; y < LINES - 5; y++) {
    for (uchar x = 10; x < COLS - 3; x++) {
      mvaddch(y, x, ' ');
    }
  }
  
  // Plot points    
  Node * nodes[axes_enabled];
  for (uchar n = 0; n < axes_enabled; n++) {
    nodes[n] = data[n] -> head;
  }
  
  short number_of_data_points_graphable = COLS - 14;
  
  for (uchar c = 0; c < number_of_data_points_graphable; c++) {
    
    if (c > data[0] -> elements - 1) break;    // not enough data to continue plotting
    
    uchar heights[axes_enabled];
    for (uchar h = 0; h < 3; h++) {
      double relative = *(double *) nodes[h] -> value;
      
      relative = (relative - min_value) / (max_value - min_value);
      
      if (relative < 0) relative = 0;
      if (relative > 1) relative = 1;
      
      relative *= LINES - 8;
      
      heights[h] = (uchar) relative;
    }

    char glyph = '*';
    
    if (axes_enabled == 1) {
      attron(COLOR_PAIR(4));
      mvaddch(LINES - 6 - heights[0], 9 + number_of_data_points_graphable - c, glyph);
      attroff(COLOR_PAIR(4));
    }
    if (axes_enabled == 3) {
      
      attron(COLOR_PAIR(4));
      mvaddch(LINES - 6 - heights[0], 9 + number_of_data_points_graphable - c, glyph);
      attroff(COLOR_PAIR(4));
      
      attron(COLOR_PAIR(2));
      mvaddch(LINES - 6 - heights[1], 9 + number_of_data_points_graphable - c, glyph);
      attroff(COLOR_PAIR(2));
      
      attron(COLOR_PAIR(6));
      mvaddch(LINES - 6 - heights[2], 9 + number_of_data_points_graphable - c, glyph);
      attroff(COLOR_PAIR(6));
      
      if (heights[0] == heights[1]) {
	attron(COLOR_PAIR(5));
	mvaddch(LINES - 6 - heights[0], 9 + number_of_data_points_graphable - c, glyph);
	attroff(COLOR_PAIR(5));
      }
      
      if (heights[1] == heights[2]) {
	attron(COLOR_PAIR(7));
	mvaddch(LINES - 6 - heights[1], 9 + number_of_data_points_graphable - c, glyph);
	attroff(COLOR_PAIR(7));
      }
      
      if (heights[0] == heights[2]) {
	attron(COLOR_PAIR(3));
	mvaddch(LINES - 6 - heights[0], 9 + number_of_data_points_graphable - c, glyph);
	attroff(COLOR_PAIR(3));
      }
      
      if (heights[0] == heights[1] && heights[0] == heights[2]) {
	attron(COLOR_PAIR(1));
	mvaddch(LINES - 6 - heights[0], 9 + number_of_data_points_graphable - c, glyph);
	attroff(COLOR_PAIR(1));
      }
    }
    
    /*for (uchar h = 0; h < axes_enabled; h++) {
      mvprintw(LINES - 7 - heights[h], number_of_data_points_graphable - c, "*");
      }*/
    
    for (uchar n = 0; n < axes_enabled; n++) {
      nodes[n] = nodes[n] -> next;
    }
  }
  
  refresh();
}

int main() {
  
  plot_title = strdup("No data");
  
  for (char l = 0; l < 3; l++) {
    data[l] = create_list(DLL, free);
    axis_names[l] = strdup("No axis name");
  }
  
  axes_enabled = 0;

  // Setup ncurses
  initscr();
  use_default_colors();
  start_color();
  curs_set(0);
  cbreak();
  noecho();
  keypad(stdscr, TRUE);
  
  // Pair ncurses colors
  init_pair(1, -1, -1);
  init_pair(2, COLOR_GREEN,   -1);
  init_pair(3, COLOR_MAGENTA, -1);
  init_pair(4, COLOR_RED,     -1);
  init_pair(5, COLOR_YELLOW,  -1);
  init_pair(6, COLOR_BLUE,  -1);
  init_pair(7, COLOR_CYAN,  -1);
  init_pair(8, COLOR_BLACK, COLOR_YELLOW);
  
  erase();

  attron(COLOR_PAIR(1));
  mvaddstr(0, 0, "Sensor Overview");
  mvprintw(0, COLS / 2 - 16, "%s", plot_title);
  mvaddstr(LINES - 3,  1, "Legend"  );
  mvaddstr(LINES - 2,  4,   "0th:"  );
  mvaddstr(LINES - 2, 10, ", 1st:"  );
  mvaddstr(LINES - 2, 18, ", 2nd:"  );
  mvaddstr(LINES - 1,  4,   "0&1:"  );
  mvaddstr(LINES - 1, 10, ", 0&2:"  );
  mvaddstr(LINES - 1, 18, ", 1&2:"  );
  mvaddstr(LINES - 1, 26, ", all: *");
  
  //mvaddstr(LINES - 3, COLS - 39, "Controls");
  //mvaddstr(LINES - 2, COLS - 36, "c: cycle graph");
  //mvaddstr(LINES - 2, COLS - 18, "r: repair graph");
  //mvaddstr(LINES - 2, COLS - 1 "f: full experiment");
  mvaddstr(LINES - 1, COLS - 36, "b: back");
  
  mvaddstr(0, COLS - 6, "FEMTA"  );
  attroff(COLOR_PAIR(1));

  attron(COLOR_PAIR(4));
  mvaddstr(LINES - 2, 9, "*");
  attroff(COLOR_PAIR(4));

  attron(COLOR_PAIR(2));
  mvaddstr(LINES - 2, 17, "*");
  attroff(COLOR_PAIR(2));

  attron(COLOR_PAIR(6));
  mvaddstr(LINES - 2, 25, "*");
  attroff(COLOR_PAIR(6));

  attron(COLOR_PAIR(5));
  mvaddstr(LINES - 1, 9, "*");
  attroff(COLOR_PAIR(5));

  attron(COLOR_PAIR(3));
  mvaddstr(LINES - 1, 17, "*");
  attroff(COLOR_PAIR(3));
  
  attron(COLOR_PAIR(7));
  mvaddstr(LINES - 1, 25, "*");
  attroff(COLOR_PAIR(7));
  
  refresh();
  
  
  char input[1024];
  
  for (;;) {
    fgets(input, 1024, stdin);
    
    if (input[0] == 'd') {
      // d:value0 value1 value2
      
      double datum[3];
      
      if      (axes_enabled == 3) sscanf(input, "d:%lf %lf %lf", &datum[0], &datum[1], &datum[2]);
      else if (axes_enabled == 2) sscanf(input, "d:%lf %lf"    , &datum[0], &datum[1]           );
      else                   	  sscanf(input, "d:%lf"        , &datum[0]                      );      
      
      for (char l = 0; l < axes_enabled; l++) {
	// remember these data
	
	// must use malloc due to size of doubles exceeding void *
	double * value = malloc(sizeof(double));    
	
	value[0] = datum[l];
	
	list_insert(data[l], value);
	
	if (data[l] -> elements > COLS - 4) {
	  list_remove(data[l], data[l] -> tail);
	}
      }
      
      update_plot();
    }
    
    else if (input[0] == 't') {
      // t:title
      
      free(plot_title);
      plot_title = strdup(input + 2);
      
      // wipe plots
      axes_enabled = 0;
    }
    
    else if (input[0] == 'a') {
      // a:axis_name
      
      free(axis_names[axes_enabled]);
      
      axis_names[axes_enabled] = strdup(input + 4);
      
      axes_enabled++;
    }
    
    else if (input[0] == 'q')    // user is quitting
      break;
  }
  
  endwin();
}
