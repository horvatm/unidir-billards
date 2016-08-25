#define __bill_h

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <ios>
#include <unistd.h>

// GUI related headers
#include <gtk/gtk.h>
#include <glade/glade.h>
#include <pthread.h>

// library for computing everythig
#include "calc.h"
#include "MersenneTwister.h"

/*

  Structure containing widgets of the GUI
*/

struct Tdata {
 
  GtkWidget *main_window,
            // configuration tab
            *text_view, 
            *open_button, 
            *saveas_button, 
            *verify_button, 
            *shape_drawingarea,
            *m_spinbutton,
            
            // trajectory tab
            *plot_button,
            *phi_spinbutton,
            *l_spinbutton,
            *vx_spinbutton,
            *vy_spinbutton,
            *t_spinbutton,
            *traj_drawingarea,
            
            // poincare tab
            *poin_start_button,
            *poin_stop_button,
            *poin_clear_button,
            *poin_drawingarea,
            *N_spinbutton,
            *n_spinbutton,
            *seed_spinbutton,
            *phi_poin_spinbutton,
            *sign_spinbutton,
            *manual_radiobutton,
            *random_radiobutton,
            *sleep_hscale,
            *b[10],

            // General
            *messages_textview;


  // colors for the poincare tab
  int color_selected;

  GdkColor  c[10];  

  GdkGC *poin_gc;

  // text buffer of the message  text_view
  GtkTextBuffer *message_buffer;

  // billiard calc class  
  Tbill *pbill;
  
  char *filename;
  
  // plotting the shape of the billiard
  double xmin,xmax,ymin, ymax, Lx, Ly;
  std::vector<Tpoint> curve[2];

  // potting trajectory
  double cross_section_phi;
  std::vector<Tpoint> traj;

  // plotting  poincare surface of section  
  GdkPixmap *pixmap;

  // threads of manual and random point selection calculations
  pthread_t poin_manual, poin_random;
  
  // flags for proper working of the threads
  bool poin_manual_running, 
       poin_manual_listening,
       poin_random_running, 
       stop_running;

} data;

/*
  Converting any variable into a string
    input: t  
    returns: string
*/
template <class T> inline std::string to_string (const T& t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

/*
  For reporting messages to the message_text_view in the main_window
*/
void message(const std::string & s){

  GtkTextIter end_start_iter; 
  gtk_text_buffer_get_end_iter (data.message_buffer, &end_start_iter);

  gtk_text_buffer_insert(data.message_buffer, &end_start_iter, s.c_str(), -1); 

  /* get end iter again */
  gtk_text_buffer_get_end_iter (data.message_buffer, &end_start_iter);
 
  /* get the current ( cursor )mark name */
  GtkTextMark *insert_mark = gtk_text_buffer_get_insert (data.message_buffer);
 
  /* move mark and selection bound to the end */
  gtk_text_buffer_place_cursor(data.message_buffer, &end_start_iter);

  /* scroll to the end view */   
  gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW (data.messages_textview),
                insert_mark, 0.0, TRUE, 0.0, 1.0); 
}


/**********************************************************
   Configuration tab
**********************************************************/

void open_button_clicked (GtkWidget *widget, gpointer user_data) {
  
  
  GtkWidget *dialog = gtk_file_chooser_dialog_new ("Open File",
				      GTK_WINDOW(data.main_window),
				      GTK_FILE_CHOOSER_ACTION_OPEN,
				      GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
				      GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
				      NULL);

  if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT) {

    if (data.filename) delete [] data.filename;
    
    data.filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    
    message("Reading filename" + to_string(data.filename)+ "\n");
    std::ifstream is(data.filename);
    
    is.seekg (0, std::ios::end);
    int len = is.tellg();
    is.seekg (0, std::ios::beg);

    char *text = new char [len];
    
    is.read (text, len);
    is.close();
    
    while (len > 0 && *(text + len -1) < 32) --len;

    gchar *gtext =
     g_convert(text, len, "UTF-8","ISO8859-1", NULL, NULL, NULL);
    
    GtkTextBuffer *text_view_buffer 
    = gtk_text_view_get_buffer(GTK_TEXT_VIEW(data.text_view));

    gtk_text_buffer_set_text (text_view_buffer, gtext, -1);
    
    delete [] text;
  }

  gtk_widget_destroy (dialog);
}


void saveas_button_clicked (GtkWidget *widget, gpointer user_data) {  


  GtkWidget *dialog = gtk_file_chooser_dialog_new ("Save File",
				        GTK_WINDOW(data.main_window),
				        GTK_FILE_CHOOSER_ACTION_SAVE,
				        GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
				        GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
				        NULL);
				        
  gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (dialog), TRUE);

  if (data.filename)
    gtk_file_chooser_set_filename (GTK_FILE_CHOOSER (dialog), data.filename);
  else {  
    gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER (dialog), "./");
    gtk_file_chooser_set_current_name (GTK_FILE_CHOOSER (dialog), "untitled.cfg");
  }

  if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT) {

    if (data.filename) delete [] data.filename;
    
    data.filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    
    message("SaveAs filename" + to_string(data.filename)+ "\n");
    GtkTextBuffer *text_view_buffer = 
    gtk_text_view_get_buffer(GTK_TEXT_VIEW(data.text_view));
    
    GtkTextIter start, end;
    
    gtk_text_buffer_get_start_iter (text_view_buffer, &start);
    gtk_text_buffer_get_end_iter (text_view_buffer, &end);

    gchar *text = 
    gtk_text_buffer_get_text (text_view_buffer, &start, &end, TRUE);

    std::ofstream os(data.filename);    
    os.write (text, strlen(text));
    
    delete [] text;
    
  }

  gtk_widget_destroy (dialog);
}

void verify_button_clicked (GtkWidget *widget, gpointer user_data) { 
  
  // get the configuration  
  GtkTextBuffer *text_view_buffer = 
  gtk_text_view_get_buffer(GTK_TEXT_VIEW(data.text_view));

  if (text_view_buffer == 0){
    message("No text buffer!\n");
    return;
  }

  GtkTextIter start, end;
    
  gtk_text_buffer_get_start_iter (text_view_buffer, &start);
  gtk_text_buffer_get_end_iter (text_view_buffer, &end);
  
  gchar *gtext = 
  gtk_text_buffer_get_text (text_view_buffer, &start, &end, FALSE);
  
  if (gtext == 0) {
    message("No text!\n");
    return;  
  }

  if (strlen(gtext) > 0) {

    // getting readable text in ISO8859-1
    std::stringstream ss (
      g_convert(gtext, strlen(gtext), "ISO8859-1", "UTF-8", 0, 0, 0), 
      std:: stringstream::in | std::stringstream::out);

    // reading and reporting on the formula
    std::string s;

    double d;  
    ss >> d;
    
    message("Billiard width: d=" + to_string(d) + "\n");

    Tmode a;
    std::vector <Tmode> coefs;
    
    while (ss >> a.s >> a.n >> a.val) {
      coefs.push_back(a);

      if (a.n != 0){

        if (a.val > 0 && s.length() > 0) s += "+";

        s += to_string(a.val);

        if (a.s) s += "sin("; else s += "cos(";

        if (a.n != 1) s += to_string(a.n) + "x)"; else s += "x)";
  
      } else if (a.s == 0) {
      
        if (a.val > 0 && s.length() > 0) s += "+";

        s += to_string(a.val);
      }
    }
    message("Inner boundary: r(x)=" + s + "\n");
    
    if (data.pbill) delete data.pbill;
      
    data.pbill = new Tbill(d, coefs);
    
    message("The billiard calculation are configured! Ready to calculate!\n");
  } else {
    message("Text is zero length!\n");
    delete [] gtext;
    return;
  }
 
  delete [] gtext;

  /* produce the curve */
  
  // getting number of peaces
  int m = int(gtk_spin_button_get_value(GTK_SPIN_BUTTON(data.m_spinbutton)));

  message("Plotting shape: m=" + to_string(m) + "\n");

  if (m > 0) {
    
    if (data.curve[0].size()) { 
      for (int i =0; i < 2; ++i) data.curve[i].clear(); 
      data.traj.clear();
    }
             
    // calculating the shape of the billiard      
    data.pbill->plot_shape(m, data.curve);

    data.xmax = data.ymax = -1e300;
    data.xmin = data.ymin = 1e300;

    std::vector<Tpoint>::iterator it;

    for (int i = 0; i < 2;++i)
      for (it = data.curve[i].begin(); it != data.curve[i].end(); ++it){
        data.xmax = std::max(data.xmax,(*it).x);
        data.xmin = std::min(data.xmin,(*it).x);

        data.ymax = std::max(data.ymax,(*it).y);
        data.ymin = std::min(data.ymin,(*it).y);
      }

  
    data.Lx = data.xmax-data.xmin+1e-10;
    data.Ly = data.ymax-data.ymin+1e-10;

    /* redraw the cairo canvas completely by exposing it */

    GdkRegion *region
     = gdk_drawable_get_clip_region (data.shape_drawingarea->window);
     
    gdk_window_invalidate_region (data.shape_drawingarea->window, region, TRUE);
    gdk_window_process_updates (data.shape_drawingarea->window, TRUE);
      
    gdk_region_destroy (region);
  }
}


gboolean on_shape_expose_event (GtkWidget *widget, 
                                GdkEventExpose *event, gpointer user_data){

  cairo_t *cr = gdk_cairo_create (widget->window);



  int scrwidth  = widget->allocation.width,
      scrheight = widget->allocation.height;
  
  cairo_set_source_rgb (cr, 0, 0, 0);
  cairo_set_line_width (cr, 1);
	cairo_rectangle (cr, 0, 0, scrwidth, scrheight);
  cairo_stroke (cr);

  if (data.curve[0].size()) {
    
    cairo_set_line_width (cr, 2.0);
    cairo_set_source_rgba (cr, 1, 0, 0, 0.6);  

    std::vector<Tpoint>::iterator it;
    
    for (int i = 0; i < 2;++i){
      it = data.curve[i].begin();

      cairo_move_to (cr, 
        scrwidth*(0.9*((*it).x - data.xmin)/data.Lx+0.05),
        scrheight*(1-(0.9*((*it).y - data.ymin)/data.Ly+0.05))
      );

      while (++it != data.curve[i].end()) {
        cairo_line_to (cr, 
          scrwidth*(0.9*((*it).x - data.xmin)/data.Lx+0.05),
          scrheight*(1-(0.9*((*it).y - data.ymin)/data.Ly+0.05))
        );
      }
      cairo_close_path(cr);
    }
    cairo_stroke (cr);
  }
  cairo_destroy (cr);

  return FALSE;
}

/**********************************************************
   Trajectory tab
**********************************************************/


void plot_button_clicked (GtkWidget *widget, gpointer user_data) {

  if (data.pbill) {
    // getting initial point
    double d = data.pbill->get_width(), 
           phi = gtk_spin_button_get_value(GTK_SPIN_BUTTON(data.phi_spinbutton)),
           l = d*gtk_spin_button_get_value(GTK_SPIN_BUTTON(data.l_spinbutton)),
           v[2]= {gtk_spin_button_get_value(GTK_SPIN_BUTTON(data.vx_spinbutton)),
                  gtk_spin_button_get_value(GTK_SPIN_BUTTON(data.vy_spinbutton))},
           t_end = gtk_spin_button_get_value(GTK_SPIN_BUTTON(data.t_spinbutton));

    data.cross_section_phi = phi;
    data.traj.clear();
           
    message("Cros-section (phi,l)=(" + 
            to_string(phi) + ","+ to_string(l)+")\n");

    message("Initial speed (v_x,v_y)=(" + 
            to_string(v[0]) + ","+ to_string(v[1])+")\n");

    message("Time t=" + to_string(t_end)+"\n");

    data.pbill->plot_traj(phi, l, v, t_end, &data.traj);

     /* redraw the cairo canvas completely by exposing it */
    message("Traj.size=" + to_string(data.traj.size())+"\n");

    GdkRegion *region
     = gdk_drawable_get_clip_region (data.traj_drawingarea->window);
     
    gdk_window_invalidate_region (data.traj_drawingarea->window, region, TRUE);
    gdk_window_process_updates (data.traj_drawingarea->window, TRUE);
      
    gdk_region_destroy (region);
  } else message("Warning:First configure the billiard\n");
}

#define EPS 1e-10
gboolean on_traj_expose_event (GtkWidget *widget, 
                                GdkEventExpose *event, gpointer user_data){
  // read data

  int scrwidth  = widget->allocation.width,
      scrheight = widget->allocation.height;

  cairo_t *cr = gdk_cairo_create (widget->window);
  
  cairo_set_source_rgb (cr, 0, 0, 0);
  cairo_set_line_width (cr, 1);
	cairo_rectangle (cr, 0, 0, scrwidth, scrheight);
  cairo_stroke (cr);

  if (data.curve && data.traj.size()) {
    double Lx = data.xmax-data.xmin+EPS, Ly = data.ymax-data.ymin+EPS;
    cairo_set_line_width (cr, 2.0);
    cairo_set_source_rgba (cr, 1, 0, 0, 0.6);

    // plotting billiard's shape
    std::vector<Tpoint>::iterator it;

    for (int i = 0; i < 2;++i){
      it = data.curve[i].begin();

      cairo_move_to (cr, 
        scrwidth*(0.9*((*it).x - data.xmin)/data.Lx+0.05),
        scrheight*(1-(0.9*((*it).y - data.ymin)/data.Ly+0.05))
      );

      while (++it != data.curve[i].end()) {
        cairo_line_to (cr, 
          scrwidth*(0.9*((*it).x - data.xmin)/data.Lx + 0.05),
          scrheight*(1-(0.9*((*it).y - data.ymin)/data.Ly + 0.05))
        );
      }
      cairo_close_path(cr);
    }
    cairo_stroke (cr);

    // plotting trajectory
    it = data.traj.begin();
    
    cairo_set_line_width (cr, 1.0);
    cairo_set_source_rgb (cr, 0, 0, 0);
    
    cairo_move_to (cr, 
      scrwidth*(0.9*((*it).x - data.xmin)/data.Lx+0.05),
      scrheight*(1-(0.9*((*it).y - data.ymin)/data.Ly+0.05))
    );
    
    while (++it != data.traj.end())
      cairo_line_to (cr, 
        scrwidth*(0.9*((*it).x - data.xmin)/data.Lx+0.05),
        scrheight*(1-(0.9*((*it).y - data.ymin)/data.Ly+0.05))
      );
    cairo_stroke (cr);
  
    // drawing the line of cross-section
    cairo_set_line_width (cr, 2.0);
    cairo_set_source_rgb (cr, 0, 0, 256);

    double x, y, nx, ny, 
           r = data.pbill->rad(data.cross_section_phi);

    sincos(data.cross_section_phi, &y, &x);

    x *= r;
    y *= r;

    cairo_move_to (cr,
        scrwidth*(0.9*(x - data.xmin)/data.Lx + 0.05),
        scrheight*(1-(0.9*(y- data.ymin)/data.Ly + 0.05))
    );

    data.pbill->normal(data.cross_section_phi, nx, ny);

    x += (nx *= data.pbill->get_width());
    y += (ny *= data.pbill->get_width());

    cairo_line_to (cr,
        scrwidth*(0.9*(x - data.xmin)/data.Lx + 0.05),
        scrheight*(1-(0.9*(y- data.ymin)/data.Ly + 0.05))
    );
    
    cairo_stroke (cr);
  }  

  cairo_destroy (cr);

  return FALSE;
}

/**********************************************************
    Poincare tab
**********************************************************/
/*
  Draw a point x,y in [0,1] rescaled on the canvas. 

  NOTE: The origin is top left, but mathematical origin should be bottom left.
*/
void draw_poin(const double & x, const double &y){
  GdkRectangle update_rect;
  
  update_rect.x = int(data.poin_drawingarea->allocation.width*x-1);
  update_rect.y = int(data.poin_drawingarea->allocation.height*y-1);

  update_rect.width = 2;
  update_rect.height = 2;

 // Draw background 
  //std::cerr << data.color_selected << '\n';
  //gdk_gc_set_foreground(data.poin_gc, &data.c[data.color_selected]); 

	gdk_gc_set_foreground(data.poin_gc, &data.c[data.color_selected]);

  gdk_draw_rectangle (data.pixmap,
          data.poin_gc,
          //data.poin_drawingarea->style->black_gc,
		      TRUE,
		      update_rect.x, update_rect.y,
		      update_rect.width, update_rect.height);

  gtk_widget_draw (data.poin_drawingarea, &update_rect);

}


struct Tpoin_random {
  double phi;
  int N, n, sign;
  unsigned long seed;
} p_rnd;


void *poin_random_simulate(void *args){
  int i, j;

  double t, d = data.pbill->get_width(), J0[2], J1[2];
  
  MTRand rnd(p_rnd.seed);

  #define BORDER 1e-4  
  for (i = 0; i < p_rnd.N && !data.stop_running; ++i) {

    J0[0] = (d - 2*BORDER)*rnd.rand53() + BORDER;
    J0[1] = 2*(1 - BORDER)*rnd.rand53() - 1 + BORDER;
    
    gdk_threads_enter ();
    draw_poin(J0[0]/d, 0.5*(1-J0[1]));   //1 - 0.5*(J0[1]+1)
    gdk_threads_leave ();

    for (j = 0; j < p_rnd.n && !data.stop_running; ++j){

    //  t = bill.poincare_map(sign, phi, J0, J1, &std::cerr);

      t = data.pbill->poincare_map(p_rnd.sign, p_rnd.phi, J0, J1);

      if (data.pbill->status != 0){  
        gdk_threads_enter ();
        switch (data.pbill->status){
          case 1: 
            message("Poincare:There was an error in INSIDE-TO-WALL\n");
          break;
          case 2: 
             message("Poincare:There was an error in WALL-TO-WALL\n");
          break;
          case 3: 
             message("Poincare:Length along Out-of-bounds " + to_string(J1[0])+ "\n"); 
          break;
          case 4: 
            message("Poincare:Time is negative \n");
          break;
        }
        message("Poincare:Point=("+to_string(J0[0])+","+to_string(J0[1])+")\n");
        gdk_threads_leave ();

        break;    
      } else {
        J0[0] = J1[0]; J0[1] = J1[1];
        gdk_threads_enter ();
        draw_poin(J0[0]/d, 0.5*(1-J0[1]));//1 - 0.5*(J0[1]+1)
        gdk_threads_leave ();
      }
    }
  }
  #undef BORDER 
  gdk_threads_enter ();
  message("End of the simulation: RANDOM-POINTS \n");
  gdk_threads_leave ();

  data.poin_random_running = false;
 
  return NULL;
}

G_LOCK_DEFINE_STATIC (p_man);

struct Tpoin_manual{
  bool change;
  double phi, J[2];
  int sign, sleep;
  
} p_man;


void *poin_manual_simulate(void *args){

  double t, d = data.pbill->get_width(), J0[2], J1[2];
  
  while (!data.stop_running) {

    if (p_man.change){    
      G_LOCK (p_man);
      J0[0] = p_man.J[0];
      J0[1] = p_man.J[1];
      p_man.change = false;
      gdk_threads_enter ();
      message("Changing the point: ("+
              to_string(J0[0]) + "," + to_string(J0[1]) + ")\n");
      gdk_threads_leave ();
      G_UNLOCK (p_man);
    }

    gdk_threads_enter ();
    draw_poin(J0[0]/d, 0.5*(1-J0[1]));
    gdk_threads_leave ();

    //  t = bill.poincare_map(sign, phi, J0, J1, &std::cerr);

    t = data.pbill->poincare_map(p_man.sign, p_man.phi, J0, J1);

    if (data.pbill->status != 0){  
      gdk_threads_enter ();
      switch (data.pbill->status){
        case 1: 
          message("Poincare:There was an error in INSIDE-TO-WALL\n");
        break;
        case 2: 
           message("Poincare:There was an error in WALL-TO-WALL\n");
        break;
        case 3: 
           message("Poincare:Length along Out-of-bounds " + to_string(J1[0])+ "\n"); 
        break;
        case 4: 
          message("Poincare:Time is negative \n");
        break;
      }
      message("Poincare:Point=("+to_string(J0[0])+","+to_string(J0[1])+")\n");
      gdk_threads_leave ();

      break;    
    } else {
      J0[0] = J1[0]; 
      J0[1] = J1[1];

      gdk_threads_enter ();
      draw_poin(J0[0]/d, 0.5*(1-J0[1]));
      gdk_threads_leave ();
    }

    #if defined(WIN32)
    if (p_man.sleep) Sleep(p_man.sleep);
    #else
    usleep(p_man.sleep);
    #endif
  }
  #undef BORDER 
  gdk_threads_enter ();
  message("End of the simulation: MANUAL-POINTS \n");
  gdk_threads_leave ();

  data.poin_manual_running = false;
  data.poin_manual_listening = false;

  return NULL;
}

double get_value(GtkWidget *w) {
 return gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
}

void poin_start_button_clicked (GtkWidget *widget, gpointer user_data){

  if (data.poin_manual_running || data.poin_random_running) {
    message("Something is already running.");
    message("Before starting another setup STOP the previous!\n");
    return;
  }

  if (data.pbill==0) {
    message("Warning:First configure the billiard!\n");
    return;
  }

  if (GTK_TOGGLE_BUTTON (data.manual_radiobutton)->active){

    p_man.phi = get_value(data.phi_poin_spinbutton);
    p_man.sign = (int)get_value(data.sign_spinbutton);
    p_man.sleep = (int)gtk_range_get_value(GTK_RANGE(data.sleep_hscale));
    p_man.change = false;

    message("Listening to the MOUSE CURSOR for the MANUAL-POINTS simulation\n");
    message(" with parameters:\n");
    message("phi=" + to_string(p_man.phi)+"\n");
    message("sign=" + to_string(p_man.sign)+"\n");
    message("sleep=" + to_string(p_man.sleep)+"\n");
    message("Click on the canvas to start simulation\n");

    data.stop_running = false;
    data.poin_manual_listening = true;
    data.poin_manual_running = false;

  } else {
    p_rnd.phi =get_value(data.phi_poin_spinbutton);
    p_rnd.N = (int)get_value(data.N_spinbutton);
    p_rnd.n = (int)get_value(data.n_spinbutton);
    p_rnd.sign = (int)get_value(data.sign_spinbutton);
    p_rnd.seed = (int)get_value(data.seed_spinbutton);

    message("Running the RANDOM-POINTS simulation with parameters:\n");
    message("phi=" + to_string(p_rnd.phi)+"\n");
    message("sign=" + to_string(p_rnd.sign)+"\n");
    message("seed=" + to_string(p_rnd.seed)+"\n");
    message("N=" + to_string(p_rnd.N)+"\n");
    message("n=" + to_string(p_rnd.n)+"\n");

    data.stop_running = false;
    data.poin_random_running = true;

    pthread_create (&data.poin_random, NULL, poin_random_simulate, NULL);

   /*
    GdkRegion *region
     = gdk_drawable_get_clip_region (data.traj_drawingarea->window);
     
    gdk_window_invalidate_region (data.traj_drawingarea->window, region, TRUE);
    gdk_window_process_updates (data.traj_drawingarea->window, TRUE);
      
    gdk_region_destroy (region);
   */
  }
}

void on_sleep_value_changed (GtkWidget *widget, gpointer user_data){

 if (data.poin_manual_running)  { 
   G_LOCK(p_man);
   p_man.sleep = (int) gtk_range_get_value(GTK_RANGE(data.sleep_hscale));
   G_UNLOCK(p_man);
   message("Sleep time changed, sleep="+to_string(p_man.sleep)+ "\n");
 }

}

void poin_stop_button_clicked (GtkWidget *widget, gpointer user_data){
 data.stop_running = true; 
}

void poin_clear_button_clicked (GtkWidget *widget, gpointer user_data){
  GdkRectangle update_rect;
  
  update_rect.x = 0;
  update_rect.y = 0;

  update_rect.width = data.poin_drawingarea->allocation.width;
  update_rect.height = data.poin_drawingarea->allocation.height;

  gdk_draw_rectangle (data.pixmap,
		      data.poin_drawingarea->style->white_gc,
		      TRUE,
		      update_rect.x, update_rect.y,
		      update_rect.width, update_rect.height);

  gtk_widget_draw (data.poin_drawingarea, &update_rect);  
}

gboolean on_poin_expose_event (GtkWidget *widget, 
                                GdkEventExpose *event, gpointer user_data) {
  gdk_draw_pixmap(widget->window,
		  widget->style->fg_gc[GTK_WIDGET_STATE (widget)],
		  data.pixmap,
		  event->area.x, event->area.y,
		  event->area.x, event->area.y,
		  event->area.width, event->area.height);

  return FALSE;
}

/* Create a new backing pixmap of the appropriate size */
gboolean  on_poin_configure_event(GtkWidget *widget,GdkEventConfigure *event){

  if (data.poin_gc == 0) {
    GdkColormap *colormap = gtk_widget_get_colormap(widget);
    data.poin_gc = gdk_gc_new(widget->window);

    for (int i = 0; i < 10; ++i) 
      gdk_colormap_alloc_color(colormap, &data.c[i], FALSE, TRUE);      
	}

	
  if (data.pixmap)  gdk_pixmap_unref(data.pixmap);

  data.pixmap = gdk_pixmap_new(widget->window,
			  widget->allocation.width,
			  widget->allocation.height,
			  -1);

  gdk_draw_rectangle (data.pixmap,
		      widget->style->white_gc,
		      TRUE,
		      0, 0,
		      widget->allocation.width,
		      widget->allocation.height);

  return TRUE;
}

gboolean  poin_button_press_event( GtkWidget *widget, GdkEventButton *event ){

  if (data.pbill &&
      data.poin_manual_listening && 
      event->button == 1 && 
      data.pixmap != NULL){

    G_LOCK(p_man); 
    p_man.change = true;
    p_man.J[0] = data.pbill->get_width()*double(event->x)/widget->allocation.width;
    p_man.J[1] = 2*(1-double(event->y)/widget->allocation.height)-1;
    G_UNLOCK(p_man); 

    if (!data.poin_manual_running) {
      data.poin_manual_running = true;
      pthread_create (&data.poin_manual, NULL, poin_manual_simulate, NULL);
    }
  }

  return TRUE;
}

void on_color_toggled (GtkWidget *widget, gpointer user_data){

  int tmp = atoi(gtk_widget_get_name(widget)+12)-1;
 
  if (data.color_selected != tmp){
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(data.b[data.color_selected]), FALSE);  
    data.color_selected = tmp;
  }
 // gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), TRUE);  
 // std::cerr << tmp << "\n";
}

#undef __bill_h
