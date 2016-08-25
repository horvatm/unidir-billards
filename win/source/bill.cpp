#include "bill.h"

int main (int argc, char *argv[]) {
 
   /* init threads */
  g_thread_init (NULL);
  gdk_threads_init ();
  gdk_threads_enter ();

  gtk_init (&argc, &argv);
  
  GladeXML  *data_handler;

  data.pbill = 0;
  data.filename = 0;
  data.pixmap = 0;
  data.poin_manual_running = false;
  data.poin_random_running = false;
  data.stop_running = false;
  data.poin_manual_listening = false;
  data.poin_gc = 0;
  /* 
    load the interface 
  */
  data_handler = glade_xml_new ("bill.glade", NULL, NULL);

  /* 
    references to all important widgets 
  */

  data.main_window = glade_xml_get_widget (data_handler, "main_window");
  
  // configuration tab
  data.text_view = glade_xml_get_widget (data_handler, "text_view");
  data.saveas_button = glade_xml_get_widget (data_handler, "saveas_button");
  data.open_button = glade_xml_get_widget (data_handler, "open_button");
  data.shape_drawingarea = glade_xml_get_widget (data_handler, "shape_drawingarea");
  data.verify_button = glade_xml_get_widget (data_handler, "verify_button");
  data.m_spinbutton = glade_xml_get_widget (data_handler, "m_spinbutton");
    
  // trajectory tab
  data.plot_button = glade_xml_get_widget (data_handler, "plot_button");
  data.phi_spinbutton =glade_xml_get_widget (data_handler,"phi_spinbutton");
  data.l_spinbutton = glade_xml_get_widget (data_handler, "l_spinbutton");
  data.vx_spinbutton = glade_xml_get_widget (data_handler, "vx_spinbutton");
  data.vy_spinbutton = glade_xml_get_widget (data_handler, "vy_spinbutton");
  data.t_spinbutton = glade_xml_get_widget (data_handler, "t_spinbutton");
  data.traj_drawingarea = glade_xml_get_widget (data_handler, "traj_drawingarea");


  // poincare tab
  data.poin_start_button = glade_xml_get_widget (data_handler, "poin_start_button");
  data.poin_stop_button = glade_xml_get_widget (data_handler, "poin_stop_button");
  data.poin_clear_button = glade_xml_get_widget (data_handler, "poin_clear_button");
  data.poin_drawingarea = glade_xml_get_widget (data_handler, "poin_drawingarea");
  data.N_spinbutton = glade_xml_get_widget (data_handler, "N_spinbutton");
  data.n_spinbutton = glade_xml_get_widget (data_handler, "n_spinbutton");
  data.seed_spinbutton = glade_xml_get_widget (data_handler, "seed_spinbutton");
  data.phi_poin_spinbutton = glade_xml_get_widget (data_handler, "phi_poin_spinbutton");
  data.sign_spinbutton = glade_xml_get_widget (data_handler, "sign_spinbutton");
  data.manual_radiobutton = glade_xml_get_widget (data_handler, "manual_radiobutton");
  data.random_radiobutton = glade_xml_get_widget (data_handler, "random_radiobutton");
  data.sleep_hscale = glade_xml_get_widget (data_handler, "sleep_hscale");
  
  data.b[0] = glade_xml_get_widget (data_handler, "togglebutton1");
  data.b[1] = glade_xml_get_widget (data_handler, "togglebutton2");
  data.b[2] = glade_xml_get_widget (data_handler, "togglebutton3");
  data.b[3] = glade_xml_get_widget (data_handler, "togglebutton4");
  data.b[4] = glade_xml_get_widget (data_handler, "togglebutton5");
  data.b[5] = glade_xml_get_widget (data_handler, "togglebutton6");
  data.b[6] = glade_xml_get_widget (data_handler, "togglebutton7");
  data.b[7] = glade_xml_get_widget (data_handler, "togglebutton8");
  data.b[8] = glade_xml_get_widget (data_handler, "togglebutton9");
  data.b[9] = glade_xml_get_widget (data_handler, "togglebutton10");

  #if defined(WIN32)
  gtk_label_set_text(
    GTK_LABEL(glade_xml_get_widget (data_handler, "label_sleep")),
   "Sleep [miliseconds]");
  gtk_range_set_range(GTK_RANGE(data.sleep_hscale), 0, 100);
  gtk_range_set_value(GTK_RANGE(data.sleep_hscale), 0);
  #endif

  // general
  data.messages_textview = glade_xml_get_widget (data_handler, "messages_textview");
  data.message_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(data.messages_textview));

	g_object_unref (data_handler);

  /* 
    connect the signals in the interface 
  */
  g_signal_connect (G_OBJECT (data.main_window), 
                    "delete_event",  G_CALLBACK (gtk_main_quit), NULL);
	g_signal_connect (G_OBJECT (data.main_window), 
                    "destroy", G_CALLBACK (gtk_main_quit), NULL);

   // configuration tab
  g_signal_connect (G_OBJECT (data.open_button), 
                    "clicked", G_CALLBACK (open_button_clicked), NULL);
  g_signal_connect (G_OBJECT (data.saveas_button), 
                    "clicked", G_CALLBACK(saveas_button_clicked),NULL);  
  g_signal_connect (G_OBJECT (data.verify_button), 
                      "clicked", G_CALLBACK (verify_button_clicked), NULL);
  g_signal_connect (G_OBJECT (data.shape_drawingarea), 
                    "expose-event", G_CALLBACK (on_shape_expose_event), NULL);
  // trajectory tab
  g_signal_connect (G_OBJECT (data.plot_button), 
                      "clicked", G_CALLBACK (plot_button_clicked), NULL);
                      
  g_signal_connect (G_OBJECT (data.traj_drawingarea), 
                    "expose-event", G_CALLBACK (on_traj_expose_event), NULL);

  // poincare tab 
  g_signal_connect (G_OBJECT (data.poin_start_button), 
                      "clicked", G_CALLBACK (poin_start_button_clicked), NULL);
  g_signal_connect (G_OBJECT (data.poin_stop_button), 
                      "clicked", G_CALLBACK (poin_stop_button_clicked), NULL);
  g_signal_connect (G_OBJECT (data.poin_clear_button), 
                      "clicked", G_CALLBACK (poin_clear_button_clicked), NULL);
                    
  g_signal_connect (G_OBJECT (data.poin_drawingarea), 
                    "expose-event", G_CALLBACK (on_poin_expose_event), NULL);

  gtk_signal_connect (GTK_OBJECT (data.poin_drawingarea), "button-press-event",
		      G_CALLBACK (poin_button_press_event), NULL);

  gtk_signal_connect (GTK_OBJECT(data.poin_drawingarea),"configure-event",
		      G_CALLBACK(on_poin_configure_event), NULL);

  g_signal_connect (G_OBJECT(data.sleep_hscale), "value-changed",
           G_CALLBACK (on_sleep_value_changed), NULL);

  for (int i = 0; i < 10; ++i)
    g_signal_connect (G_OBJECT(data.b[i]), "pressed",  
      G_CALLBACK (on_color_toggled), NULL);

  /*
    Setting the color buttons - palette
 */
  // from rgb.txt for X11
  const char *color_str[10]= {"yellow", "orange", "pink", "purple","red",
                               "cyan","green", "magenta", "blue", "black"};

  for (int i = 0; i < 10; ++i) {
    gdk_color_parse (color_str[i], &data.c[i]);

    #if defined(WIN32)
    GdkPixmap* pixmap = gdk_pixmap_new(data.main_window->window, 10, 16, -1);       
   
    GdkGC  *gc = gdk_gc_new (pixmap); 
    gdk_gc_set_rgb_fg_color(gc, &data.c[i]);
    gdk_draw_rectangle (pixmap, gc, TRUE, 0, 0, 10, 16);
  
    GtkWidget *image = gtk_image_new_from_pixmap(pixmap, NULL); 
    gtk_button_set_image(GTK_BUTTON(data.b[i]), image);

    #else
    gtk_widget_modify_bg (data.b[i], GTK_STATE_NORMAL, &data.c[i]);
    gtk_widget_modify_bg (data.b[i], GTK_STATE_ACTIVE, &data.c[i]);
    gtk_widget_modify_bg (data.b[i], GTK_STATE_PRELIGHT, &data.c[i]);
    #endif
  }

  data.color_selected = 0;
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(data.b[data.color_selected]),TRUE);

	gtk_widget_show_all (data.main_window);


  /* 
    start the event loop 
  */
  gtk_main ();
  
  if (data.pbill) delete data.pbill;
  if (data.filename) delete [] data.filename;
  
  gdk_threads_leave ();

  return EXIT_SUCCESS;
}  
