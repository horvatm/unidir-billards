AC_INIT(bill.cpp, calc.cpp)

AM_INIT_AUTOMAKE(bill,0.1)
AM_INIT_AUTOMAKE(calc,0.1)

AC_CONFIG_HEADERS(config.h)

PKG_CHECK_MODULES(LIBGLADE, libglade-2.0 >= 2.0.0, true,
   AC_MSG_ERROR([you need libglade >= 2.0.0 - http://www.jamesh.id.au/software/libglade ]))
AC_SUBST(LIBGLADE_CFLAGS)
AC_SUBST(LIBGLADE_LIBS)

PKG_CHECK_MODULES(GTK, gtk+-2.0 >= 2.6, true, 
  AC_MSG_ERROR([you need gtk+ >= 2.6 - Please see http://www.gtk.org ]))
AC_SUBST(GTK_CFLAGS)
AC_SUBST(GTK_LIBS)

PKG_CHECK_MODULES(GTHREADS, gthread-2.0)
AC_SUBST(GTHREADS_CFLAGS)
AC_SUBST(GTHREADS_LIBS)

AC_PROG_CC

AC_PROG_CXX

AC_PROG_INSTALL

AC_OUTPUT(Makefile)
