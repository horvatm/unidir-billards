rem Setting variables  (compiler and paths)

set comp="C:\Dev-Cpp\bin\g++.exe"
set gtk=..\GTK
set pthreads=..\Pre-built.2

rem Compiling

%comp% "bill.cpp" -o "bill.exe" -O3 -I"%gtk%\INCLUDE"  -I"%gtk%\INCLUDE\GTK-2.0"  -I"%gtk%\INCLUDE\GLIB-2.0"  -I"%gtk%\INCLUDE\PANGO-1.0"  -I"%gtk%\INCLUDE\CAIRO"  -I"%gtk%\INCLUDE\ATK-1.0"  -I"%gtk%\INCLUDE\GTKGLEXT-1.0"  -I"%gtk%\LIB\GTK-2.0\INCLUDE"  -I"%gtk%\LIB\GLIB-2.0\INCLUDE"  -I"%gtk%\LIB\GTKGLEXT-1.0\INCLUDE"  -I"%gtk%\INCLUDE\LIBGLADE-2.0"  -I"%gtk%\INCLUDE\LIBXML2"  -I"%pthreads%\include" -L"%pthreads%\lib" -L"%gtk%\lib" -mms-bitfields -lgtk-win32-2.0 -lgdk-win32-2.0 -lgthread-2.0 -lgdi32 -lole32 -luuid -latk-1.0  -lgdk_pixbuf-2.0 -lpangowin32-1.0 -lgdi32 -lpango-1.0 -lgobject-2.0  -lgmodule-2.0 -lglib-2.0 -lintl -liconv -latk-1.0 -lglade-2.0 -lpthreadGCE2 -lcairo



rem Compiling ./calc

%comp% "calc.cpp" -o "calc.exe" -O3
