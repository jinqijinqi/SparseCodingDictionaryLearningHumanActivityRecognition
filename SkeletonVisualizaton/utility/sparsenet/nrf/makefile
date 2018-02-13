VCHOME = D:\applications\compilers\VC98
VCINCLUDE=$(VCHOME)\include 
INCLUDE = $(VCINCLUDE)
LIB = $(VCHOME)\lib

CC= cl
CFLAGS=	-Ox /G6 /c  -I$(VCINCLUDE)


MODULES=	frprmn.c linmin.c brent.c mnbrak.c nrutil.c

OBJECTS = $(MODULES:.c=.obj)



nrc: $(OBJECTS) 
	link /OUT:nrfopt.dll /DLL /MAP /IMPLIB:nrfopt.lib $(VCHOME)\lib\BINMODE.OBJ $(OBJECTS)
	cp nrfopt.* ..


