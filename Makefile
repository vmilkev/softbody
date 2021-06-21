shell = /bin/sh

all:
	cd graphics && $(MAKE)
	cd soil_compression && $(MAKE)
	cd regular_cluster && $(MAKE)

clean:
	cd graphics && $(MAKE) clean
	cd soil_compression && $(MAKE) clean
	cd regular_cluster && $(MAKE) clean

