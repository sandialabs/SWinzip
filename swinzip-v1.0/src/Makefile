SHELL = /bin/sh

include CURRENT_defines


SUBDIRS = Utils Alpert Sampling StOMP
export
	
all:
	for i in $(SUBDIRS) ;\
	  do cd $$i; $(MAKE); cd .. ;\
	 done ;

clean:
	rm -f ../lib/*.*;\
	rm -f ../include/*.*;\
	for i in $(SUBDIRS) ;\
	 do rm -f $$i/*.a $$i/*.o ;\
	done ;

