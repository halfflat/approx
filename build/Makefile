.PHONY: clean all realclean test

demos:=run_exp
tests:=

all:: $(demos)

top=..
sources:=$(wildcard $(top)/src/*.cc)
-include $(patsubst %.cc,obj/%.d,$(notdir $(sources)))

CXXFLAGS+=-std=c++11 -g -O0 -MMD -MP
CPPFLAGS+=-I$(top)/src

VPATH = $(top)/src

define cxx-compile
$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $<
endef

define cxx-link
$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)
endef

# build demos
obj/%.o: %.cc
	$(cxx-compile)

#run_log: obj/run_log.o
#	$(cxx-link)

run_exp: obj/run_exp.o obj/common_opt.o
	$(cxx-link)

# clean up

clean:
	rm -f obj/*.o

realclean: clean
	rm -f $(demos) obj/*.d

