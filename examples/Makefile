CFLAGS=-Wall -Wextra -std=c11 $(SANFLAGS)
CXXFLAGS=-Wall -Wextra $(SANFLAGS)

INCS=-I ../include
LIBS=-L ../src
LINK=-lalign -lpthread

all: nw_example sw_example nw_cpp sw_cpp

sw_example: sw_example.c ../src/libalign.a ../include/nucleics_align.h
	$(CC) $< -o $@ $(CFLAGS) $(INCS) $(LIBS) $(LINK)

nw_example: nw_example.c ../src/libalign.a ../include/nucleics_align.h
	$(CC) $< -o $@ $(CFLAGS) $(INCS) $(LIBS) $(LINK)

nw_cpp: nw_example.cpp ../src/libalign.a ../include/nucleics_align.h
	$(CXX) $< -o $@ $(CXXFLAGS) $(INCS) $(LIBS) $(LINK)

sw_cpp: sw_example.cpp ../src/libalign.a ../include/nucleics_align.h
	$(CXX) $< -o $@ $(CXXFLAGS) $(INCS) $(LIBS) $(LINK)

clean:
	rm -rf nw_example sw_example nw_cpp sw_cpp *.greg *.dSYM

.PHONY: all clean
