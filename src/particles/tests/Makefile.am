check_PROGRAMS = test-pbp

AM_CPPFLAGS= -Wall -DFSF_MAINTAINER=\"@ILP_MAINTAINER@\" -std=gnu++11
AM_LDFLAGS =  -Wall -pedantic -ansi

test_pbp_SOURCES = testPBP.cc 
test_pbp_CPPFLAGS = $(AM_CPPFLAGS) \
        $(DUNE_CFLAGS) \
        -I $(top_srcdir)

