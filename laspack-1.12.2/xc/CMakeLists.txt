include_directories(..)

set(xc_sources
    getopts.c
    xstring.c
)

add_library(xc ${xc_sources})
