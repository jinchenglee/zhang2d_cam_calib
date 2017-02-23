# Source file list
SET (SRC_LIST cminpack.c m3.c main.c matrix.c minimizing.c zhang2dcalib.c)

# Add an executable binary file
ADD_EXECUTABLE( main ${NCURSES_LIBRARIES} ${SRC_LIST}) 

# Link math and ncurses library
TARGET_LINK_LIBRARIES (main m ncurses)

