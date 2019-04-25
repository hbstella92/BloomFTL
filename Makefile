CC = gcc

LIBS = -lm\
	   #-ljemalloc\

CFLAGS = \
		-O3\
        -march=armv8-a+crypto\
		#-mcpu=cortex-a53\
		#-g\
#-DPF\

OBJS = $(patsubst %.c,%.o,$(SRCS))
SRCS = main.c bloomfilter.c \
        sha256-arm.c\

TARGET = simulationFTL

all : $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

clean :
	rm -rf $(OBJS) $(TARGET)
