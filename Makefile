CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -Ofast -Wno-volatile
TARGET_EXEC := raytracer.bin
BUILD_DIR := ./build
SRC_DIRS := ./src
INCLUDE_DIRS := ./include ./libs
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.c' -or -name '*.s')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)
INC_DIRS := $(shell find $(SRC_DIRS) -type d) $(shell find $(INCLUDE_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
# -MMD and -MP flags together generate Makefiles having .d instead of .o
CPPFLAGS := $(INC_FLAGS) -MMD -MP
LDFLAGS =

all: $(TARGET_EXEC)

# The final build step.
$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

## Build step for C source
#$(BUILD_DIR)/%.c.o: %.c
#	mkdir -p $(dir $@)
#	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Build step for C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -r $(BUILD_DIR) $(TARGET_EXEC)

# Include the .d makefiles. The - at the front suppresses the errors of missing
# Makefiles. Initially, all the .d files will be missing, and we don't want those
# errors to show up.
-include $(DEPS)

run: $(TARGET_EXEC)
	./$(TARGET_EXEC)
