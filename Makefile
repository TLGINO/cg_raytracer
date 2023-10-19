# vars
CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -Ofast -Wno-volatile
SRC_DIR = .
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_DIR = obj
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
OUTPUT_BIN = raytracer.bin
RESULT_PPM = result.ppm

# targets
.PHONY: all clean run


all: $(OUTPUT_BIN)

$(OUTPUT_BIN): $(OBJ_FILES)
	$(CXX) -o $@ $^

# change to SRC_FILES
$(OBJ_DIR)/%.o: $(SRC_FILES) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

run: all
	./$(OUTPUT_BIN)
	eog $(RESULT_PPM) &

clean:
	rm -rf $(OUTPUT_BIN) $(RESULT_PPM) $(OBJ_DIR)
