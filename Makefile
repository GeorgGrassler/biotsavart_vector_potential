BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all ninja clean
all: ninja

$(BUILD_NINJA):
	cmake -G Ninja -B $(BUILD_DIR)

ninja: $(BUILD_NINJA)
	cd $(BUILD_DIR) && ninja

clean:
	rm -rf $(BUILD_DIR)
