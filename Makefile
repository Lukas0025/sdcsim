cmake: build
	cd build && \
	make

build:
	mkdir build
	cd build && \
	cmake ..

test: cmake
	cd build && \
	./unitTests

clean:
	rm -rf build
