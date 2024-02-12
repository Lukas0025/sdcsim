cmake: build
	cd build && \
	make

build:
	mkdir build
	cd build && \
	cmake -D CMAKE_BUILD_TYPE=Release ..

debug_build:
	mkdir debug_build
	cd debug_build && \
	cmake -D CMAKE_BUILD_TYPE=Debug ..

debug: debug_build
	cd debug_build && \
	make

test: cmake
	cd build && \
	./unitTests

clean:
	rm -rf build
	rm -rf debug_build
