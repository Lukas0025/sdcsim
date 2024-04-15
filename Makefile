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
	@echo "Starting E2E"
	for f in ./test/e2e/*; do \
		if [ -d "$$f" -a $$(echo -n "$$f" | tail -c 1) != "-" ]; then \
			mkdir ./test-submission/ && \
			cp -r ./build/* ./test-submission/ && \
			cp -r $$f/*     ./test-submission/ && \
			cd ./test-submission && $(MAKE) -s test && \
			cd ../ && \
			rm -rf ./test-submission/; \
		fi \
	done
	@echo "[info] All tests PASS \e[33m ^_^ \e[39m"

plots: cmake
	cd plots && \
	./make.sh

clean:
	rm -rf build
	rm -rf debug_build
	rm -rf ./test-submission/
