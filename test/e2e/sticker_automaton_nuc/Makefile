NAME="Sticker automaton nucleotides level"

test:
	for it in $$(seq 1 100); do \
		./sdcsim in.sdcasm -f assembly -a --silent -n -T 25 > ./test_out 2>&1 && \
		diff ./expected ./test_out || continue && \
		echo "[info] test $(NAME) PASS in $$it iteration \e[33m ^_^ \e[39m" && \
		exit 0 ; \
	done && \
	echo "[info] test $(NAME) FAIL \e[31m -_- \e[39m" && \
	exit 1