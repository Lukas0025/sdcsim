NAME="Sticker automaton"

test:
	./sdcsim in.sdcasm -f assembly -a --silent > ./test_out 2>&1
	diff ./expected ./test_out || (echo "[info] test $(NAME) FAIL \e[31m -_- \e[39m" && exit 1)
	@echo "[info] test $(NAME) PASS \e[33m ^_^ \e[39m"