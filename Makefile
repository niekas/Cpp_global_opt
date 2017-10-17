do: run

compile:
	@echo "==================================================================================="
	g++ gkls.c rnd_gen.c main.cpp -o libre.out

run: compile
	clear
	./libre.out --func_name=ep1 --alpha=0.1
	
test:
	g++ test.cpp -o test.out
	clear
	./test.out

cat:
	cat *.sh.po*
	cat *.sh.o*

queue:
	qstat -f
	# qdel

num:
	ls results/Disimpl-v/ | wc -l

mem_check:  compile
	valgrind --tool=memcheck --leak-check=full --track-origins=yes -v ./libre.out --func_cls=1 --func_id=1

profiler:  compile
	valgrind --tool=callgrind ./libre.out --func_cls=1 --func_id=1
	# git clone https://github.com/jrfonseca/gprof2dot bin/gprof2dot
	./bin/gprof2dot/gprof2dot.py -f callgrind callgrind.out.* | dot -Tsvg -o profile.svg
