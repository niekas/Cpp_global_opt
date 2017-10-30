do: run

compile:
	@echo "==================================================================================="
	g++ gkls.c rnd_gen.c main.cpp -o libre.out -I/usr/include/python2.7 -lpython2.7

run: compile
	clear
	./libre.out --func_name=zdt1 --alpha=0.1 --max_calls=15000
	
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
	valgrind --tool=memcheck --leak-check=full --track-origins=yes -v ./libre.out --func_name=zdt1 --alpha=0.1

profiler:  compile
	valgrind --tool=callgrind ./libre.out --func_name=zdt1 --alpha=0.1
	# git clone https://github.com/jrfonseca/gprof2dot bin/gprof2dot
	./bin/gprof2dot/gprof2dot.py -f callgrind callgrind.out.* | dot -Tsvg -o profile.svg
