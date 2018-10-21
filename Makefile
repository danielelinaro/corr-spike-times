
CFLAGS := -I/usr/local/include
CPPFLAGS := -I/usr/local/include
LDFLAGS := -L/usr/local/lib -l gsl -l gslcblas

generate_correlated_spike_times: erfinv.o correlatedspikes.o Correlations.o
	g++ $(LDFLAGS) -o $@ $^

%.o: %.c
	gcc $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	g++ $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o

