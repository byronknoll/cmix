CC = g++
CFLAGS = -std=c++11 -Wall -c
LFLAGS = -std=c++11 -Wall

OBJS = build/preprocessor.o build/encoder.o build/decoder.o build/predictor.o build/logistic.o build/mixer-input.o build/mixer.o build/sse.o build/manager.o build/direct.o build/direct-hash.o build/indirect.o build/nonstationary.o build/run-map.o build/byte-run.o build/match.o build/dmc.o build/ppm.o build/paq8l.o build/paq8hp.o build/context-hash.o build/sparse.o build/indirect-hash.o build/bit-context.o

all: CFLAGS += -Ofast -march=native -s
all: LFLAGS += -Ofast -march=native -s
all: cmix

debug: CFLAGS += -ggdb
debug: LFLAGS += -ggdb
debug: cmix

cmix: $(OBJS) src/runner.cpp
	$(CC) $(LFLAGS) $(OBJS) src/runner.cpp -o cmix

build/preprocessor.o: src/preprocess/preprocessor.h src/preprocess/preprocessor.cpp src/preprocess/textfilter.cpp src/predictor.h
	$(CC) $(CFLAGS) src/preprocess/preprocessor.cpp -o build/preprocessor.o

build/encoder.o: src/coder/encoder.h src/coder/encoder.cpp src/predictor.h
	$(CC) $(CFLAGS) src/coder/encoder.cpp -o build/encoder.o

build/decoder.o: src/coder/decoder.h src/coder/decoder.cpp src/predictor.h
	$(CC) $(CFLAGS) src/coder/decoder.cpp -o build/decoder.o

build/predictor.o: src/predictor.h src/predictor.cpp src/mixer/mixer-input.h src/mixer/mixer.h src/sse.h src/models/model.h src/models/direct.h src/models/direct-hash.h src/models/indirect.h src/models/byte-run.h src/models/match.h src/models/dmc.h src/models/ppm.h src/models/facade.h src/models/paq8l.h src/models/paq8hp.h src/manager.h src/contexts/context-hash.h src/contexts/sparse.h src/contexts/indirect-hash.h src/contexts/bit-context.h src/models/facade.h src/mixer/logistic.h
	$(CC) $(CFLAGS) src/predictor.cpp -o build/predictor.o

build/logistic.o: src/mixer/logistic.h src/mixer/logistic.cpp
	$(CC) $(CFLAGS) src/mixer/logistic.cpp -o build/logistic.o

build/mixer-input.o: src/mixer/mixer-input.h src/mixer/mixer-input.cpp src/mixer/logistic.h
	$(CC) $(CFLAGS) src/mixer/mixer-input.cpp -o build/mixer-input.o

build/mixer.o: src/mixer/mixer.h src/mixer/mixer.cpp src/mixer/logistic.h
	$(CC) $(CFLAGS) src/mixer/mixer.cpp -o build/mixer.o

build/sse.o: src/sse.h src/sse.cpp
	$(CC) $(CFLAGS) src/sse.cpp -o build/sse.o

build/manager.o: src/manager.h src/manager.cpp src/contexts/context.h src/contexts/bit-context.h
	$(CC) $(CFLAGS) src/manager.cpp -o build/manager.o

build/direct.o: src/models/direct.h src/models/direct.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/direct.cpp -o build/direct.o

build/direct-hash.o: src/models/direct-hash.h src/models/direct-hash.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/direct-hash.cpp -o build/direct-hash.o

build/indirect.o: src/models/indirect.h src/models/indirect.cpp src/states/state.h src/models/model.h
	$(CC) $(CFLAGS) src/models/indirect.cpp -o build/indirect.o

build/byte-run.o: src/models/byte-run.h src/models/byte-run.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/byte-run.cpp -o build/byte-run.o

build/match.o: src/models/match.h src/models/match.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/match.cpp -o build/match.o

build/dmc.o: src/models/dmc.h src/models/dmc.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/dmc.cpp -o build/dmc.o

build/ppm.o: src/models/ppm.h src/models/ppm.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/ppm.cpp -o build/ppm.o

build/paq8l.o: src/models/paq8l.h src/models/paq8l.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/paq8l.cpp -o build/paq8l.o

build/paq8hp.o: src/models/paq8hp.h src/models/paq8hp.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/paq8hp.cpp -o build/paq8hp.o

build/paq8pxd.o: src/models/paq8pxd.h src/models/paq8pxd.cpp src/models/model.h
	$(CC) $(CFLAGS) src/models/paq8pxd.cpp -o build/paq8pxd.o

build/nonstationary.o: src/states/nonstationary.h src/states/nonstationary.cpp src/states/state.h
	$(CC) $(CFLAGS) src/states/nonstationary.cpp -o build/nonstationary.o

build/run-map.o: src/states/run-map.h src/states/run-map.cpp src/states/state.h
	$(CC) $(CFLAGS) src/states/run-map.cpp -o build/run-map.o

build/context-hash.o: src/contexts/context-hash.h src/contexts/context-hash.cpp src/contexts/context.h
	$(CC) $(CFLAGS) src/contexts/context-hash.cpp -o build/context-hash.o

build/indirect-hash.o: src/contexts/indirect-hash.h src/contexts/indirect-hash.cpp src/contexts/context.h
	$(CC) $(CFLAGS) src/contexts/indirect-hash.cpp -o build/indirect-hash.o

build/sparse.o: src/contexts/sparse.h src/contexts/sparse.cpp src/contexts/context.h
	$(CC) $(CFLAGS) src/contexts/sparse.cpp -o build/sparse.o

build/bit-context.o: src/contexts/bit-context.h src/contexts/bit-context.cpp
	$(CC) $(CFLAGS) src/contexts/bit-context.cpp -o build/bit-context.o

clean:
	rm -f -r build/* cmix
