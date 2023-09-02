CC = clang++
LFLAGS = -std=c++11 -Wall

all: LFLAGS += -Ofast -march=native
all: cmix enwik9-preproc

debug: LFLAGS += -ggdb
debug: cmix enwik9-preproc

cmix: src/coder/decoder.cpp src/coder/decoder.h src/coder/encoder.cpp src/coder/encoder.h src/context-manager.cpp src/context-manager.h src/contexts/bit-context.cpp src/contexts/bit-context.h src/contexts/bracket-context.cpp src/contexts/bracket-context.h src/contexts/combined-context.cpp src/contexts/combined-context.h src/contexts/context-hash.cpp src/contexts/context-hash.h src/contexts/context.h src/contexts/indirect-hash.cpp src/contexts/indirect-hash.h src/contexts/interval-hash.cpp src/contexts/interval-hash.h src/contexts/interval.cpp src/contexts/interval.h src/contexts/sparse.cpp src/contexts/sparse.h src/mixer/byte-mixer.cpp src/mixer/byte-mixer.h src/mixer/lstm-layer.cpp src/mixer/lstm-layer.h src/mixer/lstm.cpp src/mixer/lstm.h src/mixer/mixer-input.cpp src/mixer/mixer-input.h src/mixer/mixer.cpp src/mixer/mixer.h src/mixer/sigmoid.cpp src/mixer/sigmoid.h src/mixer/sse.cpp src/mixer/sse.h src/models/bracket.cpp src/models/bracket.h src/models/byte-model.cpp src/models/byte-model.h src/models/direct-hash.cpp src/models/direct-hash.h src/models/direct.cpp src/models/direct.h src/models/indirect.cpp src/models/indirect.h src/models/fxcmv1.cpp src/models/fxcmv1.h src/models/match.cpp src/models/match.h src/models/model.h src/models/paq8.cpp src/models/paq8.h src/models/ppmd.cpp src/models/ppmd.h src/predictor.cpp src/predictor.h src/preprocess/dictionary.cpp src/preprocess/dictionary.h src/preprocess/preprocessor.cpp src/preprocess/preprocessor.h src/runner.cpp src/states/nonstationary.cpp src/states/nonstationary.h src/states/run-map.cpp src/states/run-map.h src/states/state.h
	$(CC) $(LFLAGS) src/coder/decoder.cpp src/coder/encoder.cpp src/context-manager.cpp src/contexts/bit-context.cpp src/contexts/bracket-context.cpp src/contexts/combined-context.cpp src/contexts/context-hash.cpp src/contexts/indirect-hash.cpp src/contexts/interval-hash.cpp src/contexts/interval.cpp src/contexts/sparse.cpp src/mixer/byte-mixer.cpp src/mixer/lstm-layer.cpp src/mixer/lstm.cpp src/mixer/mixer-input.cpp src/mixer/mixer.cpp src/mixer/sigmoid.cpp src/mixer/sse.cpp src/models/bracket.cpp src/models/byte-model.cpp src/models/direct-hash.cpp src/models/direct.cpp src/models/indirect.cpp src/models/fxcmv1.cpp src/models/match.cpp src/models/paq8.cpp src/models/ppmd.cpp src/predictor.cpp src/preprocess/dictionary.cpp src/preprocess/preprocessor.cpp src/runner.cpp src/states/nonstationary.cpp src/states/run-map.cpp -o cmix

enwik9-preproc: src/enwik9-preproc/article_reorder.h src/enwik9-preproc/main.cpp src/enwik9-preproc/misc.h src/enwik9-preproc/phda9_preprocess.h
	$(CC) $(LFLAGS) src/enwik9-preproc/main.cpp -o enwik9-preproc

clean:
	rm -f cmix enwik9-preproc
