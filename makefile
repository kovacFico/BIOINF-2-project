CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
INCLUDES = -Iinclude -Isrc
BIN = bin

# ===============================
# Source files
# ===============================

PREPROCESS_SRC = \
	./apps/preprocess.cpp \
	./preprocesing/genome_preprocesing.cpp

HMM_INIT_SRC = \
	./apps/hmm_params_init.cpp \
	./hmm/hmm.cpp \
	./hmm/hmm_io.cpp

TRAIN_SRC = \
	./apps/train.cpp \
	./hmm/hmm.cpp \
	./hmm/hmm_io.cpp \
	./algorithms/baum_welch.cpp \
	./algorithms/forward_backward.cpp \
	./train_functions/train_func.cpp

DECODE_SRC = \
	./apps/decode_and_evaluation.cpp \
	./hmm/hmm.cpp \
	./hmm/hmm_io.cpp \
	./algorithms/forward_backward.cpp \
	./postprocesing/decoded_postprocesing.cpp \
	./evaluation/evaluation.cpp

LAUNCHER_SRC = ./main.cpp

# ===============================
# Targets
# ===============================

all: dirs preprocess hmm_init train decode launcher

dirs:
	mkdir -p $(BIN)

preprocess:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(PREPROCESS_SRC) -o $(BIN)/preprocess

hmm_init:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(HMM_INIT_SRC) -o $(BIN)/hmm_init

train:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TRAIN_SRC) -o $(BIN)/train

decode:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DECODE_SRC) -o $(BIN)/decode_and_evaluation

launcher:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LAUNCHER_SRC) -o $(BIN)/launcher

clean:
	rm -rf $(BIN)
