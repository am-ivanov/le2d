CC = nvcc
CFLAGS = #--ptxas-options=--verbose

SMS = 20 30 35 37 50

$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))

HIGHEST_SM := $(lastword $(sort $(SMS)))
GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)

all: le_core.cu main.cu le_core.h
	$(CC) $(CFLAGS) $(GENCODE_FLAGS) le_core.cu main.cu -DUSE_DOUBLE -O2 -o cuda_le2d_d -lm
	$(CC) $(CFLAGS) $(GENCODE_FLAGS) le_core.cu main.cu -O2 -o cuda_le2d_f -lm
	$(CC) $(CFLAGS) $(GENCODE_FLAGS) -use_fast_math le_core.cu main.cu -O2 -o cuda_le2d_f_fm -lm

debug: le_core.cu main.cu le_core.h
	$(CC) $(CFLAGS) le_core.cu main.cu -g -DDEBUG -o cuda_le2d_debug_f -lm
	$(CC) $(CFLAGS) le_core.cu main.cu -DUSE_DOUBLE -g -DDEBUG -o cuda_le2d_debug_d -lm

clean: 
	rm cuda_le2d*
