# --- SYSTEM ---

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
BASISDIR    = /opt/ibm/ILOG

# --- DIRECTORIES ---

CCC = g++ -std=gnu++11 -no-pie -Iincludes
BASISILOG  = $(shell find $(BASISDIR) -maxdepth 1 -type d -name "CPLEX_Studio*" | sort -V | tail -1)
CONCERTDIR = $(BASISILOG)/concert
CPLEXDIR   = $(BASISILOG)/cplex

# --- FLAGS ---

CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -DIL_STD -Wno-deprecated-declarations -Wno-ignored-attributes
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

# --- OPTIMIZATION FLAGS ---

DEBUG_OPT = -DNDEBUG -O3
#DEBUG_OPT = -g3 -O0
#PROF = -pg
PROF =

CFLAGS += $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I./include $(DEBUG_OPT) -c $(PROF)

LDFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -ldl

# ---- COMPILE  ----
SRC_DIR_cnv   := src/src_cnv
OBJ_DIR_cnv   := obj/obj_cnv

SRC_DIRS_cnv  := $(shell find $(SRC_DIR_cnv) -type d)
OBJ_DIRS_cnv  := $(addprefix $(OBJ_DIR_cnv)/,$(SRC_DIRS_cnv))

SOURCES_cnv   := $(shell find $(SRC_DIR_cnv) -name '*.cpp')
OBJ_FILES_cnv := $(addprefix $(OBJ_DIR_cnv)/, $(SOURCES_cnv:.cpp=.o))

SRC_DIR_aff   := src/src_aff
OBJ_DIR_aff   := obj/obj_aff

SRC_DIRS_aff  := $(shell find $(SRC_DIR_aff) -type d)
OBJ_DIRS_aff  := $(addprefix $(OBJ_DIR_aff)/,$(SRC_DIRS_aff))

SOURCES_aff   := $(shell find $(SRC_DIR_aff) -name '*.cpp')
OBJ_FILES_aff := $(addprefix $(OBJ_DIR_aff)/, $(SOURCES_aff:.cpp=.o))

vpath %.cpp $(SRC_DIRS_cnv)
vpath %.cpp $(SRC_DIRS_aff)

# ---- TARGETS ----

EXECUTABLE1 = msa_cnv 
EXECUTABLE2 = msa_aff

EXECUTABLES = $(EXECUTABLE1) $(EXECUTABLE2)

all: $(EXECUTABLES)

$(EXECUTABLE1): makedir $(SOURCES_cnv) $(OBJ_FILES_cnv) 
	$(CCC) $(OBJ_FILES_cnv) $(LDFLAGS) $(PROF) -o $@

$(EXECUTABLE2): makedir $(SOURCES_aff) $(OBJ_FILES_aff) 
	$(CCC) $(OBJ_FILES_aff) $(LDFLAGS) $(PROF) -o $@

$(OBJ_DIR_cnv)/%.o: %.cpp
	$(CCC) $(CFLAGS) $< -o $@
$(OBJ_DIR_aff)/%.o: %.cpp
	$(CCC) $(CFLAGS) $< -o $@

makedir: $(OBJ_DIRS_cnv)
makedir: $(OBJ_DIRS_aff)

$(OBJ_DIRS_cnv):
	@mkdir -p $@
$(OBJ_DIRS_aff):
	@mkdir -p $@

clean:
	@rm -rf $(OBJ_DIR_cnv)
	@rm -rf $(OBJ_DIR_aff)
	@rm -rf $(EXECUTABLE1)
	@rm -rf $(EXECUTABLE2)


