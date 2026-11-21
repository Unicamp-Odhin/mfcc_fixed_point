CC = ccache gcc
AR = ar

# Flags padrão
CFLAGS = -O3 -march=native
CPPFLAGS = -Iinclude -Isrc
LDFLAGS =
LDLIBS = -lm

# Diretórios
BUILD_DIR = build
SRC_DIR = src
INCLUDE_DIR = include

# Cores
BLUE := \033[1;34m
GREEN := \033[1;32m
YELLOW := \033[1;33m
RED := \033[1;31m
NC := \033[0m

# Flags adicionais (valores padrão)
DEBUG ?= 0
CREATE_DATABANK ?= 1
LOG ?= 1
VERBOSE ?= 1

# Debug tem prioridade máxima
ifeq ($(DEBUG),1)
    CFLAGS = -g -O0
else
    CFLAGS = -O3 -march=native
endif

# Adicionar flags condicionais
ifeq ($(CREATE_DATABANK),1)
    CFLAGS += -DCONFIG_CREATE_DATABANK
endif

ifeq ($(LOG),1)
    CFLAGS += -DCONFIG_LOG
endif

ifeq ($(VERBOSE),1)
    CFLAGS += -DCONFIG_VERBOSE
endif

# Mostrar flags finais
$(info [INFO] CFLAGS: $(CFLAGS))

# Lista de fontes da biblioteca
LIB_SOURCES = \
    $(SRC_DIR)/mfcc.c \
    $(SRC_DIR)/dct.c \
    $(SRC_DIR)/fft_fp.c \
    $(SRC_DIR)/mel.c \
    $(SRC_DIR)/process.c

LIB_OBJECTS = \
    $(LIB_SOURCES:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

LIB_NAME = libmfcc.a

all: buildFolder $(BUILD_DIR)/$(LIB_NAME)

.PHONY: all clean buildFolder lib

lib: buildFolder $(BUILD_DIR)/$(LIB_NAME)

# Criar biblioteca estática
$(BUILD_DIR)/$(LIB_NAME): $(LIB_OBJECTS)
	@printf "$(GREEN)[AR]$(NC) Criando biblioteca estática\n"
	$(AR) rcs $@ $^

# Compilação
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@printf "$(BLUE)[CC]$(NC) Compilando $<\n"
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

# Criar diretório
buildFolder:
	@mkdir -p $(BUILD_DIR)

clean:
	@printf "$(RED)[CLEAN]$(NC) Removendo arquivos gerados\n"
	rm -rf $(BUILD_DIR)
	rm -rf "$(REF_C_DIR)/dumps/"
	rm -rf "$(TESTS_DIR)/ref_vectors/"

# Paralelização
MAKEFLAGS += -j$(shell nproc)
