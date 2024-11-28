MAIN_MAKEFILE_DIR = LIB
MAIN_MAKEFILE = $(MAIN_MAKEFILE_DIR)/libmakefile
include $(MAIN_MAKEFILE)

EXEC = PROG
.INTERMEDIATE: $(EXEC)
.DEFAULT_GOAL := $(EXEC)

SRC_PATH ?= SRC
SRC_SOURCES = $(wildcard $(SRC_PATH)/*.c)
SRC_OBJECTS = $(SRC_SOURCES:.c=.o)
.INTERMEDIATE: $(SRC_OBJECTS)

ifeq ($(wildcard $(GSL_LIB_PATH)/libgsl.a),)
    $(info )
    $(info =============================================)
    $(info ==> Warning: GSL_LIB/lib not found)
    $(info --> Searching for system-wide GSL installation paths...)
    $(info --> Tested paths:)
    $(info ..... 1. $(GSL_LIB_PATH) (not found)) 
    $(info ..... 2. /usr/local/lib (default system path))
    $(info ..... 3. /usr/lib (default system path))
    $(info =============================================)
    GSL_LIB_FLAG =
    LIBRARY_LOCATION := $(shell \
        if echo 'int main() { return 0; }' > gsl_test.c && \
           $(CC) gsl_test.c -o gsl_test -L/usr/local/lib -lgsl -lgslcblas -lm >/dev/null 2>&1; then \
           echo "/usr/local/lib"; \
        elif $(CC) gsl_test.c -o gsl_test -L/usr/lib -lgsl -lgslcblas -lm >/dev/null 2>&1; then \
           echo "/usr/lib"; \
        else \
           echo "Error: No GSL library found"; \
        fi; \
        rm -f gsl_test gsl_test.c)
    ifneq ($(findstring Error, $(LIBRARY_LOCATION)),)
        $(error $(LIBRARY_LOCATION))
    else
        $(info --> System-wide GSL library found at $(LIBRARY_LOCATION).)
    endif
else
    $(info )
    $(info =============================================)
    $(info ==> Using GSL libraries from $(GSL_LIB_PATH))
    $(info --> Found GSL at $(GSL_LIB_PATH).)
    $(info =============================================)
    GSL_LIB_FLAG = -L$(GSL_LIB_PATH)
    LIBRARY_LOCATION := $(GSL_LIB_PATH)
endif

$(SRC_PATH)/%.o: $(SRC_PATH)/%.c
	@echo "..... --> Compiling $< with libraries from $(LIBRARY_LOCATION)."
	@$(CC) $(CFLAGS) $(INCLUDE_ALL_HEADERS) -c $< -o $@

.PHONY: $(EXEC)

UNAME_S := $(shell uname)
ifeq ($(UNAME_S), Darwin)
    TIME_CMD = @/usr/bin/time -p
else
    TIME_CMD = 
endif

$(EXEC): $(SRC_OBJECTS) $(PERSONNAL_MERGEDLIB_FILE)
	@echo "============================================="
	@echo "..... ==> Linking $@..."
	$(CC) $(CFLAGS) $(INCLUDE_ALL_HEADERS) $^ -o $@ $(INCLUDE_ALL_LIBS)
	@echo "..... ==> Build complete!"
	@echo "..... ==> Running $@..."
	@echo "============================================="
	@$(TIME_CMD) ./$(EXEC)

.PHONY: clean
clean:
	$(RM) $(SRC_OBJECTS) $(EXEC) $(PERSONNAL_MERGEDLIB_FILE)
	
.PHONY: help
help:
	@echo "Available commands:"
	@echo "  make           - Build the executable and run it."
	@echo "  make clean     - Remove all generated files."
	@echo "  make help      - Show this help message."