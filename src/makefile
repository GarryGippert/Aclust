BIN = ../bin
LFLAGS = -lm
EXE = aclust

.SUFFIXES:
.SUFFIXES: .c .o
 
.c:
	$(CC) $? -o $@ $(LIB) $(LFLAGS)
 
.c.o:   $<
	$(CC) -c $<

all: $(EXE)

clean:
	/bin/rm -rf *.o *~ core* $(EXE)

install: $(EXE)
	@-files="$(EXE)"; \
	for f in $$files; do \
	  target=$(BIN)/`basename $$f`; \
	  install -m 775 $$f $$target; \
	  /bin/rm -f $$f; \
	  echo $$f installed at $$target and removed locally; \
	done
